# 2021-01-07 15:08:58

Let me try the usual DGE analysis without cropping the number of genes by HUGO.

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]

data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))
data$Diagnosis = factor(data$Diagnosis, levels=c('Control', 'Case'))
data$substance_group = factor(data$substance_group)
data$comorbid_group = factor(data$comorbid_group_update)
data$evidence_level = factor(data$evidence_level)

library(caret)
pp_order = c('zv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)

# I'll need some annotations here to keep only the autosomes
# just because I know Derek annotated based on hg38
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/gencode.v36.annotation.sqlite')
```

This is still not covering everything. After further inspection of what was
missing, it looks like Derek used (GRCh38.p12) for the columns of the count
matrix, and we're currently in (GRCh38.p13) in the ensemble website. So, some of
the genes can only be found in the archive for p12. 

I then downloaded files from:

http://jul2019.archive.ensembl.org/info/data/ftp/index.html

Specifically, the GTF matching the p12 version, and created a DB with it:

```r
library(GenomicFeatures)
txdb <- makeTxDbFromGFF('~/Downloads/Homo_sapiens.GRCh38.97.gtf.gz')
saveDb(txdb, '~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
```

And that covers everything!

```r
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
noversion = data.frame(GENEID = substr(rownames(geneCounts), 1, 15))
noversion = merge(noversion, txdf, by='GENEID', sort=F)
imautosome = which(noversion$TXCHROM != 'X' &
                   noversion$TXCHROM != 'Y' &
                   noversion$TXCHROM != 'MT')
geneCounts = geneCounts[imautosome, ]
noversion = noversion[imautosome, ]
# there are no duplicates here for ACC

library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
geneCountsExpr = geneCounts[isexpr,]
genesExpr = noversion[isexpr,]

genes = DGEList( geneCountsExpr, genes=genesExpr ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
data.pm = cbind(data, mydata)
rownames(data.pm) = data$hbcc_brain_id
```

Now, let's see first if our main PM results are still there:

```r
num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data.pm[, x], data.pm[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(data.pm[, y], data.pm[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}

print(which(num_pvals < .01, arr.ind = T))
print(which(categ_pvals < .01, arr.ind = T))
```

```
                        row col
clusters                  2   1
RINe                      4   1
PMI                       5   1
RINe                      4   2
PMI                       5   2
clusters                  2   4
pcnt_optical_duplicates   1   7
clusters                  2   7
Age                       3   8
Age                       3   9
C4                        9   9
                row col
batch             1   1
batch             1   2
batch             1   7
MoD               3   9
substance_group   4   9
```

We had a few changes in the PCs: 1, 2, 4, 7, 8, 9.

```r
form = ~ Diagnosis + PC1 + PC2 + PC4 + PC7 + PC8 + PC9
design = model.matrix( form, data.pm)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

Let's do a quick run of GSEA:

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'
eval(parse(text=sprintf('res = rnaseq_%s', region)))

ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(geneid=res$GENEID, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

# my own GMTs
db = sprintf('my_%s_sets', region)
cat(region, db, '\n')
db_file = sprintf('~/data/post_mortem/%s.gmt', db)
enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                    organism="hsapiens",
                    enrichDatabaseFile=db_file,
                    enrichDatabaseType="genesymbol",
                    interestGene=tmp2,
                    outputDirectory = data_dir,
                    interestGeneType="ensembl_gene_id",
                    sigMethod="top", topThr=150000,
                    minNum=3,
                    isOutput=F, isParallel=T,
                    nThreads=ncpu, perNum=10000, maxNum=800))
out_fname = sprintf('%s/WG_ENSID_%s_%s_10K.csv', data_dir, region, db)
write.csv(enrichResult, file=out_fname, row.names=F)

DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5,
                                isOutput=F, isParallel=T,
                                nThreads=ncpu, perNum=10000)
    out_fname = sprintf('%s/WG_ENSID_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

Results are there, but FDR for GO is not as significant. Maybe if I remove some
genes using nz? Maybe the overlap with the fusion results will remain?

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]
data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))
data$Diagnosis = factor(data$Diagnosis, levels=c('Control', 'Case'))
data$substance_group = factor(data$substance_group)
data$comorbid_group = factor(data$comorbid_group_update)
data$evidence_level = factor(data$evidence_level)
library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
noversion = data.frame(GENEID = substr(rownames(geneCounts), 1, 15))
noversion = merge(noversion, txdf, by='GENEID', sort=F)
imautosome = which(noversion$TXCHROM != 'X' &
                   noversion$TXCHROM != 'Y' &
                   noversion$TXCHROM != 'MT')
geneCounts = geneCounts[imautosome, ]
noversion = noversion[imautosome, ]
# there are no duplicates here for ACC
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
geneCountsExpr = geneCounts[isexpr,]
genesExpr = noversion[isexpr,]
genes = DGEList( geneCountsExpr, genes=genesExpr ) 
genes = calcNormFactors( genes)
lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)
library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
data.pm = cbind(data, mydata)
rownames(data.pm) = data$hbcc_brain_id

num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data.pm[, x], data.pm[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}
categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(data.pm[, y], data.pm[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}

print(which(num_pvals < .01, arr.ind = T))
print(which(categ_pvals < .01, arr.ind = T))
```

```
                        row col
clusters                  2   1
RINe                      4   1
PMI                       5   1
RINe                      4   2
PMI                       5   2
clusters                  2   4
pcnt_optical_duplicates   1   7
clusters                  2   7
Age                       3   8
Age                       3   9
C4                        9   9
                row col
batch             1   1
batch             1   2
batch             1   7
MoD               3   9
substance_group   4   9
```

No changes for the PCs, even though I have a few more genes than before (22203
vs 22197). I checked and the code is working, but the filterByExp part throws
out different genes, and we end up with more.

```r
form = ~ Diagnosis + PC1 + PC2 + PC4 + PC7 + PC8 + PC9
design = model.matrix( form, data.pm)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

And we redo GSEA:

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'
eval(parse(text=sprintf('res = rnaseq_%s', region)))

ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(geneid=res$GENEID, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

# my own GMTs
db = sprintf('my_%s_sets', region)
cat(region, db, '\n')
db_file = sprintf('~/data/post_mortem/%s.gmt', db)
enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                    organism="hsapiens",
                    enrichDatabaseFile=db_file,
                    enrichDatabaseType="genesymbol",
                    interestGene=tmp2,
                    outputDirectory = data_dir,
                    interestGeneType="ensembl_gene_id",
                    sigMethod="top", topThr=150000,
                    minNum=3,
                    isOutput=F, isParallel=T,
                    nThreads=ncpu, perNum=10000, maxNum=800))
out_fname = sprintf('%s/WG_ENSIDnzv_%s_%s_10K.csv', data_dir, region, db)
write.csv(enrichResult, file=out_fname, row.names=F)

DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5,
                                isOutput=F, isParallel=T,
                                nThreads=ncpu, perNum=10000)
    out_fname = sprintf('%s/WG_ENSIDnzv_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

No, using nzv makes me keep only the GABA result... not sure if it's worth it.

Well, none of this is worth it if there's no relationship with the TWAS results.
Let's make a few TWAS cut-offs just to play with the data a bit:

```r
# bw
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_ACC.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)

wei = read.delim('~/data/expression_impute/fusion_twas-master/WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.P01.pos')
cnt = 75
wei$geneid = substring(wei$WGT, cnt, cnt+14)
mf = merge(fusion, wei, by='ID', all.x=T, all.y=F, sort=F)
load('~/tmp/ensid2.rdata')
mf$adjPval = p.adjust(mf$TWAS.P, method='fdr')

library(ActivePathways)
gmt = read.GMT('~/data/post_mortem/acc_developmental.gmt')
junk = gmt[1:6]
g = mf[mf$adjPval < .05, 'geneid']
name = 'FDR 0.05'
junk[[1]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$adjPval < .1, 'geneid']
name = 'FDR 0.1'
junk[[2]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .05/nrow(mf), 'geneid']
name = 'Bonferroni'
junk[[3]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .001, 'geneid']
name = 'P < .001'
junk[[4]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .005, 'geneid']
name = 'P < .005'
junk[[5]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .01, 'geneid']
name = 'P < .01'
junk[[6]] = list(id = name, genes = unique(g), name = name)
gmt_name = '~/data/post_mortem/FUSION_ACC_sets.gmt'
write.GMT(junk, gmt_name)
```

And let's see if GSEA works here:

```r
# bw
library(WebGestaltR)
data_dir = '~/data/rnaseq_derek/'
ncpu=6
region='acc'
eval(parse(text=sprintf('res = rnaseq_%s', region)))
ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(geneid=res$GENEID, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]
db_file = '~/data/post_mortem/FUSION_ACC_sets.gmt'
enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                    organism="hsapiens",
                    enrichDatabaseFile=db_file,
                    enrichDatabaseType="ensembl_gene_id",
                    interestGene=tmp2,
                    outputDirectory = data_dir,
                    interestGeneType="ensembl_gene_id",
                    sigMethod="top", topThr=150000,
                    minNum=3,
                    isOutput=F, isParallel=T,
                    nThreads=ncpu, perNum=10000, maxNum=800))
out_fname = sprintf('%s/WG_ENSIDnzv_FUSION_%s_10K.csv', data_dir, region)
write.csv(enrichResult, file=out_fname, row.names=F)
```

Nope, nothing even nominally... well, it was a valiant effort.

Should we try overrepresentation again?

```r
library(GeneOverlap)
both_res = merge(rnaseq_acc, mf, by.x='GENEID', by.y='geneid',
                 all.x=F, all.y=F)
thresh = c(.05, .01, .005, .001)
imp_nums = vector(length=length(thresh), mode='numeric')
rna_nums = vector(length=length(thresh), mode='numeric')
pvals = matrix(data=NA, nrow=length(thresh), ncol=length(thresh))
rownames(pvals) = sapply(thresh, function(x) sprintf('PM_%.3f', x))
colnames(pvals) = sapply(thresh, function(x) sprintf('IMP_%.3f', x))
inter_nums = pvals
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$TWAS.P < thresh[ti], 'GENEID']
    imp_nums[ti] = length(imp_genes)
    for (tr in 1:length(thresh)) {
        rna_genes = both_res[both_res$P.Value < thresh[tr], 'GENEID']
        rna_nums[tr] = length(rna_genes)
        go.obj <- newGeneOverlap(imp_genes, rna_genes,
                                 genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(imp_genes, rna_genes)
        pval = getPval(go.obj)
        pvals[tr, ti] = pval
        inter_nums[tr, ti] = length(intersect(rna_genes, imp_genes))
    }
}
print(pvals)
print(inter_nums)
```

```
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050 0.1016351 0.0665071 0.1886333 0.1907347
PM_0.010 0.3747058 0.1910814 0.6381977 0.3264215
PM_0.005 0.7689247 1.0000000 1.0000000 1.0000000
PM_0.001 1.0000000 1.0000000 1.0000000 1.0000000
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050        25        10         6         3
PM_0.010         6         3         1         1
PM_0.005         2         0         0         0
PM_0.001         0         0         0         0
```

Still nothing.

## Caudate

Let's repeat the same step for the Caudate now:

```r
myregion = 'Caudate'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]

data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))
data$Diagnosis = factor(data$Diagnosis, levels=c('Control', 'Case'))
data$substance_group = factor(data$substance_group)
data$comorbid_group = factor(data$comorbid_group_update)
data$evidence_level = factor(data$evidence_level)

library(caret)
pp_order = c('zv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
noversion = data.frame(GENEID = substr(rownames(geneCounts), 1, 15))
noversion = merge(noversion, txdf, by='GENEID', sort=F)
imautosome = which(noversion$TXCHROM != 'X' &
                   noversion$TXCHROM != 'Y' &
                   noversion$TXCHROM != 'MT')
geneCounts = geneCounts[imautosome, ]
noversion = noversion[imautosome, ]
# there are no duplicates here for Caudate

library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
geneCountsExpr = geneCounts[isexpr,]
genesExpr = noversion[isexpr,]

genes = DGEList( geneCountsExpr, genes=genesExpr ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
data.pm = cbind(data, mydata)
rownames(data.pm) = data$hbcc_brain_id
```

Now, let's see first if our main PM results are still there:

```r
num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data.pm[, x], data.pm[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(data.pm[, y], data.pm[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}

print(which(num_pvals < .01, arr.ind = T))
print(which(categ_pvals < .01, arr.ind = T))
```

```
                        row col
clusters                  2   1
RINe                      4   1
PMI                       5   1
RINe                      4   2
PMI                       5   2
clusters                  2   4
pcnt_optical_duplicates   1   7
clusters                  2   7
Age                       3   8
Age                       3   9
C4                        9   9
                row col
batch             1   1
batch             1   2
batch             1   7
MoD               3   9
substance_group   4   9
```

We had a few changes in the PCs: 1, 3, 5, 6, 8, 9.

```r
form = ~ Diagnosis + PC1 + PC3 + PC5 + PC6 + PC8 + PC9
design = model.matrix( form, data.pm)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_caudate = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

Let's do a quick run of GSEA:

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='caudate'
eval(parse(text=sprintf('res = rnaseq_%s', region)))

ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(geneid=res$GENEID, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

# my own GMTs
db = sprintf('my_%s_sets', region)
cat(region, db, '\n')
db_file = sprintf('~/data/post_mortem/%s.gmt', db)
enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                    organism="hsapiens",
                    enrichDatabaseFile=db_file,
                    enrichDatabaseType="genesymbol",
                    interestGene=tmp2,
                    outputDirectory = data_dir,
                    interestGeneType="ensembl_gene_id",
                    sigMethod="top", topThr=150000,
                    minNum=3,
                    isOutput=F, isParallel=T,
                    nThreads=ncpu, perNum=10000, maxNum=800))
out_fname = sprintf('%s/WG_ENSID_%s_%s_10K.csv', data_dir, region, db)
write.csv(enrichResult, file=out_fname, row.names=F)

DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5,
                                isOutput=F, isParallel=T,
                                nThreads=ncpu, perNum=10000)
    out_fname = sprintf('%s/WG_ENSID_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

And we do the same thing with the Caudate FUSION results:

```r
# bw
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_Caudate.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)

wei = read.delim('~/data/expression_impute/fusion_twas-master/WEIGHTS/Brain_Caudate_basal_ganglia.P01.pos')
cnt = 57
wei$geneid = substring(wei$WGT, cnt, cnt+14)
mf = merge(fusion, wei, by='ID', all.x=T, all.y=F, sort=F)
load('~/tmp/ensid3.rdata')
mf$adjPval = p.adjust(mf$TWAS.P, method='fdr')

library(ActivePathways)
gmt = read.GMT('~/data/post_mortem/acc_developmental.gmt')
junk = gmt[1:6]
g = mf[mf$adjPval < .05, 'geneid']
name = 'FDR 0.05'
junk[[1]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$adjPval < .1, 'geneid']
name = 'FDR 0.1'
junk[[2]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .05/nrow(mf), 'geneid']
name = 'Bonferroni'
junk[[3]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .001, 'geneid']
name = 'P < .001'
junk[[4]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .005, 'geneid']
name = 'P < .005'
junk[[5]] = list(id = name, genes = unique(g), name = name)
g = mf[mf$TWAS.P < .01, 'geneid']
name = 'P < .01'
junk[[6]] = list(id = name, genes = unique(g), name = name)
gmt_name = '~/data/post_mortem/FUSION_Caudate_sets.gmt'
write.GMT(junk, gmt_name)
```

And let's see if GSEA works here:

```r
# bw
library(WebGestaltR)
data_dir = '~/data/rnaseq_derek/'
ncpu=6
region='caudate'
eval(parse(text=sprintf('res = rnaseq_%s', region)))
ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(geneid=res$GENEID, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]
db_file = '~/data/post_mortem/FUSION_Caudate_sets.gmt'
enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                    organism="hsapiens",
                    enrichDatabaseFile=db_file,
                    enrichDatabaseType="ensembl_gene_id",
                    interestGene=tmp2,
                    outputDirectory = data_dir,
                    interestGeneType="ensembl_gene_id",
                    sigMethod="top", topThr=150000,
                    minNum=3,
                    isOutput=F, isParallel=T,
                    nThreads=ncpu, perNum=10000, maxNum=800))
out_fname = sprintf('%s/WG_ENSIDnzv_FUSION_%s_10K.csv', data_dir, region)
write.csv(enrichResult, file=out_fname, row.names=F)
```

We get P < .005 set nominally significant at .02, but that's it. For
over-representation...

```r
library(GeneOverlap)
both_res = merge(rnaseq_caudate, mf, by.x='GENEID', by.y='geneid',
                 all.x=F, all.y=F)
thresh = c(.05, .01, .005, .001)
imp_nums = vector(length=length(thresh), mode='numeric')
rna_nums = vector(length=length(thresh), mode='numeric')
pvals = matrix(data=NA, nrow=length(thresh), ncol=length(thresh))
rownames(pvals) = sapply(thresh, function(x) sprintf('PM_%.3f', x))
colnames(pvals) = sapply(thresh, function(x) sprintf('IMP_%.3f', x))
inter_nums = pvals
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$TWAS.P < thresh[ti], 'GENEID']
    imp_nums[ti] = length(imp_genes)
    for (tr in 1:length(thresh)) {
        rna_genes = both_res[both_res$P.Value < thresh[tr], 'GENEID']
        rna_nums[tr] = length(rna_genes)
        go.obj <- newGeneOverlap(imp_genes, rna_genes,
                                 genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(imp_genes, rna_genes)
        pval = getPval(go.obj)
        pvals[tr, ti] = pval
        inter_nums[tr, ti] = length(intersect(rna_genes, imp_genes))
    }
}
print(pvals)
print(inter_nums)
```

         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050 0.2476946 0.5458204  0.918706 0.5810571
PM_0.010 0.8615554 1.0000000  1.000000 1.0000000
PM_0.005 1.0000000 1.0000000  1.000000 1.0000000
PM_0.001 1.0000000 1.0000000  1.000000 1.0000000
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050        13         4         1         1
PM_0.010         1         0         0         0
PM_0.005         0         0         0         0
PM_0.001         0         0         0         0

Again nothing...

I'm still losing about 500 out of 2600 genes in ACC in the switch from hg19 to
hg38 (FUSION to PM). In Caudate it's about the same, losing about 700 out of 3500. 

Let's give this a break for now because I might be able to show better things
with DTE/DTU instead.

# TODO
 * Try same thing for the caudate... maybe better luck there?
 * Maybe try all genes, instead of just P01?
 * anything else we can do in the default FUSION prediction parameters?