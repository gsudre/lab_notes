# 2021-01-12 17:13:33

DESeq2 has some interesting capabilities, such as the inde[endent filtering and
IHW. Could we use it for DGE and DTE as well? I'm mostly focusing on:

http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html

and

https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

Let's try DGE first:

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

set.seed(42)
lcpm.pca <- prcomp(t(geneCounts), scale=TRUE)

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
pcnt_optical_duplicates   1   1
clusters                  2   1
C8                       13   1
RINe                      4   2
RINe                      4   3
pcnt_optical_duplicates   1   4
PMI                       5   4
C6                       11   7
               row col
batch            1   1
batch            1   2
batch            1   4
evidence_level   8   4
```

Debatable whether we should include PC7... let's leave it out for now:

```r
countdata = round(geneCounts)
colnames(countdata) = rownames(data.pm)
form = ~ Diagnosis + PC1 + PC2 + PC3 + PC4

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = data.pm,
                              design = form)
# removing rows of the DESeqDataSet that have no counts, or only a single count
# across all samples. Additional filtering will be applied later
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05)
```

```
r$> summary(res)                                                                        

out of 44290 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 6, 0.014%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

```r
library(IHW)
resIHW <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05,
                  filterFun=ihw)
summary(resIHW)
```

```
r$> summary(resIHW)                                                                     

out of 44290 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2, 0.0045%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting
```

I don't think this is particularly bad, assuming we can keep the 6 genes (if
they make sense), and we still have some meaningful GSEA results? There are 11
if we go for .1 FDR btw.

```
r$> rownames(res[res$padj < .05,])                                                      
[1] "ENSG00000103995.14" "ENSG00000135245.10" "ENSG00000240758.2"  "ENSG00000258864.1" 
[5] "ENSG00000258890.7"  "ENSG00000284084.1" 

r$> rownames(res[res$padj < .1,])                                                       
 [1] "ENSG00000002016.17" "ENSG00000078401.7"  "ENSG00000103995.14" "ENSG00000135245.10"
 [5] "ENSG00000169245.6"  "ENSG00000196584.3"  "ENSG00000240758.2"  "ENSG00000258864.1" 
 [9] "ENSG00000258890.7"  "ENSG00000268100.1"  "ENSG00000284084.1" 
```

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'

ranks = -log(res$pvalue) * sign(res$log2FoldChange)
tmp2 = data.frame(geneid=substring(rownames(res), 1, 15), rank=ranks)
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
                    nThreads=ncpu, perNum=1000, maxNum=800))
out_fname = sprintf('%s/WG_DESeq2_ENSID_%s_%s_1K.csv', data_dir, region, db)
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
                                nThreads=ncpu, perNum=1000)
    out_fname = sprintf('%s/WG_DESeq2_ENSID_%s_%s_1K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

Or maybe just using the log change as the rank is enough?

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'

ranks = res$log2FoldChange
tmp2 = data.frame(geneid=substring(rownames(res), 1, 15), rank=ranks)
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
                    nThreads=ncpu, perNum=1000, maxNum=800))
out_fname = sprintf('%s/WG_DESeq2LFC_ENSID_%s_%s_1K.csv', data_dir, region, db)
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
                                nThreads=ncpu, perNum=1000)
    out_fname = sprintf('%s/WG_DESeq2LFC_ENSID_%s_%s_1K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

In the meanwhile, let's plot the 6 good genes to see if they are just a fluke:

```r
quartz()
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Diagnosis"))
```

![](images/2021-01-12-17-47-43.png)

No, this is a no-go. It doesn't mean I can't go with DESSeq2, but we need to
make sure we remove these weird genes. Maybe using nzv?

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

set.seed(42)
lcpm.pca <- prcomp(t(geneCounts), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
data.pm = cbind(data, mydata)
rownames(data.pm) = data$hbcc_brain_id
```

Now, let's see how the PCs change after nzv:

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
pcnt_optical_duplicates   1   1
clusters                  2   1
C8                       13   1
RINe                      4   2
pcnt_optical_duplicates   1   4
PMI                       5   4
C6                       11   7
               row col
batch            1   1
batch            1   4
evidence_level   8   4
```

I still won't include 7:

```r
countdata = round(geneCounts)
colnames(countdata) = rownames(data.pm)
form = ~ Diagnosis + PC1 + PC2 + PC4

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = data.pm,
                              design = form)
# removing rows of the DESeqDataSet that have no counts, or only a single count
# across all samples. Additional filtering will be applied later
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

dds <- DESeq(dds)
resNZV <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05)
```

```
r$> summary(resNZV)                                                                     

out of 42245 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 7, 0.017%
LFC < 0 (down)     : 1, 0.0024%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 0)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

```r
quartz()
topGene <- rownames(resNZV)[which.min(resNZV$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Diagnosis"))
```

![](images/2021-01-12-17-53-45.png)

Still... but hard to tell if this is just by chance though... some stats will be
needed here.

# 2021-01-13 07:08:47

What if we do the same filtering as before, using edgeR?

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

set.seed(42)
lcpm.pca <- prcomp(t(geneCountsExpr), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
data.pm = cbind(data, mydata)
rownames(data.pm) = data$hbcc_brain_id
```

Now, let's see how the PCs change after nzv and filterByExpr:

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
pcnt_optical_duplicates   1   1
clusters                  2   1
C8                       13   1
RINe                      4   2
pcnt_optical_duplicates   1   3
PMI                       5   3
               row col
batch            1   1
batch            1   3
evidence_level   8   3
```

```r
countdata = round(geneCountsExpr)
colnames(countdata) = rownames(data.pm)
form = ~ Diagnosis + PC1 + PC2 + PC3

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = data.pm,
                              design = form)
dds <- DESeq(dds)
resNZV <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05)
```

Nothing survives. But if I use IHW, it does:

```r
library(IHW)
resIHW <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05,
                  filterFun=ihw)
```

```
    summary(resIHW)                                                                        

Attaching package: ‘IHW’

The following object is masked from ‘package:ggplot2’:

    alpha


out of 22203 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2, 0.009%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting
```

And there are 3 at q < .1. Let's see how the .05 look like:

```r
quartz()
topGene <- rownames(resIHW)[which.min(resIHW$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Diagnosis"))
topGene2 <- rownames(resIHW)[sort(resIHW$padj, index.return=T)$ix[2]]
r$> plotCounts(dds, gene = topGene2, intgroup=c("Diagnosis")) 
```

![](images/2021-01-13-07-17-56.png)
![](images/2021-01-13-07-19-45.png)

That is a bit better than all those zero counts from before. But do we keep out
GSEA results?

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'

ranks = -log(resIHW$pvalue) * sign(resIHW$log2FoldChange)
tmp2 = data.frame(geneid=substring(rownames(resIHW), 1, 15), rank=ranks)
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
                    nThreads=ncpu, perNum=1000, maxNum=800))
out_fname = sprintf('%s/WG_DESeq2_IHW_ENSID_%s_%s_1K.csv', data_dir, region, db)
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
                                nThreads=ncpu, perNum=1000)
    out_fname = sprintf('%s/WG_DESeq2_IHW_ENSID_%s_%s_1K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

Or maybe just using the log change as the rank is enough?

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'

ranks = resIHW$log2FoldChange
tmp2 = data.frame(geneid=substring(rownames(resIHW), 1, 15), rank=ranks)
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
                    nThreads=ncpu, perNum=1000, maxNum=800))
out_fname = sprintf('%s/WG_DESeq2LFC_IHW_ENSID_%s_%s_1K.csv', data_dir, region, db)
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
                                nThreads=ncpu, perNum=1000)
    out_fname = sprintf('%s/WG_DESeq2LFC_IHW_ENSID_%s_%s_1K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

These were actually the best results so far! And using logP was better than LFC.
I'll change the scripts above to run 10K then just so we have better results,
including the plots:

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
ncpu=6

region='acc'

ranks = -log(resIHW$pvalue) * sign(resIHW$log2FoldChange)
tmp2 = data.frame(geneid=substring(rownames(resIHW), 1, 15), rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

# my own GMTs
db = sprintf('my_%s_sets', region)
cat(region, db, '\n')
db_file = sprintf('~/data/post_mortem/%s.gmt', db)
project_name = sprintf('WG_DESeq2_IHW_ENSID_%s_%s_10K', region, db)
enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                    organism="hsapiens",
                    enrichDatabaseFile=db_file,
                    enrichDatabaseType="genesymbol",
                    interestGene=tmp2,
                    outputDirectory = data_dir,
                    interestGeneType="ensembl_gene_id",
                    sigMethod="top", topThr=150000,
                    minNum=3, projectName=project_name,
                    isOutput=T, isParallel=T,
                    nThreads=ncpu, perNum=10000, maxNum=800))
out_fname = sprintf('%s/WG_DESeq2_IHW_ENSID_%s_%s_10K.csv', data_dir, region, db)
write.csv(enrichResult, file=out_fname, row.names=F)

DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    project_name = sprintf('WG_DESeq2_IHW_ENSID_%s_%s_10K', region, db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5,
                                isOutput=T, isParallel=T,
                                projectName=project_name,
                                nThreads=ncpu, perNum=10000)
    out_fname = sprintf('%s/WG_DESeq2_IHW_ENSID_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

## DTE

If this pipeline works, then we should try it for DTE as well:

```r
load('~/data/isoforms/tximport_rsem_DTE.RData')
txi = rsem
myregion = 'ACC'

data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)
samples = data[, !grepl(colnames(data), pattern='^ENS')]
keep_me = colnames(txi$counts) %in% samples$submitted_name
for (i in c('abundance', 'counts', 'length')) {
    txi[[i]] = txi[[i]][, keep_me]
}
# sort samples to match order in rsem
rownames(samples) = samples$submitted_name
samples = samples[colnames(txi$counts), ]

samples$POP_CODE = as.character(samples$POP_CODE)
samples[samples$POP_CODE=='WNH', 'POP_CODE'] = 'W'
samples[samples$POP_CODE=='WH', 'POP_CODE'] = 'W'
samples$POP_CODE = factor(samples$POP_CODE)
samples$Individual = factor(samples$hbcc_brain_id)
samples[samples$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
samples[samples$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
samples$MoD = factor(samples$Manner.of.Death)
samples$batch = factor(as.numeric(samples$run_date))
samples$Diagnosis = factor(samples$Diagnosis, levels=c('Control', 'Case'))
samples$substance_group = factor(samples$substance_group)
samples$comorbid_group = factor(samples$comorbid_group_update)
samples$evidence_level = factor(samples$evidence_level)

# going back to the workflow... let's first determine our PCs
# remove rows that has zero counts across all samples
cts = round(txi$counts)
cts = cts[rowSums(cts) > 0,]
library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(cts), method = pp_order)
X = predict(pp, t(cts))
txCounts = t(X)

library(edgeR)
isexpr <- filterByExpr(txCounts, group=data$Diagnosis)
txCountsExpr = txCounts[isexpr,]

pca <- prcomp(t(txCountsExpr), scale=TRUE)

library(nFactors)
eigs <- pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(pca$x[, keep_me])
data.pm = cbind(samples, mydata)
rownames(data.pm) = data$hbcc_brain_id
```

Now, let's see which PCs to remove for DTE:

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
               'comorbid_group', 'POP_CODE', 'Sex')
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
pcnt_optical_duplicates   1   1
clusters                  2   1
C8                       13   1
RINe                      4   2
pcnt_optical_duplicates   1   3
PMI                       5   3
RINe                      4   4
C6                       11   4
      row col
batch   1   1
batch   1   2
batch   1   3
```

```r
countdata = round(txCountsExpr)
colnames(countdata) = rownames(data.pm)
form = ~ Diagnosis + PC1 + PC2 + PC3 + PC4

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = data.pm,
                              design = form)
dds <- DESeq(dds)
res <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05)
```

```
r$> summary(res)                                                                           

out of 70359 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 1, 0.0014%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 7)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results
```

Nothing survives. But if I use IHW, it does:

```r
library(IHW)
resIHW <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05,
                  filterFun=ihw)
```

```
r$> summary(resIHW)                                                                        

out of 70359 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 3, 0.0043%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 0, 0%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting
```

```r
quartz()
topGene <- rownames(resIHW)[which.min(resIHW$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Diagnosis"))
topGene2 <- rownames(resIHW)[sort(resIHW$padj, index.return=T)$ix[2]]
plotCounts(dds, gene = topGene2, intgroup=c("Diagnosis")) 
```

![](images/2021-01-13-07-32-57.png)
![](images/2021-01-13-07-33-16.png)

At least the results are not being driven by zeros.


# TODO
 * what if I added the PCs to the filtering function?
 * Caudate?
 * do DX*region under this framework
 * try Kallisto
 * try Salmon
 * check Science paper to see what other analysis we can do to co bine the two results
 * make -log10(p) for DGE against DTU and DGE against DTE
 * make Volcano plots for each modality? see other forms of plots from main tutorial