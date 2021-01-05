# 2020-12-29 14:03:08

After plotting several individual genes in the PM data I started getting
concerned that we have too many outliers in a per-gene basis. How do our results
change if I remove outliers, either by imputation or by winsorizing?

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
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]
# we're down from 60K to 38K samples by only looking at the ones with hgnc symbol. We might be losing too much here, so it's a step to reconsider in the future

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
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
geneCountsExpr = geneCounts[isexpr,]
genesExpr = G_list2[isexpr,]

mu = apply(geneCountsExpr, 1, mean)
sigma = apply(geneCountsExpr, 1, sd)
over = geneCountsExpr > (mu+3*sigma)
under = geneCountsExpr < (mu-3*sigma)
geneCountsExpr[over | under] = NA

library(VIM)
geneCountsExprImp = irmi(t(geneCountsExpr), imp_var=F) 



# this is taking forever... maybe irmi is not the best way to go here... maybe 
# winsorizing, medianImpute, or something else?

```

# 2020-12-31 12:35:42

Let's redo this using caret for imputation:

```r
out_Z = 5 # 3

mu = apply(count_matrix, 1, mean)
sigma = apply(count_matrix, 1, sd)
over = count_matrix > (mu + out_Z*sigma)
under = count_matrix < (mu - out_Z*sigma)
clean_count_matrix = count_matrix
clean_count_matrix[over | under] = NA

library(caret)
pp_order = c('zv', 'nzv', 'medianImpute')
pp = preProcess(t(clean_count_matrix), method = pp_order)
X = predict(pp, t(clean_count_matrix))
geneCounts = t(X)
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
geneCountsExpr = geneCounts[isexpr,]
genesExpr = G_list2[isexpr,]

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

# Grabbing PRS
fname = '~/data/post_mortem/genotyping/1KG/merged_PM_1KG_PRS_12032020.csv'
prs = read.csv(fname)
prs$hbcc_brain_id = sapply(prs$IID,
                          function(x) {
                              br = strsplit(x, '_')[[1]][2];
                              as.numeric(gsub(br, pattern='BR',
                                              replacement=''))})
imWNH = data$C1 > 0 & data$C2 < -.075
wnh_brains = data[which(imWNH),]$hbcc_brain_id

# using the most appropriate PRS, make sure we don't switch subject order
m = merge(data.pm, prs, by='hbcc_brain_id', sort=F)
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
m[, prs_names] = NA
keep_me = m$hbcc_brain_id %in% wnh_brains
m[keep_me, prs_names] = m[keep_me, 75:86]
m[!keep_me, prs_names] = m[!keep_me, 63:74]
data.prs = m[, c(1:61, 87:98)]

genes2 = genes[, data.pm$hbcc_brain_id %in% data.prs$hbcc_brain_id]
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
        res = cor.test(data[, x], mydata[, y])
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
        res = kruskal.test(mydata[, y], data[, x])
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
C6                       11   2
clusters                  2   3
C7                       12   3
pcnt_optical_duplicates   1   6
clusters                  2   6
pcnt_optical_duplicates   1   7
Age                       3   9
C4                        9   9
                row col
batch             1   1
batch             1   2
batch             1   6
batch             1   7
MoD               3   9
substance_group   4   9
```

We had a few changes in the PCs for Z=3: 1, 2, 3, 6, 7, 9.

```r
form = ~ Diagnosis + PC1 + PC2 + PC3 + PC6 + PC7 + PC9
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
tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
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
                    interestGeneType="genesymbol",
                    sigMethod="top", topThr=150000,
                    minNum=3, projectName=project_name,
                    isOutput=F, isParallel=T,
                    nThreads=ncpu, perNum=1000, maxNum=800))
if (class(enrichResult) != "try-error") {
    out_fname = sprintf('%s/WGrobust_%s_%s_1K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5, projectName=project_name,
                                isOutput=F, isParallel=T,
                                nThreads=ncpu, perNum=1000)
    out_fname = sprintf('%s/WGrobust_%s_%s_1K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

That wiped out all results... maybe defining outliers a bit more conservatively?
Say, Z = 5? Then, our PCs become: 1, 2, 3, 7, 8, 9.

```r
form = ~ Diagnosis + PC1 + PC2 + PC3 + PC7 + PC8 + PC9
design = model.matrix( form, data.pm)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

This worked better with a few of the usual hits coming back. Let's increase a
bit more? With Z = 10, our PCs become: 1, 2, 7, 8, 9.

```r
form = ~ Diagnosis + PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data.pm)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

Yeah, Z=10 basically gave me the same (good) results as before. 

Well, maybe this is not as big of an issue as I though... 

