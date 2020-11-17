# 2020-11-17 09:53:22

Looking at that paper Philip sent:

"Genome-wide changes in lncRNA, splicing, and regional gene expression patterns
in autism" 

gave me the idea to look at long non-coding RNA in our data. In fact, at al the
non-coding RNA. Maybe there's something there? That's the step where we go from
60 to 38K markers, so we might be losing information here.

If anything, we'll need to get some good hits statistically because I'm not sure
if we'll be able to recoup much stuff using gene sets here.

Let's run the same analysis as before, but only keep the non-anotated genes instead.

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
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
keep_me = is.na(G_list0$hgnc_symbol) | G_list0$hgnc_symbol==''
G_list <- G_list0[keep_me,]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
has_anno = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[has_anno, ]
# we're now working with 22K markers

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

library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)
# down to 15K transcripts
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]
# down to 14K transcripts
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)
# down to 4K genes

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
```

We need to do the same PCA investigation here. But I did check randomly some of
the ENSG entries in our current genes matrix, and they do indeed map to lncRNA
according to GeneCards.

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
               'comorbid_group', 'POP_CODE', 'Sex')
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
```

```
r$> which(num_pvals < .01, arr.ind = T)
                         row col
clusters                  2   1
Age                       3   1
PMI                       5   1
RINe                      4   2
C6                       11   2
pcnt_optical_duplicates   1   6
clusters                  2   6
C1                        6   9
C2                        7   9
C3                        8   9
C4                        9   9
Age                       3  10
C1                        6  10
C2                        7  10
C3                        8  10
C4                        9  10
Age                       3  13
C7                       12  13

r$> which(categ_pvals < .01, arr.ind = T)
                row col
batch             1   1
batch             1   6
substance_group   4   9
POP_CODE          6   9
Diagnosis         2  13
MoD               3  13
```

PC13 is related to Diagnosis and MoD, Age and C7. Let's keep it in.

```r
data2 = cbind(data, mydata)
form = ~ Diagnosis + PC1 + PC2 + PC6 + PC9 + PC10
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
ncrna_acc = rnaseq_acc[order(rnaseq_acc$P.Value),]
```

We do get a result under q< .05, but it'sa novel protein... maybe we could find
some sort of list for these markers? Maybe a developmental list like in the
paper?

Let's repeat it for Caudate:

```r
myregion = 'Caudate'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 58 samples for Caudate
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
keep_me = is.na(G_list0$hgnc_symbol) | G_list0$hgnc_symbol==''
G_list <- G_list0[keep_me,]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
has_anno = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[has_anno, ]
# we're now working with 22K markers

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

library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)
# down to 15K transcripts
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]
# down to 14K transcripts
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)
# down to 5K genes

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])

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
               'comorbid_group', 'POP_CODE', 'Sex')
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
```

```
r$> which(num_pvals < .01, arr.ind = T)
         row col
clusters   2   1
RINe       4   1
PMI        5   1
RINe       4   3
clusters   2   6
clusters   2   8
C7        12   8
C1         6  10
C4         9  10
C1         6  12
C2         7  12
C3         8  12
C4         9  12
C9        14  12
C10       15  12
Age        3  13

r$> which(categ_pvals < .01, arr.ind = T)
         row col
batch      1   1
batch      1   6
batch      1   7
batch      1   8
MoD        3   8
POP_CODE   6  10
POP_CODE   6  12
```

PC13 is related to Diagnosis at about .02, so I'll keep it.

```r
data2 = cbind(data, mydata)
form = ~ Diagnosis + PC1 + PC3 + PC6 + PC8 + PC10 + PC12 + PC7
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_caudate = topTable(fit2, coef='DiagnosisCase', number=Inf)
ncrna_caudate = rnaseq_caudate[order(rnaseq_caudate$P.Value),]
```

Nothing significant here though.

# TODO
 * overlaps between caudate and ACC?
 * compare to some sort of developmental expression list like in the paper?