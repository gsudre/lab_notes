# 2020-06-22 19:54:17

Let's try a few things I chatted with Philip and wrote down in 120.

## Re-coding ACC outlier

That's the simplest thing to do. Let's see if our WNH results are still there if
instead of removing that one sample we just recorde it as Caudate.

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# change obvious outlier that's likely caudate labeled as ACC
change_me = which(rownames(data) %in% c('68080'))
data[change_me, 'Region'] = 'Caudate'

imWNH = which(data$C1 > 0 & data$C2 < -.075)
data = data[imWNH, ]

# ending up with 62 data points across both brain regions (before we had 61)
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
# data matrix goes on a diet...
data = data[, !grepl(colnames(data), pattern='^ENS')]
# remove that weird .num after ENSG
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
# remove any genes without a HUGOID
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
# remove genes that appear more than once
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
# keep only gene counts for genes that we have information
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]

# some data variables modifications
data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))

library(caret)
set.seed(42)
# remove genes with zero or near zero variance so we can run PCA
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)

# match gene counts to gene info
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'

# keep only autosomal genes
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]

library(edgeR)
isexpr = rowSums(cpm(geneCounts)>1) >= 0.1*ncol(geneCounts)

genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

library(variancePartition)
library(BiocParallel)
param = SnowParam(32, "SOCK", progressbar=TRUE)
register(param)

get_enrich_order2 = function( res, gene_sets ){
  if( !is.null(res$z.std) ){
    stat = res$z.std
  }else if( !is.null(res$F.std) ){
    stat = res$F.std
  }else if( !is.null(res$t) ){
    stat = res$t
  }else{
    stat = res$F
  }
  names(stat) = res$hgnc_symbol
  stat = stat[!is.na(names(stat))]
  # print(head(stat))
  index = ids2indices(gene_sets, names(stat))
  cameraPR( stat, index )
}
load('~/data/rnaseq_derek/adhd_genesets_philip.RDATA')

form = ~ (1|Region) + (1|Individual) + (1|batch)
vobjMM = voomWithDreamWeights( genes, form, data, plot=FALSE)
form = ~ 0 + Region + Region:Diagnosis + (1|Individual) + (1|batch) + Sex + scale(RINe) + scale(PMI) + scale(Age) + MoD
L = getContrast( vobjMM, form, data, "RegionCaudate:DiagnosisControl")
L['RegionACC:DiagnosisControl'] = 1
fitmm = dream( vobjMM, form, data, L)
resMM = topTable(fitmm, coef="L1", number=Inf) 
adhd_dream_camera = get_enrich_order2( resMM, t2 ) 
load('~/data/rnaseq_derek/dev_gene_sets.RData')
dev_dream_camera = get_enrich_order2( resMM, genes_unique )
load('~/data/rnaseq_derek/c5_gene_sets.RData')
c5_dream_camera = get_enrich_order2( resMM, c5_all)
load('~/data/rnaseq_derek/brain_disorders_gene_sets.RData')
dis_dream_camera = get_enrich_order2( resMM, disorders)
```

That's interesting... that didn't look like a Caudate sample. None of the
results held by a mile when adding it as Caudate!!!

# Adding more WNH

Let's play a bit with the WNH threshold to see if we can include a few more.
Basically, we'll raise the y axis limit a bit to first only leave out 5 dots,
then only 2 (the two highest ones in the Y axis).

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier that's likely caudate labeled as ACC
rm_me = rownames(data) %in% c('68080')
ggplot(data, aes(x=C1, y=C2, col=POP_CODE)) + geom_point() + geom_hline(yintercept=-.07, linetype="dashed", color = "black") + geom_vline(xintercept=0, linetype="dashed", color = "black")
```

![](images/2020-06-22-20-20-54.png)

```r
ggplot(data, aes(x=C1, y=C2, col=POP_CODE)) + geom_point() + geom_hline(yintercept=-.05, linetype="dashed", color = "black") + geom_vline(xintercept=0, linetype="dashed", color = "black")
```

![](images/2020-06-22-20-21-42.png)

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier that's likely caudate labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]

imWNH = which(data$C1 > 0 & data$C2 < -.07)
data = data[imWNH, ]

# ending up with 62 data points across both brain regions (before we had 61)
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
# data matrix goes on a diet...
data = data[, !grepl(colnames(data), pattern='^ENS')]
# remove that weird .num after ENSG
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
# remove any genes without a HUGOID
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
# remove genes that appear more than once
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
# keep only gene counts for genes that we have information
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]

# some data variables modifications
data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))

library(caret)
set.seed(42)
# remove genes with zero or near zero variance so we can run PCA
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)

# match gene counts to gene info
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'

# keep only autosomal genes
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]

library(edgeR)
isexpr = rowSums(cpm(geneCounts)>1) >= 0.1*ncol(geneCounts)

genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

library(variancePartition)
library(BiocParallel)
param = SnowParam(32, "SOCK", progressbar=TRUE)
register(param)

get_enrich_order2 = function( res, gene_sets ){
  if( !is.null(res$z.std) ){
    stat = res$z.std
  }else if( !is.null(res$F.std) ){
    stat = res$F.std
  }else if( !is.null(res$t) ){
    stat = res$t
  }else{
    stat = res$F
  }
  names(stat) = res$hgnc_symbol
  stat = stat[!is.na(names(stat))]
  # print(head(stat))
  index = ids2indices(gene_sets, names(stat))
  cameraPR( stat, index )
}
load('~/data/rnaseq_derek/adhd_genesets_philip.RDATA')

form = ~ (1|Region) + (1|Individual) + (1|batch)
vobjMM = voomWithDreamWeights( genes, form, data, plot=FALSE)
form = ~ 0 + Region + Region:Diagnosis + (1|Individual) + (1|batch) + Sex + scale(RINe) + scale(PMI) + scale(Age) + MoD
L = getContrast( vobjMM, form, data, "RegionCaudate:DiagnosisControl")
L['RegionACC:DiagnosisControl'] = 1
fitmm = dream( vobjMM, form, data, L)
resMM = topTable(fitmm, coef="L1", number=Inf) 
adhd_dream_camera = get_enrich_order2( resMM, t2 ) 
load('~/data/rnaseq_derek/dev_gene_sets.RData')
dev_dream_camera = get_enrich_order2( resMM, genes_unique )
load('~/data/rnaseq_derek/c5_gene_sets.RData')
c5_dream_camera = get_enrich_order2( resMM, c5_all)
load('~/data/rnaseq_derek/brain_disorders_gene_sets.RData')
dis_dream_camera = get_enrich_order2( resMM, disorders)
```

![](images/2020-06-22-21-08-06.png)

![](images/2020-06-22-21-11-00.png)

With more WNH I still get the TWAS result, just not as significant. I do lose
the developmental result though. The c5 results are not that good either,
loosing all FDR < .05 and the top results are not even that good.

![](images/2020-06-22-21-09-24.png)

![](images/2020-06-22-21-11-57.png)

With less WNH (still more than the original more strict set) is the same effect,
but even worse results. I still have some c5 results using the more strict set,
but it's not great.

What if I use the same model in voom? Does it change anything?

```r
form = ~ 0 + Region + Region:Diagnosis + (1|Individual) + (1|batch) + Sex + scale(RINe) + scale(PMI) + scale(Age) + MoD
vobjMM = voomWithDreamWeights( genes, form, data, plot=FALSE)
L = getContrast( vobjMM, form, data, "RegionCaudate:DiagnosisControl")
L['RegionACC:DiagnosisControl'] = 1
fitmm = dream( vobjMM, form, data, L)
resMM = topTable(fitmm, coef="L1", number=Inf) 
adhd_dream_camera = get_enrich_order2( resMM, t2 ) 
load('~/data/rnaseq_derek/dev_gene_sets.RData')
dev_dream_camera = get_enrich_order2( resMM, genes_unique )
load('~/data/rnaseq_derek/c5_gene_sets.RData')
c5_dream_camera = get_enrich_order2( resMM, c5_all)
load('~/data/rnaseq_derek/brain_disorders_gene_sets.RData')
dis_dream_camera = get_enrich_order2( resMM, disorders)
```

