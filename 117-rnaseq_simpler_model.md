# 2020-06-17 08:00:25

Let's see if we get similar gene sets results with a slightly simpler model.

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ]  # removing ACC outlier
rownames(data) = data$submitted_name  # just to ensure compatibility later
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
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

# Standard usage of limma/voom
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)
data$Individual = factor(data$hbcc_brain_id)
data$batch = factor(data$run_date)

design = model.matrix( ~ Region + batch , data)
vobj_tmp = voom( genes, design, plot=FALSE)
# apply duplicateCorrelation 
dupcor <- duplicateCorrelation(vobj_tmp,design,block=data$Individual)
# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
vobj = voom( genes, design, plot=FALSE, block=data$Individual,
             correlation=dupcor$consensus)

form =  ~ 0 + Region + Region:Diagnosis + batch + Sex + RINe + PMI + Age
design = model.matrix(form, data)
# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj, design, block=data$Individual)
# But this step uses only the genome-wide average for the random effect
fit <- lmFit(vobj, design, block=data$Individual, correlation=dupcor$consensus)

Lc = matrix(0, ncol=ncol(design))
colnames(Lc) = colnames(design)
# make the 2 contrast terms positive
Lc[length(Lc):(length(Lc)-1)] = 1
fitDupCor = contrasts.fit( fit, t(Lc))

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )
```

Now we check if our gene set results are still there:

```r
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
load('~/data/rnaseq_derek/enrich.RDATA')
load('~/data/rnaseq_derek/adhd_genesets_philip.RDATA')

resDC = topTable(fitDupCor, number=Inf) 
enrich_dupCor_camera = get_enrich_order2( resDC, geneSetsCombined ) 
adhd_dupCor_camera = get_enrich_order2( resDC, t2 ) 
```

## Going back to old school grouping

What if we go back to following the limma pipeline, like:

https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ]  # removing ACC outlier
rownames(data) = data$submitted_name  # just to ensure compatibility later
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
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

# Standard usage of limma/voom
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)
data$Individual = factor(data$hbcc_brain_id)
data$batch = as.factor(as.numeric(factor(data$run_date)))
DX2 = sapply(1:nrow(data), function(x) sprintf('%s_%s', data[x, 'Diagnosis'],
                                                data[x, 'Region']))
data$group = factor(DX2)

design = model.matrix( ~ group + batch , data)
vobj_tmp = voom( genes, design, plot=FALSE)
# apply duplicateCorrelation 
dupcor <- duplicateCorrelation(vobj_tmp,design,block=data$Individual)
# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
vobj = voom( genes, design, plot=FALSE, block=data$Individual,
             correlation=dupcor$consensus)

form =  ~ 0 + group + batch + Sex + RINe + PMI + Age
design = model.matrix(form, data)
# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj, design, block=data$Individual)
# But this step uses only the genome-wide average for the random effect
fit <- lmFit(vobj, design, block=data$Individual, correlation=dupcor$consensus)

contr <- makeContrasts((groupCase_ACC + groupCase_Caudate) -
                       (groupControl_ACC + groupControl_Caudate),
                       levels = colnames(coef(fit)))
fitDupCor = contrasts.fit( fit, contr )

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )

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
load('~/data/rnaseq_derek/enrich.RDATA')
load('~/data/rnaseq_derek/adhd_genesets_philip.RDATA')

resDC = topTable(fitDupCor, number=Inf) 
enrich_dupCor_camera = get_enrich_order2( resDC, geneSetsCombined ) 
adhd_dupCor_camera = get_enrich_order2( resDC, t2 ) 
```

I get a very similar results to running Gabriel's model... not exactly, probably
because I'm not modeling the interaction, and my voom is different. I also
didn't run the results for the developmental gene sets.

We can also just check whether our results are still there if we add the
intercept back in the old analysis:

```r
design = model.matrix( ~ Region + batch , data)
vobj_tmp = voom( genes, design, plot=FALSE)
# apply duplicateCorrelation 
dupcor <- duplicateCorrelation(vobj_tmp,design,block=data$Individual)
# run voom considering the duplicateCorrelation results
# in order to compute more accurate precision weights
vobj = voom( genes, design, plot=FALSE, block=data$Individual,
             correlation=dupcor$consensus)

form =  ~ Region + Region:Diagnosis + batch + Sex + RINe + PMI + Age
design = model.matrix(form, data)
# Estimate linear mixed model with a single variance component
# Fit the model for each gene, 
dupcor <- duplicateCorrelation(vobj, design, block=data$Individual)
# But this step uses only the genome-wide average for the random effect
fit <- lmFit(vobj, design, block=data$Individual, correlation=dupcor$consensus)

Lc = matrix(0, ncol=ncol(design))
colnames(Lc) = colnames(design)
# make the 2 contrast terms positive
Lc[length(Lc):(length(Lc)-1)] = 1
fitDupCor = contrasts.fit( fit, t(Lc))

# Fit Empirical Bayes for moderated t-statistics
fitDupCor <- eBayes( fitDupCor )
resDC = topTable(fitDupCor, number=Inf) 
enrich_dupCor_camera = get_enrich_order2( resDC, geneSetsCombined ) 
adhd_dupCor_camera = get_enrich_order2( resDC, t2 ) 
```

Results were exactly the same!

Let me check if the results using DX2 are only affected because of voom:

```r
DX2 = sapply(1:nrow(data), function(x) sprintf('%s_%s', data[x, 'Diagnosis'],
                                                data[x, 'Region']))
data$group = factor(DX2)

design = model.matrix( ~ Region + batch , data)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=data$Individual)
vobj = voom( genes, design, plot=FALSE, block=data$Individual,
             correlation=dupcor$consensus)

form =  ~ 0 + group + batch + Sex + RINe + PMI + Age
design = model.matrix(form, data)
dupcor <- duplicateCorrelation(vobj, design, block=data$Individual)
fit <- lmFit(vobj, design, block=data$Individual, correlation=dupcor$consensus)

contr <- makeContrasts((groupCase_ACC + groupCase_Caudate) -
                       (groupControl_ACC + groupControl_Caudate),
                       levels = colnames(coef(fit)))
fitDupCor = contrasts.fit( fit, contr )
fitDupCor <- eBayes( fitDupCor )
resDC = topTable(fitDupCor, number=Inf) 
enrich_dupCor_camera = get_enrich_order2( resDC, geneSetsCombined ) 
adhd_dupCor_camera = get_enrich_order2( resDC, t2 ) 
```

They're not exactly the same, but very similar. Like, differences like .0348 and
0.0379.

We could also test multiple coefficients at the same time and look at the F
stats. For example, use two separate contrasts:

```r
DX2 = sapply(1:nrow(data), function(x) sprintf('%s_%s', data[x, 'Diagnosis'],
                                                data[x, 'Region']))
data$group = factor(DX2)

design = model.matrix( ~ Region + batch , data)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp,design,block=data$Individual)
vobj = voom( genes, design, plot=FALSE, block=data$Individual,
             correlation=dupcor$consensus)

form =  ~ 0 + group + batch + Sex + RINe + PMI + Age
design = model.matrix(form, data)
dupcor <- duplicateCorrelation(vobj, design, block=data$Individual)
fit <- lmFit(vobj, design, block=data$Individual, correlation=dupcor$consensus)

contrACC <- makeContrasts(groupCase_ACC - groupControl_ACC,
                          levels = colnames(coef(fit)))
contrCaudate <- makeContrasts(groupCase_Caudate - groupControl_Caudate,
                              levels = colnames(coef(fit)))
contr = cbind(contrACC, contrCaudate)
fitDupCor = contrasts.fit( fit, contr )
fitDupCor <- eBayes( fitDupCor )
resDC = topTable(fitDupCor, number=Inf) 
enrich_dupCor_camera = get_enrich_order2( resDC, geneSetsCombined ) 
adhd_dupCor_camera = get_enrich_order2( resDC, t2 ) 
```

This works and I get a F-value because I'm running two coefficients/comparisons,
but it's the end result is not great.

ADHD set result is somewhat there if I look only at the ACC comparison, but not
as strong as using the uber contrast, or Gabriel's mode.

What if I try something like example 3.5 in 
https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

```r
# design <- model.matrix(~Individual + batch + Sex + RINe + PMI + Age, data=data)
design <- model.matrix(~Individual + batch, data=data)
ACC.Case <- data$Diagnosis=="Case" & data$Region=="ACC"
Caudate.Case <- data$Diagnosis=="Case" & data$Region=="Caudate"
design <- cbind(design, ACC.Case, Caudate.Case)

y <- estimateDisp(genes, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef="ACC.Case")
topTags(qlf)
```

Nevermind... the design matrix is not full rank. Is it because of the people
that don't have ACC?

```r
idx = data$Individual %in% c(1530, 1683, 2485, 3003, 3006, 3006, 3023, 3055)
data2 = data[!idx, ]
data2$Individual = factor(data2$Individual)
design <- model.matrix(~Individual + batch, data=data2)
ACC.Case <- data2$Diagnosis=="Case" & data2$Region=="ACC"
Caudate.Case <- data2$Diagnosis=="Case" & data2$Region=="Caudate"
design <- cbind(design, ACC.Case, Caudate.Case)
y <- estimateDisp(genes[, !idx], design)
```

Still not singular... oh well.

But it's also important to read 3.2.4 to see if other contrasts would make more
sense to run with limma and random subjects.