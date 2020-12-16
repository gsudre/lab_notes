# 2020-12-16 14:03:55

Even though it's a bit circular, we could also check whether the PRS is
predictive of the imputed data. 

```r
data = read.delim('~/data/expression_impute/results/NCR_v3_ACC_summary_1KG_mashr.txt')
```









```r
a = read.table('~/data/expression_impute/results/NCR_v3_ACC_predict_1KG_mashr.txt', header=1)
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)

imwnh = imp_data[imp_data$PC01<0 & imp_data$PC02>-.02,]$IID
data_dir = '~/data/expression_impute/'
phenotypes = list(ACC=c('res_ACC_thickness'),
                  Caudate=c('res_Caudate_volume'))
for (region in c('ACC', 'Caudate')) {
   for (my_phen in phenotypes[[region]]) {
       print(my_phen)
      data3 = data2[data2$SID %in% imwnh, ]
      data3 = data3[data3$bestInFamily==T, ]
      data3 = data3[, c('SID', my_phen, 'Sex...Subjects')]
      colnames(data3)[1] = 'IID'
      colnames(data3)[2] = 'phen'
      colnames(data3)[3] = 'sex'
      data3$sex = as.numeric(as.factor(data3$sex))
      data3 = data3[order(data3$IID), ]
      # it expects no more than the number of people we have in the phenotypes
      a = read.table(sprintf('%s/results/NCR_v3_%s_predict_1KG_en.txt',
                             data_dir, region), header=1)
    #    a = readRDS(sprintf('%s/results/NCR_v3_%s_1KG_mashr.rds', data_dir,
    #                        region))
      # remove FAMID from IID
      iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
      iid3 = gsub(x=iid2, pattern='SID.', replacement='')
      a$IID = as.numeric(iid3)
      b = a[a$IID %in% data3$IID, ]
      b = b[order(b$IID), ]
      data3$FID = b$FID # they're both sorted on IID
      write.table(b, file=sprintf('%s/ANAT_cropped_imp_EN_%s.tab', data_dir,
                                  region), row.names=F, quote=F, sep='\t')
    #   write.table(b, file=sprintf('%s/ANAT_cropped_imp_MASHR_%s.tab', data_dir,
    #                               region), row.names=F, quote=F, sep='\t')
      write.table(data3, file=sprintf('%s/phen_%s.tab', data_dir, my_phen),
                  row.names=F, quote=F, sep='\t')
   }
}
```



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
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])

data2 = cbind(data, mydata)
# we don't residualize Diagnosis!

form = ~ PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
resids = residuals(fit2, genes)
```

Let's first break down the 3 different sets of PRS we will try:

```r
fname = '~/data/post_mortem/genotyping/1KG/merged_PM_1KG_PRS_12032020.csv'
prs = read.csv(fname)
prs$hbcc_brain_id = sapply(prs$IID,
                          function(x) {
                              br = strsplit(x, '_')[[1]][2];
                              as.numeric(gsub(br, pattern='BR',
                                              replacement=''))})
imWNH = data$C1 > 0 & data$C2 < -.075
wnh_brains = data[which(imWNH),]$hbcc_brain_id

# using whole population PRS
m = merge(data, prs, by='hbcc_brain_id')
m = m[, 1:63]
colnames(m)[52:63] = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                              .5, .4, .3, .2),
                            function(x) sprintf('PRS%f', x))
data_whole = m

# WNH samples only
m = merge(data, prs, by='hbcc_brain_id')
m = m[, c(1:51, 64:75)]
colnames(m)[52:63] = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                              .5, .4, .3, .2),
                            function(x) sprintf('PRS%f', x))
keep_me = m$hbcc_brain_id %in% wnh_brains
data_wnh = m[keep_me, ]

# using the most appropriate PRS
m = merge(data, prs, by='hbcc_brain_id')
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
m[, prs_names] = NA
keep_me = m$hbcc_brain_id %in% wnh_brains
m[keep_me, prs_names] = m[keep_me, 64:75]
m[!keep_me, prs_names] = m[!keep_me, 52:63]
data_app = m[, c(1:51, 76:87)]
```

And of course we'll need to try this for all possible PRS values. Let's try the
PCA approach first, as it's easier to quantify the results.

Now that I'm re-reading the paper, they only took the first PC in Fig 3 to
compare the different modalities. For PRS, they used the formula in the
supplement. So, let's just go with that. I can use the formula to obtain the
best PRS threshold, and start examining that, but I'll likely eventually run all
of them just for comparison.

Actually, I think it's both. To come up with a single R2, you need one summary
value, so in that case it has to be the first PC.

```r
library('fmsb')
library(lmtest)
mod.baseline = glm(Diagnosis ~ Sex + Age, family=binomial, data=data_whole)
mod.full = glm(Diagnosis ~ PRS0.000100 + Sex + Age, family=binomial, data=data_whole)
adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
prs.significance = lrtest(mod.baseline, mod.full)
```

So, if we are goig to try all of them:

```r
mydata = data_wnh

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
for (prs in prs_names) {
    # fm_root = 'Diagnosis ~ %s Sex + Age'
    fm_root = 'Diagnosis ~ %s Sex + Age + C1 + C2 + C3 + C4 + C5'
    fm_str = sprintf(fm_root, '')
    mod.baseline = glm(as.formula(fm_str), family=binomial, data=mydata)
    fm_str = sprintf(fm_root, sprintf('%s +', prs))
    mod.full = glm(as.formula(fm_str), family=binomial, data=mydata)
    adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
    prs.significance = lrtest(mod.baseline, mod.full)
    myp = prs.significance$"Pr(>Chisq)"[2] 
    cat(sprintf('%s: pval = %.3f\n', fm_str, myp))
}
```

```
(whole)
Diagnosis ~ PRS0.000100 + Sex + Age: pval = 0.980
Diagnosis ~ PRS0.001000 + Sex + Age: pval = 0.424
Diagnosis ~ PRS0.010000 + Sex + Age: pval = 0.226
Diagnosis ~ PRS0.100000 + Sex + Age: pval = 0.706
Diagnosis ~ PRS0.000050 + Sex + Age: pval = 0.933
Diagnosis ~ PRS0.000500 + Sex + Age: pval = 0.289
Diagnosis ~ PRS0.005000 + Sex + Age: pval = 0.285
Diagnosis ~ PRS0.050000 + Sex + Age: pval = 0.209
Diagnosis ~ PRS0.500000 + Sex + Age: pval = 0.715
Diagnosis ~ PRS0.400000 + Sex + Age: pval = 0.834
Diagnosis ~ PRS0.300000 + Sex + Age: pval = 0.753
Diagnosis ~ PRS0.200000 + Sex + Age: pval = 0.501

(appropriate)
Diagnosis ~ PRS0.000100 + Sex + Age: pval = 0.519
Diagnosis ~ PRS0.001000 + Sex + Age: pval = 0.069
Diagnosis ~ PRS0.010000 + Sex + Age: pval = 0.150
Diagnosis ~ PRS0.100000 + Sex + Age: pval = 0.273
Diagnosis ~ PRS0.000050 + Sex + Age: pval = 0.062
Diagnosis ~ PRS0.000500 + Sex + Age: pval = 0.178
Diagnosis ~ PRS0.005000 + Sex + Age: pval = 0.041
Diagnosis ~ PRS0.050000 + Sex + Age: pval = 0.107
Diagnosis ~ PRS0.500000 + Sex + Age: pval = 0.681
Diagnosis ~ PRS0.400000 + Sex + Age: pval = 0.601
Diagnosis ~ PRS0.300000 + Sex + Age: pval = 0.837
Diagnosis ~ PRS0.200000 + Sex + Age: pval = 0.842

(WNH)
Diagnosis ~ PRS0.000100 + Sex + Age: pval = 0.670
Diagnosis ~ PRS0.001000 + Sex + Age: pval = 0.910
Diagnosis ~ PRS0.010000 + Sex + Age: pval = 0.802
Diagnosis ~ PRS0.100000 + Sex + Age: pval = 0.744
Diagnosis ~ PRS0.000050 + Sex + Age: pval = 0.624
Diagnosis ~ PRS0.000500 + Sex + Age: pval = 0.854
Diagnosis ~ PRS0.005000 + Sex + Age: pval = 0.586
Diagnosis ~ PRS0.050000 + Sex + Age: pval = 0.618
Diagnosis ~ PRS0.500000 + Sex + Age: pval = 0.968
Diagnosis ~ PRS0.400000 + Sex + Age: pval = 0.934
Diagnosis ~ PRS0.300000 + Sex + Age: pval = 0.999
Diagnosis ~ PRS0.200000 + Sex + Age: pval = 0.941
```

It looks like using the "Appropriate" dataset has better results in predicting
Diagnosis. But there are other cofounders there too, such as comorbidity and
substance abuse, that are not being taken into this analysis. Let's at least add
the population PCs just to mirror the Science paper a bit more and see if they
have any difference. Just the first 5:

```
(whole)
Diagnosis ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.435
Diagnosis ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.380
Diagnosis ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.117
Diagnosis ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.092
Diagnosis ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.542
Diagnosis ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.585
Diagnosis ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.054
Diagnosis ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.194
Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.015
Diagnosis ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.018
Diagnosis ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.045
Diagnosis ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.138

(appropriate)
Diagnosis ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.676
Diagnosis ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.471
Diagnosis ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.140
Diagnosis ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.057
Diagnosis ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.405
Diagnosis ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.528
Diagnosis ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.013
Diagnosis ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.125
Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.004
Diagnosis ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.009
Diagnosis ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.019
Diagnosis ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.065

(WNH)
Diagnosis ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.356
Diagnosis ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.208
Diagnosis ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.645
Diagnosis ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.244
Diagnosis ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.662
Diagnosis ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.276
Diagnosis ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.019
Diagnosis ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.201
Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.436
Diagnosis ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.537
Diagnosis ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.559
Diagnosis ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.493
```

Nothing to write home about... maybe there is something there in WNH-only, but
that's only 30 subjects. If we're doing a WNH analysis, I think we'd need to do
it from scratch, including the initial removal of PCs?

Just because I'll wonder about this later, doing resids(fit) = resids(fit2), so
we don't need to worry about doing it in the results of lmFit or Bayes.

But let's check on the individual genes. Anything interesting there? Let's go
with our best PRS and then I can try other things:

```r
lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
rownames(mydata) = data$hbcc_brain_id

data2 = merge(data_app, mydata, by.x='hbcc_brain_id', by.y=0)
genes2 = genes[, data$hbcc_brain_id %in% data2$hbcc_brain_id]
form = ~ PRS0.500000 + PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes2, design, plot=FALSE)
prs.fit <- lmFit(vobj, design)
prs.fit2 <- eBayes( prs.fit )
res = topTable(prs.fit2, coef='PRS0.500000', number=Inf)
```

OK, now I have lots of results. What to do? Maybe first we should check if there
is a significant overlap between these genes and the PM_ACC genes?

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
for (t in c(.05, .01, .005, .001)) {
    prs_genes = res[res$P.Value < t, 'hgnc_symbol']
    dx_genes = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
    go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
    go.obj <- testGeneOverlap(go.obj)
    inter = intersect(prs_genes, dx_genes)
    pval = getPval(go.obj)
    cat(sprintf('t=%.3f, prs=%d, pm=%d, in=%d, p=%f\n', t,
                length(prs_genes), length(dx_genes), length(inter), pval))
}
```

This is a potentially cool result. 

```
t=0.050, prs=1504, pm=1335, in=202, p=0.000000
t=0.010, prs=305, pm=325, in=19, p=0.000004
t=0.005, prs=168, pm=165, in=10, p=0.000004
t=0.001, prs=42, pm=38, in=1, p=0.086516
```

So, these are quite cool. Anything in the caudate?

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
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
rownames(mydata) = data$hbcc_brain_id

data2 = merge(data_app, mydata, by.x='hbcc_brain_id', by.y=0)
genes2 = genes[, data$hbcc_brain_id %in% data2$hbcc_brain_id]
form = ~ PRS0.500000 + PC1 + PC3 + PC5 + PC6 + PC8
design = model.matrix( form, data2)
vobj = voom( genes2, design, plot=FALSE)
prs.fit <- lmFit(vobj, design)
prs.fit2 <- eBayes( prs.fit )
res2 = topTable(prs.fit2, coef='PRS0.500000', number=Inf)
```

Then we run the same over-representation analysis:

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
for (t in c(.05, .01, .005, .001)) {
    prs_genes = res2[res2$P.Value < t, 'hgnc_symbol']
    dx_genes = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
    go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
    go.obj <- testGeneOverlap(go.obj)
    inter = intersect(prs_genes, dx_genes)
    pval = getPval(go.obj)
    cat(sprintf('t=%.3f, prs=%d, pm=%d, in=%d, p=%f\n', t,
                length(prs_genes), length(dx_genes), length(inter), pval))
}
```

And this is the Caudate... interesting:

```
t=0.050, prs=1339, pm=1063, in=72, p=0.860254
t=0.010, prs=264, pm=220, in=1, p=0.964251
t=0.005, prs=125, pm=117, in=0, p=1.000000
t=0.001, prs=27, pm=26, in=0, p=1.000000
```

# 2020-12-04 09:49:11

Might as well save the results for all PRS thresholds:

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC2 + PC7 + PC8 + PC9', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t, 'hgnc_symbol']
        dx_genes = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval = getPval(go.obj)
        this_res = c(p, t, length(prs_genes), length(dx_genes), length(inter),
                     pval)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('PRS', 'nomPvalThresh', 'PRsgenes', 'PMgenes',
                      'overlap', 'pval')
write.csv(all_res, file='~/data/post_mortem/all_acc_prs_overlap_results.csv',
          row.names=F)
```

And run the same thing for Caudate:

```r
all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC3 + PC5 + PC6 + PC8', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t, 'hgnc_symbol']
        dx_genes = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval = getPval(go.obj)
        this_res = c(p, t, length(prs_genes), length(dx_genes), length(inter),
                     pval)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('PRS', 'nomPvalThresh', 'PRsgenes', 'PMgenes',
                      'overlap', 'pval')
write.csv(all_res, file='~/data/post_mortem/all_caudate_prs_overlap_results.csv',
          row.names=F)
```

# TODO
* do WNH analysis from scratch?