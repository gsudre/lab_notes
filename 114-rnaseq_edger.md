# 2020-06-02 07:16:45

I want to try a few examples using the edgeR user-guide, just to be safe that
we're running the right commands. No COMBAT, at least not for now.

https://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf

I'll try them in order, based on what I identified in note 112:

## 3.4.2

This considers all subjects as independent. Not true, but it's a start. Here's
something from the user guide about their filtering function:

```
The filtering should be based on the grouping factors or treatment factors that will be involved in the differential expression tested for, rather than on blocking variables that are not of scientific interest in themselves. For example, consider a paired comparison experiment in which the same treatment regimes applied to each of a number of subjects or patients. In this design, Patient is included in the design matrix to correct for baseline differences between the Patients, but we will not be testing for differential expression between the Patients. The filtering should therefore be based soley Treatment rather than on Patient
```

That's a similar situation to what we have, meaning that we should filter on
Diagnosis but not Region. Of course I'll try different filtering later, to
assess for result robustness, but for now I'll just go with their suggestion. 

```r
library(edgeR)
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ]  # removing ACC outlier
rownames(data) = data$submitted_name  # just to ensure compatibility later

grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
# remove that weird .num after ENSG
id_num = sapply(grex_vars,
                function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num

dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]
imautosome = which(G_list$chromosome_name != 'X' &
                   G_list$chromosome_name != 'Y' &
                   G_list$chromosome_name != 'MT')
count_matrix = count_matrix[imautosome, ]
G_list = G_list[imautosome, ]

y <- DGEList(count_matrix, genes=G_list, group=data$Diagnosis)

DX <- factor(data$Diagnosis)
Region <- factor(data$Region)
design <- model.matrix(~Region + DX)

keep <- filterByExpr(y)  # doing it based on group
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)
fit <- glmQLFit(y, design)
```

![](images/2020-06-02-08-14-52.png)

To detect genes that are differentially expressed between Case and Control,
adjusting for brain region differences:

```r
qlf <- glmQLFTest(fit, coef=3)
topTags(qlf)
```

![](images/2020-06-02-08-18-39.png)

I somewhat remember this result... I'll take another look, but I think it was
something related to these genes being all zero in one of the brain regions?

![](images/2020-06-05-07-00-19.png)

That was the top gene... I think they're all like that?

![](images/2020-06-05-07-02-40.png)

Well, not all of them, but a good amount of the 10 that are below .05 (and .1).

## 3.4.3

Here we'l be correcting for batch effects as well, which we know to exist based
on our old plots. We could just try COMBAT and then the analysis above again.
But let's see how edgeR would do it.

The challenge here is that we have batching and region within batches. So, it'd
be more like example 3.5? This might be way too confusing... why not jst adjust
for batch effects with COMBAT, which seems to be working, and use the model
above?

```r
library(sva)
library(edgeR)
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ]  # removing ACC outlier
rownames(data) = data$submitted_name  # just to ensure compatibility later

grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
# remove that weird .num after ENSG
id_num = sapply(grex_vars,
                function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num

dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]
library('biomaRt')
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]
imautosome = which(G_list$chromosome_name != 'X' &
                   G_list$chromosome_name != 'Y' &
                   G_list$chromosome_name != 'MT')
count_matrix = count_matrix[imautosome, ]
G_list = G_list[imautosome, ]

x <- DGEList(count_matrix, genes=G_list, group=data$Diagnosis)
batch = factor(data$run_date)
covar_mat = cbind(data$Diagnosis, data$Region)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=NULL,
                              covar_mod=covar_mat)

y <- DGEList(adjusted_counts, genes=G_list, group=data$Diagnosis)

DX <- factor(data$Diagnosis)
Region <- factor(data$Region)
design <- model.matrix(~Region + DX)

keep <- filterByExpr(y)  # doing it based on group
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

y <- estimateDisp(y, design, robust=TRUE)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=3)
topTags(qlf)
```

The results are similar to the ones from before:

![](images/2020-06-05-07-10-19.png)

Except that now, after the batch adjustments, our counts don't look that off:

![](images/2020-06-05-07-13-10.png)

What happens if we run a model with interactions?

```r
design <- model.matrix(~Region + DX + Region:DX)
y2 <- estimateDisp(y, design, robust=TRUE)
fit2 <- glmQLFit(y2, design)
qlf2 <- glmQLFTest(fit2, coef=3)
topTags(qlf2)
```

![](images/2020-06-05-07-15-58.png)

In that case my results mostly go away. Is there anything good in the
interaction?

![](images/2020-06-05-07-18-14.png)

Let's take a closer look at that NEURO gene...

![](images/2020-06-05-07-21-12.png)

Yeah, probably not expressed in one of the regions. Let me focus on this result
a bit more:

```r
qlf3 = glmQLFTest(fit2, coef=4)
summary(decideTests(qlf3))
plotMD(qlf3)
```

![](images/2020-06-05-07-30-16.png)

Maybe these are driven by outlier data? In any case, let's see how it looks when
we narrow it down to some meaningful log fold change:

```r
tr <- glmTreat(fit2, coef=4, lfc=log2(1.2))
topTags(tr)
```

That didn't help much. So, let's see what's going on with that gene:

```r
g = 'ENSG00000164600'
library(ggplot2)
lcpm = cpm(y2, log=TRUE)
plot_data = data.frame(x=1:ncol(y2), y=lcpm[g, ], region=data$Region,
                       group=data$Diagnosis)
ggplot(plot_data, aes(x=x, y=y, shape=region, color=group)) + geom_point()
```

![](images/2020-06-05-07-40-11.png)

We also have that one outlier there... not sure if we can recode it as ACC, or
is it just a bad Caudate sample?

Also, we need to remember that the coefficient DXControl represents the
difference in mean expression between Control and the reference level (Case),
for the ACC, which is the reference level for Region.

Similarly, the coefficient RegionCaudate represents the difference in mean
expression between Caudate and ACC, for the Cases (reference in DX).

The interaction term asks if the change in expression between case and control
is the same for ACC as it is for Caudate. 

In such a small dataset, outliers like these will make a difference. Let's plot
these data again. First, before COMBAT:

```r
library(pca3d)
library(caret)
set.seed(42)
# remove genes with zero or near zero variance so we can run PCA
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
pca <- prcomp(X, scale=TRUE)
pca3d(pca, group=data$Diagnosis)
```

![](images/2020-06-05-07-58-28.png)

There are the two clear groups there, and a couple outliers. 

Let's focus on the 2 dimensions first:

```r
pca2d(pca, group=data$run_date, shape=as.numeric(data$Region))
```

![](images/2020-06-05-08-03-12.png)

Both batch and Region are very evident in the first 2 PCs. So, let's run COMBAT
(on this cleaned up data) and see what remains:

```r
library(sva)
batch = factor(data$run_date)
covar_mat = cbind(data$Diagnosis, data$Region)
adjusted_counts <- ComBat_seq(t(X), batch=batch, group=NULL,
                              covar_mod=covar_mat)
pca2 <- prcomp(t(adjusted_counts), scale=TRUE)
pca2d(pca2, group=data$run_date, shape=as.numeric(data$Region))
```

![](images/2020-06-05-08-09-43.png)

Hum... nothing happened.

```r
adjusted_counts2 <- ComBat_seq(t(X), batch=batch, group=data$Diagnosis)
pca3 <- prcomp(t(adjusted_counts2), scale=TRUE)
dev.new()
pca2d(pca3, group=data$run_date, shape=as.numeric(data$Region))
```

![](images/2020-06-09-08-12-43.png)

It's still not taking care of it... colors are clearly split in PC1.

Is it brain bank?

```r
dev.new()
pca2d(pca3, group=data$bainbank, shape=as.numeric(data$Region))
```

![](images/2020-06-09-08-24-54.png)

No, batch colors it better. How about using that one package for QC?

https://bmcresnotes.biomedcentral.com/articles/10.1186/s13104-019-4179-2

https://www.bioconductor.org/packages/release/bioc/html/BatchQC.html

# 2020-06-09 08:03:10

What if we only keep the group variable?


## 4.4

This analysis doesn't look that different from the other ones...


(possible alternatives to edgeR?)
https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0067019
https://www.datanovia.com/en/lessons/mixed-anova-in-r/ 
