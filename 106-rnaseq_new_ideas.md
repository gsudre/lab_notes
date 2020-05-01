# 2020-04-29 08:10:01

Let's implement some of the ideas from note 104 TODO list. First, let's grab
Derek's original data to make sure there are no misalignments, and we keep the
ENS number to later map it to genes.

```r
df = read.csv('~/data/rnaseq_derek/UPDATED_file_for_derek_add_cause_of_death.csv')
df = df[!duplicated(df$submitted_name),]
data = read.table('~/data/rnaseq_derek/logCPM_rnaseq.txt')
data = t(data)
sn = gsub(x=rownames(data), pattern='X', replacement='')
data = cbind(as.numeric(sn), data)
colnames(data)[1] = 'submitted_name'
m = merge(df, data, by='submitted_name', all.x=F, all.y=T)
pop_code = read.csv('~/data/rnaseq_derek/file_pop.csv')
m2 = merge(m, pop_code, by='hbcc_brain_id')
pcs = read.table('~/data/rnaseq_derek/HM3_b37mds.mds', header=1)
myids = sapply(1:nrow(pcs), function(x) as.numeric(gsub('BR', '',
                                                        strsplit(as.character(pcs[x,'IID']), '_')[[1]][1])))
pcs$numids = myids
m3 = merge(m2, pcs, by.x='hbcc_brain_id', by.y='numids', all.x=T, all.y=F)
```

Now we go ahead with the binomial idea for gene filtering. Out of the 119
observations, 56 are ACC and 63 are caudate. Actually, there are 56 unique ACC
observations, and only 58 for caudate. So, let's choose which ones in the
caudate to take, for the repeated subjects. Looking at the UPDATED file, it's
unclear how to choose it. It could be done based on pcnt_optical_duplicates or
clusters. DEFINITELY A QUESTION FOR DEREK.

I could also just check the number of non-expressed genes in each of them, and
just pick the one with the fewest. Or, just choose based on the ACC measurement
that's being used for that subject?

```r
library(caret)
cdata = m3[m3$Region=='Caudate',]
grex_names = colnames(m3)[grepl(colnames(m3), pattern='^ENS')]
pp = preProcess(cdata[, grex_names],
                method=c('zv', 'nzv', 'range'), rangeBounds=c(0,1))
a = predict(pp, cdata[, grex_names])
n0 = rowSums(a==0)
idx = which(cdata$hbcc_brain_id==2877)
idx

[1] 25 26 27 28 29 30

n0[idx]
  49   50   51   52   53   54 
3980 3712 3599 3794 4081 2347 
```

So, row 54 has the fewest genes with zero transcription counts. We'll take that.

```r
cdata = cdata[-c(25:29, ]
m4 = rbind(cdata, m3[m3$Region=='ACC',])
saveRDS(m4, file='~/data/rnaseq_derek/complete_data_04292020.rds')
```

So, under a binomial distribution with 56 tries, how many should I get positive
to consider that not happening by chance?

```r
> x=38
> binom.test(c(x, 56-x), p = .5)

	Exact binomial test

data:  c(x, 56 - x)
number of successes = 38, number of trials = 56, p-value = 0.01045
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.5403638 0.7971455
sample estimates:
probability of success 
             0.6785714 

> x=39
> binom.test(c(x, 56-x), p = .5)

	Exact binomial test

data:  c(x, 56 - x)
number of successes = 39, number of trials = 56, p-value = 0.004562
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.5590326 0.8122013
sample estimates:
probability of success 
             0.6964286 
```

So, for ACC (56 subjects / coin tosses), I need at least 39 subjects with
success (transcription reads above 0), for that to be considered a expressed
gene with more than chance probability (.5) at p<.01.

For Caudate, our number is 58 tries instead:

```r
> x=39
> binom.test(c(x, 58-x), p = .5)

	Exact binomial test

data:  c(x, 58 - x)
number of successes = 39, number of trials = 58, p-value = 0.01193
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.5365938 0.7899462
sample estimates:
probability of success 
             0.6724138 

> x=40
> binom.test(c(x, 58-x), p = .5)

	Exact binomial test

data:  c(x, 58 - x)
number of successes = 40, number of trials = 58, p-value = 0.005355
alternative hypothesis: true probability of success is not equal to 0.5
95 percent confidence interval:
 0.5545582 0.8046135
sample estimates:
probability of success 
             0.6896552 
```

So we'll go with 40 for the Caudate. Let's also make some PCA plots just to
have some idea of the number of components to keep across dimensionality
reduction methods.

```r
library(caret)
data = readRDS('~/data/rnaseq_derek/complete_data_04292020.rds')
data = data[data$Region=='ACC',]
grex_names = colnames(data)[grepl(colnames(data), pattern='^ENS')]
pp = preProcess(data[, grex_names],
                method=c('zv', 'nzv', 'range'), rangeBounds=c(0,1))
a = predict(pp, data[, grex_names])
n0 = colSums(a==0)
imbad = names(n0)[n0> 15]  # ACC
imbad = names(n0)[n0> 18]  # Caudate
good_grex = grex_names[!(grex_names %in% imbad)]
res.pca <- prcomp(data[, good_grex])
library(factoextra)
fviz_eig(res.pca, ncp=40)

saveRDS(data[, good_grex],
        file='~/data/rnaseq_derek/goodgrexACC_binp01_04292020.rds')
```

![](images/2020-04-29-14-50-45.png)

It looks like 10 is a nice round number of PCs here. We could go nuts on this
analysis, but I just want to have a general idea. For completeness, this is the
plot for the Caudate:

![](images/2020-04-29-14-53-18.png)

So, 10 should be alright here as well.

Let me run a series of dimensionality reduction methods and check if anything
looks good.

I got an interesting plots with t-SNE, but the biggest issue is that it doesn't
provide any sort of variable weights, so it'd be hard to assgign back importance
to each gene. But there are other methods that could do that, so the idea itself
is not bad.

So, I started playing with this package:
https://cran.r-project.org/web/packages/dimRed/vignettes/dimensionality-reduction.pdf

And I could see right away hat there was something funky with the first PC:

![](images/2020-04-29-20-33-09.png)

and I got that both using PCA and Isomaps. So, before I try getting fancy with
the methods, I need to see if there isn't something funky witht the data that's
carrying all this variance.

But the code I used was:

```r
data = readRDS('~/data/rnaseq_derek/goodgrexACC_binp01_04292020.rds')
library(Rtsne)
X = normalize_input(as.matrix(data))
library(dimRed)
data_emb <- lapply(embed_methods, function(x) embed(X, x))
names(data_emb) <- embed_methods
plot_R_NX(data_emb)
plot(data_emb[[2]], type='2vars')
```

The R_NX plot was very similar for Isomap and PCA. If I were to optimize the
parameter for Isomap, then this would work:

```r
kk <- floor(seq(5, 40, length.out = 20))
emb <- lapply(kk, function(x) embed(X, "Isomap", knn = x))
qual <- sapply(emb, function(x) quality(x, "Q_local"))
ind_max <- which.max(qual)
k_max <- kk[ind_max]
```

which gives me about 23. And we could run all available models this way:

```r
embed_methods <- dimRedMethodList()
embed_methods = embed_methods[-c(1, 4, 11)] # remove AutoEncoder, ICA, PCAL1 
quality_methods <- c("Q_local", "Q_global", "AUC_lnK_R_NX",
                     "cophenetic_correlation")
quality_results <- matrix(
  NA, length(embed_methods), length(quality_methods),
  dimnames = list(embed_methods, quality_methods)
)
embedded_data <- list()
for (em in embed_methods) {
    print(sprintf('Trying %s', em))
    embedded_data[[em]] <- embed(X, em)
    for (q in quality_methods) {
        print(sprintf('%s (%s)', em, q))
        try(quality_results[em, q] <- quality(embedded_data[[em]], q))
    }
}
```

AutoEncoder needed tensorflow and I didn't want to install in my MacAir. ICA was
taking forever. PCA_L1 did nt compile. NNMF requires matrix with only non-negative
entries! UMAP was having issues installing in the MacbookAir... maybe all of
these will be options in the cluster. We'll see. For now, we need to optimize
the paremeters anyways. So

But of course we should optimize eahc one individually before running this! Each
model has its own set of parameters, and it maybe the same parameter maximizes
one quality metric but not others? Ideally, I could keep the number of
dimensions constant, and make bar graphs for varying Ks. Each group of bars
would be a K, and four bars per group, one per metric. I could add more
dimensions later, but I think staring with 2 is a good one for now, and we can
run some classifiers or visualizations on that. If I script it, it'll be trivial
to increase the number of dimensions later.

So, before we go any further, let's see if that weird PC1 by PC2 plot was
indicative of something I should pay attention to...

# 2020-05-01 08:10:20

So, let's start with Isomaps, because we can optimize knn and dimensions.

```r
data = readRDS('~/data/rnaseq_derek/goodgrexACC_binp01_04292020.rds')
library(Rtsne)
X = normalize_input(as.matrix(data))
library(dimRed)
nd = 2
kk <- floor(seq(5, 40, length.out = 20))
emb <- lapply(kk, function(x) embed(X, "Isomap", knn = x, ndim= nd))
names(emb) = kk
quality_methods <- c("Q_local", "Q_global", "AUC_lnK_R_NX",
                     "cophenetic_correlation")
quality_results <- matrix(
  NA, length(kk), length(quality_methods),
  dimnames = list(kk, quality_methods)
)
for (k in kk) {
    print(sprintf('Evaluating k=%d',k))
    kname = as.character(k)
    for (q in quality_methods) {
        try(quality_results[kname, q] <- quality(emb[[kname]], q))
    }
}

#plotting
library(ggplot2)
nneighbors = c()
qual_metric = c()
value = c()
for (i in 1:nrow(quality_results)) {
    for (j in 1:ncol(quality_results)) {
        nneighbors = c(nneighbors, as.numeric(rownames(quality_results)[i]))
        qual_metric = c(qual_metric, colnames(quality_results)[j])
        value = c(value, quality_results[i, j])
    }
}
plot_data <- data.frame(qual_metric,nneighbors,value)
ggplot(plot_data, aes(fill=qual_metric, y=value, x=nneighbors)) +
    geom_bar(position="dodge", stat="identity")
```

# TODO

* keep it to only metrics that give an inverse transform! (i.e. no tsne)
* rank average the different metrics to choose best combination, including
  number of dimensions
* add covariates... do they help?
* play with Caudate
