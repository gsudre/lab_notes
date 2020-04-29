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

# TODO
* add covariates... do they help?
* play with Caudate
