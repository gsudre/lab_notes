# 2019-01-23 15:04:55

The more I played with the logistic regression results, the less I liked them.
The comparison to NV doesn't make that much sense to me.

So, let's take a couple hours to go back to the prediction framework, this time
using a few simpler algorithms.

```r
if (Sys.info()['sysname'] == 'Darwin') {
  max_mem = '16G'
  base_name = '~/data/'
} else {
  max_mem = paste(Sys.getenv('SLURM_MEM_PER_NODE'),'m',sep='')
  base_name = '/data/NCR_SBRB/'
}
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/struct_thickness_11142018_260timeDiff12mo.RData.gz', base_name))
df = merge(clin, data, by='MRN')
df$OLS_inatt_categ = NULL
df[df$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
df[df$OLS_inatt_slope > -.33 & df$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
df[df$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
df = df[df$DX != 'NV', ]
df$OLS_inatt_categ = as.factor(df$OLS_inatt_categ)
x = colnames(df)[grepl(pattern = '^v_rh', colnames(df))]

library(caret)
set.seed(42)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the models
set.seed(42)
modelLinearSvm <- train(df[, x], df$OLS_inatt_categ, method="svmLinear", trControl=control)
set.seed(42)
modelNB <- train(df[, x], df$OLS_inatt_categ, method="nb", trControl=control)


```

Philip doens't want to go this route yet.

# 2019-01-25 14:05:00

OK, but given that the massive univariate (MU) results are not doing great for
classification, it's worth trying a MVPA approach, even if just for comparison.
Let's do it for structural first, as it runs quite fast.

So, we start with scaling the voxels. Then, we can do a PCA, Gaussianize, and
use classifiers (GNB, SVM, Tree), or just use the latter 2 after scaling (not
fair to use GNB with so many correlated variables).

Let's create a quick script for that:

```r
if (Sys.info()['sysname'] == 'Darwin') {
  max_mem = '16G'
  base_name = '~/data/'
} else {
  max_mem = paste(Sys.getenv('SLURM_MEM_PER_NODE'),'m',sep='')
  base_name = '/data/NCR_SBRB/'
}
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/struct_thickness_11142018_260timeDiff12mo.RData.gz', base_name))

# first cleanup
# remove constant variables that screw up PCA and univariate tests
print(sprintf('Starting with %d variables', ncol(data)-1))
print('Removing variables with near zero or zero variance')
nNAs = colSums(is.na(data))  # number of NAs in each variable
# remove variables that are all NAs
data = data[, nNAs < nrow(data)]
library(caret)
nzv <- nearZeroVar(data)
data = data[, -nzv]
print(sprintf('Features remaining: %d (%d with NAs)', ncol(data)-1, sum(nNAs>0)))

df = merge(clin, data, by='MRN')
df$OLS_inatt_categ = NULL
df[df$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
df[df$OLS_inatt_slope > -.33 & df$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
df[df$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
df = df[df$DX != 'NV', ]
df$OLS_inatt_categ = as.factor(df$OLS_inatt_categ)
x = colnames(df)[grepl(pattern = '^v', colnames(df))]

pp <- preProcess(df[, x], method = c("center", "scale"))
df[, x] = predict(pp, df[, x])


set.seed(42)
control <- trainControl(method="repeatedcv", number=10, repeats=10)
# train the models
set.seed(42)
grid <- expand.grid(C = c(10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1))
modelLinearSvm <- train(df[, x], df$OLS_inatt_categ, method="svmLinear",
                        trControl=control, tuneGrid=grid)

grid <- expand.grid(k = seq(10, 40, 5))
modelKNN <- train(df[, x], df$OLS_inatt_categ, method="knn", trControl=control,
                  tuneGrid=grid)


# if using PCA
preproc='kaiser'
print('Running PCA')
# quick hack to use na.action on prcomp
fm_str = sprintf('~ %s', paste0(x, collapse='+ ', sep=' '))
pca = prcomp(as.formula(fm_str), df[, x], scale=T, na.action=na.exclude)
eigs <- pca$sdev^2
if (grepl(pattern='elbow', preproc)) {
    library(nFactors)
    nS = nScree(x=eigs)
    keep_me = 1:nS$Components$noc
} else if (grepl(pattern='kaiser', preproc)) {
    library(nFactors)
    nS = nScree(x=eigs)
    keep_me = 1:nS$Components$nkaiser
} else {
    keep_me = 1:nrow(df)
}
df2 = data.frame(pca$x[, keep_me])
df2[, 'MRN'] = df$MRN
df2[, 'OLS_inatt_categ'] = df$OLS_inatt_categ
df = df2
x = colnames(df)[grepl(pattern = '^PC', colnames(df))]
pp <- preProcess(df[, x], method = c("center", "scale", "YeoJohnson"))
df[, x] = predict(pp, df[, x])

set.seed(42)
grid <- expand.grid(C = c(10^-4, 10^-3, 10^-2, 10^-1, 1, 10, 100))
modelNB <- train(df[, x], df$OLS_inatt_categ, method="nb", trControl=control)
```

Oh well, still nothing. The best I can do here is .45 accuracy, and although
better than chance (.33), it's not great. 

But what if the question were: will my kid get better or not?

```r
if (Sys.info()['sysname'] == 'Darwin') {
  max_mem = '16G'
  base_name = '~/data/'
} else {
  max_mem = paste(Sys.getenv('SLURM_MEM_PER_NODE'),'m',sep='')
  base_name = '/data/NCR_SBRB/'
}
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/struct_thickness_11142018_260timeDiff12mo.RData.gz', base_name))

# first cleanup
# remove constant variables that screw up PCA and univariate tests
print(sprintf('Starting with %d variables', ncol(data)-1))
print('Removing variables with near zero or zero variance')
nNAs = colSums(is.na(data))  # number of NAs in each variable
# remove variables that are all NAs
data = data[, nNAs < nrow(data)]
library(caret)
nzv <- nearZeroVar(data)
data = data[, -nzv]
print(sprintf('Features remaining: %d (%d with NAs)', ncol(data)-1, sum(nNAs>0)))

df = merge(clin, data, by='MRN')
df$OLS_inatt_categ = NULL
df[df$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
df[df$OLS_inatt_slope > -.33 & df$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
df[df$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
df = df[df$DX != 'NV', ]
df[df$OLS_inatt_categ=='marked', 'OLS_inatt_categ'] = 'improve'
df[df$OLS_inatt_categ=='mild', 'OLS_inatt_categ'] = 'improve'
df$OLS_inatt_categ = as.factor(df$OLS_inatt_categ)
x = colnames(df)[grepl(pattern = '^v', colnames(df))]

pp <- preProcess(df[, x], method = c("center", "scale"))
df[, x] = predict(pp, df[, x])


set.seed(42)
control <- trainControl(method="repeatedcv", number=10, repeats=10)
# train the models
set.seed(42)
grid <- expand.grid(C = c(10^-6, 10^-5, 10^-4, 10^-3, 10^-2, 10^-1, 1))
modelLinearSvm <- train(df[, x], df$OLS_inatt_categ, method="svmLinear",
                        trControl=control, tuneGrid=grid)

grid <- expand.grid(k = seq(10, 40, 5))
modelKNN <- train(df[, x], df$OLS_inatt_categ, method="knn", trControl=control,
                  tuneGrid=grid)


# if using PCA
preproc='kaiser'
print('Running PCA')
# quick hack to use na.action on prcomp
fm_str = sprintf('~ %s', paste0(x, collapse='+ ', sep=' '))
pca = prcomp(as.formula(fm_str), df[, x], scale=T, na.action=na.exclude)
eigs <- pca$sdev^2
if (grepl(pattern='elbow', preproc)) {
    library(nFactors)
    nS = nScree(x=eigs)
    keep_me = 1:nS$Components$noc
} else if (grepl(pattern='kaiser', preproc)) {
    library(nFactors)
    nS = nScree(x=eigs)
    keep_me = 1:nS$Components$nkaiser
} else {
    keep_me = 1:nrow(df)
}
df2 = data.frame(pca$x[, keep_me])
df2[, 'MRN'] = df$MRN
df2[, 'OLS_inatt_categ'] = df$OLS_inatt_categ
df = df2
x = colnames(df)[grepl(pattern = '^PC', colnames(df))]
pp <- preProcess(df[, x], method = c("center", "scale", "YeoJohnson"))
df[, x] = predict(pp, df[, x])

set.seed(42)
modelNB <- train(df[, x], df$OLS_inatt_categ, method="nb", trControl=control)
```

Still, the best I could do was .56 :(