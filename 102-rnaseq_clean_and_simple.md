# 2020-04-11 13:28:09

Since the results with xgbTree didn't generalize, let's go back to the drawing
board. Let's first check how much cleaning we need to do to remove outliers. I
definitely don't want to lose subjects, so I'll check variable by variable. 

I'll take advantage that I have already scaled the data, so I can take a look at
the Z range for outliers.

```r
myregion='ACC'
fname = sprintf('~/data/rnaseq_derek/X_%snoPH_zv_nzv_center_scale.rds',
                    myregion)
X = readRDS(fname)
# at this point everything I have in X is numeric, so it's safe to proceed with this
m = apply(X, 2, max)
# setting anything above SDT standard deviations to SDT. In ACC I only had 5% of
# variables above 4, 1% above 5
SDT = 5
trim_vars = names(which(m > SDT))
for (v in trim_vars) {
    X[X[, v] > SDT, v] = SDT
}
# do the same thing below zero
m = apply(X, 2, min)
SDT = -SDT
trim_vars = names(which(m < SDT))
for (v in trim_vars) {
    X[X[, v] < SDT, v] = SDT
}
```

Now that I trimmed some of the variables, let's make sure they are as normal as
they can be:

```r
library(bestNormalize)
for (v in 1:ncol(X[, 1:20])) {
    if ((v %% 100)==0) {
        print(sprintf('%d / %d', v, ncol(X)))
    }
    bn = bestNormalize(X[, v], warn=F, k=10, r=10)
    X[, v] = bn$x.t
}
# finally, re-scale everything to make sure out mean and SD are the same across 
# variables, in case some classifiers need that (looking at you, LDA)
X2 = scale(X, center=T, scale=T)
fname = sprintf('~/data/rnaseq_derek/X_%snoPH_zv_nzv_center_scale_SDT%d_normal.rds',
                myregion, SDT)
writeRDS(X2, fname=fname)
```

## Simple models

Now that the data is cleaner, let's run some simple models:

```r
library(caret)
library(caretEnsemble)
library(doParallel)

ncores = 2
myseed = 42
clf = 'slda'
registerDoParallel(ncores)
getDoParWorkers()

# change later to clean data!
myregion='ACC'
just_target = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
if (myregion == 'both') {
    fname = '~/data/rnaseq_derek/X_noPH_zv_nzv_center_scale.rds'
    y = just_target[, 'Diagnosis']
} else {
    fname = sprintf('~/data/rnaseq_derek/X_%snoPH_zv_nzv_center_scale.rds',
                    myregion)
    y = just_target[just_target$Region==myregion, 'Diagnosis']
}
X = readRDS(fname)

fitControl <- trainControl(method = "none",
                           allowParallel = TRUE,
                           classProbs = TRUE)

y_probs = c()
y_preds = c()
best_params = c()
for (test_row in 1:nrow(X)) {
    train_rows = setdiff(1:nrow(X), test_row)
    X_train <- X[train_rows, ]
    X_test <- X[-train_rows, ]
    y_train <- y[train_rows]
    y_test <- y[-train_rows]

    print(sprintf('LOOCV %d / %d', test_row, nrow(X)))
    set.seed(myseed)
    fit <- train(X_train, y_train,
                 trControl = fitControl,
                 method = clf)
    # updated LOOCV predictions
    y_probs = rbind(y_probs, predict(fit, X_test, type='prob'))
    y_preds = c(y_preds, levels(y)[predict(fit, X_test)])

    # just some ongoing summary
    dat = cbind(data.frame(obs = y[1:nrow(y_probs)],
                           pred = factor(y_preds, levels=levels(y))),
                y_probs)
    print(twoClassSummary(dat, lev=levels(y)))
}
```
