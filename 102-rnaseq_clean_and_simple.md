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
                myregion, abs(SDT))
saveRDS(X2, file=fname)
```

## Simple models

Now that the data is cleaner, let's run some simple models:

```r
library(caret)
library(doParallel)

ncores = 32
myseed = 42
clf = 'gaussprLinear'
registerDoParallel(ncores)
getDoParWorkers()

# change later to clean data!
myregion='ACC'
just_target = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
if (myregion == 'both') {
    fname = '~/data/rnaseq_derek/X_noPH_zv_nzv_center_scale_SDT5_normal.rds'
    y = just_target[, 'Diagnosis']
} else {
    fname = sprintf('~/data/rnaseq_derek/X_%snoPH_zv_nzv_center_scale_SDT5_normal.rds',
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

slda failed due to protection for stack overflow. How about treebag? Same
error... bayesglm too. That's becomeing a big problem. I might need to go
straight for regularized models, or at least with some sort of FS. even
glmStepAIC got protected. Oh well. Let me try doing some tuning with the entire
dataset first. Going for logistic regression and cranking up lambda. Maybe there's a nice interplay of alpha and lambda here?

```r
nfolds = 5
nreps = 10
clf = 'glmnet'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(lambda=c(1, 5, 10, 50, 100),
                     alpha=c(0, .25, .75, 1))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

```
  alpha  lambda  ROC        Sens   Spec     
  0.00     1     0.6552381  0.516  0.7033333
  0.00     5     0.6533333  0.520  0.7033333
  0.00    10     0.6582857  0.508  0.7000000
  0.00    50     0.6722857  0.504  0.7066667
  0.00   100     0.6717143  0.492  0.7104762
```

How about SVM?

```r
nfolds = 5
nreps = 10
clf = 'svmLinear'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(C=c(.001, .1, 1, 5, 10, 100, 200, 500, 1000))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

```
  C      ROC        Sens   Spec     
  1e-03  0.6060000  0.384  0.7714286
  1e-01  0.6206667  0.368  0.7704762
  1e+00  0.5886667  0.288  0.7909524
  5e+00  0.6326667  0.380  0.7738095
  1e+01  0.5993333  0.336  0.7871429
  1e+02  0.5783810  0.300  0.8033333
  2e+02  0.6010476  0.328  0.7942857
  5e+02  0.6170476  0.356  0.7485714
  1e+03  0.5806667  0.292  0.7947619
```

I could also try some veeery shallow trees? For example, rf fixing mtry, or
'RRFglobal' fixing mtry but learning the regularization?

```r
nfolds = 5
nreps = 10
clf = 'rf'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(mtry=c(3, 5, 10, 15, 20))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

```
  mtry  ROC        Sens   Spec     
   3    0.6586190  0.396  0.8085714
   5    0.6491905  0.436  0.7552381
  10    0.6499524  0.444  0.7661905
  15    0.6360000  0.428  0.7423810
  20    0.6446667  0.464  0.7695238
```

```r
nfolds = 5
nreps = 10
clf = 'RRFglobal'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(mtry=c(1, 3, 5, 10, 15, 20, 30),
                     coefReg=c(.01, .1, 1))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

Another option is to try some other bagged models: AdaBag, LogicBag, cforest...
logicBag had lots of errors. AdaBag also failed... can they not handle
probabilities?


```r
nfolds = 5
nreps = 10
clf = 'cforest'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(mtry=c(1, 2, 3, 5, 10, 15, 20))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

```r
nfolds = 5
nreps = 10
clf = 'AdaBag'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(mfinal=c(1, 2, 3, 5, 10, 15, 20),
                     maxdepth=c(1, 2, 3))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

I should also play with some of the Naive Bayes variants: manb, naive_bayes,
nbDiscrete. But manb had lots of errors... nbDiscrete also didn't run...

```r
nfolds = 5
nreps = 10
clf = 'naive_bayes'
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(usekernel=c(FALSE, TRUE),
                     laplace=c(0, .5, 1),
                     adjust=c(0, .1, 1, 10))
set.seed(42)
fit <- train(X, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

```
  usekernel  ROC        Sens   Spec     
  FALSE      0.5804286  0.752  0.4028571
   TRUE      0.5562381  0.580  0.5180952

Tuning parameter 'laplace' was held constant at a value of 0
Tuning parameter 'adjust' was held constant at a value of 1
```

## Univariate filter

Why not use some univariate filtering?

```r
nfolds = 5
nreps = 10

set.seed(42)
filterCtrl <- sbfControl(functions = rfSBF, method = "repeatedcv",
                         repeats = nreps, number=nfolds)
set.seed(42)
rfWithFilter <- sbf(X, y, sbfControl = filterCtrl)
rfWithFilter
```

# TODO
* try
  http://topepo.github.io/caret/miscellaneous-model-functions.html#partial-least-squares-discriminant-analysis
  ?
* crank up MBO