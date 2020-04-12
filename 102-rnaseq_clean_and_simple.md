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

```
  mtry  coefReg  ROC        Sens   Spec     
   1    0.01     0.5501905  0.464  0.5547619
   1    0.10     0.4676667  0.416  0.5261905
   1    1.00     0.6023810  0.440  0.6985714
   3    0.01     0.5649048  0.444  0.6480952
   3    0.10     0.5070476  0.436  0.5528571
   3    1.00     0.6261429  0.448  0.6676190
   5    0.01     0.5096190  0.428  0.5795238
   5    0.10     0.5208571  0.468  0.5428571
   5    1.00     0.6324286  0.436  0.7052381
  10    0.01     0.4944762  0.436  0.5447619
  10    0.10     0.5433333  0.452  0.6214286
  10    1.00     0.6126667  0.484  0.7052381
  15    0.01     0.5173810  0.456  0.5852381
  15    0.10     0.5609048  0.456  0.5952381
  15    1.00     0.6236190  0.456  0.7004762
  20    0.01     0.5276190  0.468  0.6014286
  20    0.10     0.5375714  0.476  0.5752381
  20    1.00     0.6063333  0.456  0.6909524
  30    0.01     0.5244762  0.460  0.5433333
  30    0.10     0.5099524  0.452  0.5819048
  30    1.00     0.6055238  0.476  0.6500000

ROC was used to select the optimal model using the largest value.
The final values used for the model were mtry = 5 and coefReg = 1.
```

Another option is to try some other bagged models: AdaBag, LogicBag, cforest...
logicBag had lots of errors. AdaBag also failed... can they not handle
probabilities? cforest also failed.

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
  usekernel  laplace  adjust  ROC        Sens   Spec     
  FALSE      0.0       0.0    0.5804286  0.752  0.4028571
  FALSE      0.0       0.1    0.5804286  0.752  0.4028571
  FALSE      0.0       1.0    0.5804286  0.752  0.4028571
  FALSE      0.0      10.0    0.5804286  0.752  0.4028571
  FALSE      0.5       0.0    0.5804286  0.752  0.4028571
  FALSE      0.5       0.1    0.5804286  0.752  0.4028571
  FALSE      0.5       1.0    0.5804286  0.752  0.4028571
  FALSE      0.5      10.0    0.5804286  0.752  0.4028571
  FALSE      1.0       0.0    0.5804286  0.752  0.4028571
  FALSE      1.0       0.1    0.5804286  0.752  0.4028571
  FALSE      1.0       1.0    0.5804286  0.752  0.4028571
  FALSE      1.0      10.0    0.5804286  0.752  0.4028571
   TRUE      0.0       0.0          NaN    NaN        NaN
   TRUE      0.0       0.1    0.5110000  0.060  0.9704762
   TRUE      0.0       1.0    0.5562381  0.580  0.5180952
   TRUE      0.0      10.0    0.5683333  0.728  0.3352381
   TRUE      0.5       0.0          NaN    NaN        NaN
   TRUE      0.5       0.1    0.5110000  0.060  0.9704762
   TRUE      0.5       1.0    0.5562381  0.580  0.5180952
   TRUE      0.5      10.0    0.5683333  0.728  0.3352381
   TRUE      1.0       0.0          NaN    NaN        NaN
   TRUE      1.0       0.1    0.5110000  0.060  0.9704762
   TRUE      1.0       1.0    0.5562381  0.580  0.5180952
   TRUE      1.0      10.0    0.5683333  0.728  0.3352381
```

So, there's plenty of room for improvement by tuning values in most models. Maybe MBO?

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

```
Selection By Filter

Outer resampling method: Cross-Validated (5 fold, repeated 10 times)

Resampling performance:

Accuracy Kappa AccuracySD KappaSD
0.553     0    0.01531       0

Using the training set, 3533 variables were selected:
grex2, grex13, grex18, grex20, grex40...

During resampling, no variables were selected.
```

This didn't do too well, but I'd be more comfortable tuning this by hand
anyways. For example, something like this:

```r
p_thresh = .05
ps = sapply(1:ncol(X), function(v) t.test(X[y=='Case', v],
                                          X[y=='Control', v])$p.value)
X2 = X[, ps < p_thresh]
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
fit <- train(X2, y,
                 trControl = fitControl,
                 method = clf,
                 tuneGrid=mygrid,
                 metric='ROC')
```

```
  mtry  ROC        Sens   Spec     
   3    0.8609524  0.676  0.8347619
   5    0.8372381  0.668  0.8280952
  10    0.8360952  0.656  0.8476190
  15    0.8399048  0.676  0.8338095
  20    0.8364762  0.660  0.8504762
```

This is a considerable improvement over using the entire dataset! Let's crank up
the pvalue to .01 and do some very shallow trees:

```
  mtry  ROC        Sens   Spec     
   1    0.9053333  0.720  0.8661905
   2    0.8906190  0.728  0.8604762
   3    0.8896190  0.716  0.8676190
   5    0.8839524  0.708  0.8442857
  10    0.8888095  0.724  0.8609524
  15    0.8838095  0.712  0.8538095
  20    0.8836667  0.708  0.8471429

ROC was used to select the optimal model using the largest value.
The final value used for the model was mtry = 1.
> dim(X2)
[1]  56 911
```

This is the best result we've had so far. How does it perform in a LOOCV?

```r
fitControl <- trainControl(method = "none",
                           allowParallel = TRUE,
                           classProbs = TRUE)
p_thresh = .01
clf = 'rf'

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

    ps = sapply(1:ncol(X), function(v) t.test(X_train[y_train=='Case', v],
                                              X_train[y_train=='Control', v])$p.value)
    good_vars = ps < p_thresh
    X_train2 = X_train[, good_vars]
    X_test2 = t(as.matrix(t(as.matrix(X_test))[, good_vars]))
    mygrid = data.frame(mtry=1)

    set.seed(42)
    fit <- train(X_train2, y_train,
                    trControl = fitControl,
                    method = clf,
                    tuneGrid=mygrid,
                    metric='ROC')
    # updated LOOCV predictions
    y_probs = rbind(y_probs, predict(fit, X_test2, type='prob'))
    y_preds = c(y_preds, levels(y)[predict(fit, X_test2)])

    # just some ongoing summary
    dat = cbind(data.frame(obs = y[1:nrow(y_probs)],
                           pred = factor(y_preds, levels=levels(y))),
                y_probs)
    print(twoClassSummary(dat, lev=levels(y)))
}
```

This is what I get with p_thresh=.01 and forcing mtry to 1. 

```
[1] "LOOCV 56 / 56"
      ROC      Sens      Spec 
0.6238710 0.5200000 0.7419355 
```

If I go up to p_thresh=.05, I get:

```
ROC      Sens      Spec
0.6251613 0.4800000 0.6774194
```

What if I let mtry vary a bit?

```r
nfolds = 5
nreps = 10
clf = 'rf'
p_thresh = .05
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)
mygrid = expand.grid(mtry=c(1:5))

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

    ps = sapply(1:ncol(X), function(v) t.test(X_train[y_train=='Case', v],
                                              X_train[y_train=='Control', v])$p.value)
    good_vars = ps < p_thresh
    X_train2 = X_train[, good_vars]
    X_test2 = t(as.matrix(t(as.matrix(X_test))[, good_vars]))
    mygrid = data.frame(mtry=1)

    set.seed(42)
    fit <- train(X_train2, y_train,
                    trControl = fitControl,
                    method = clf,
                    tuneGrid=mygrid,
                    metric='ROC')

    # updated LOOCV predictions
    y_probs = rbind(y_probs, predict(fit, X_test2, type='prob'))
    y_preds = c(y_preds, levels(y)[predict(fit, X_test2)])

    # just some ongoing summary
    dat = cbind(data.frame(obs = y[1:nrow(y_probs)],
                           pred = factor(y_preds, levels=levels(y))),
                y_probs)
    print(twoClassSummary(dat, lev=levels(y)))
}
```

For p<.01:

```
[1] "LOOCV 56 / 56"
      ROC      Sens      Spec 
0.6367742 0.5200000 0.7741935 
```

And for p < .05:

```
[1] "LOOCV 56 / 56"
      ROC      Sens      Spec 
0.6483871 0.4800000 0.7419355 
```

It's clear no generalizing well. Let me try a method that is less needy of
generalization, and then I can also try plugging in the population PCs and see
if they help.


## simple PLS analysis

I could also go with a simpler PLS analysis:
  https://www.ncbi.nlm.nih.gov/pubmed/20656037 like what I used to do for
  neuroimaging...
  
```r
df = data.frame(X)
df$y = y
plsfit <- plsr(y ~ ., ncomp = 10, data=df, validation = "none")
```

Nope, got protection overflow... but I could go ahead and use the caret methods
from
http://topepo.github.io/caret/miscellaneous-model-functions.html#partial-least-squares-discriminant-analysis.
Then, I could run the same sort of permutation and bootstrap analysis in the
results of those functions? If doing them in a CV fashion doesn't help...

```r
plsFit <- plsda(X, y, ncomp = 10)
plsFitBayes <- plsda(X, y, ncomp = 10, probMethod="Bayes")
splsFit <- splsda(X, y, eta=.1, K=10)
splsFitBayes <- splsda(X, y, ncomp = 10, probMethod="Bayes")
```

I need to get to the SVD, so let's derive everything according to the paper. I
should actually go back to the original data, but for now I'll do it with the
pre-processed data just to see if the scripts work:

```r
Xca = X[y=='Case',]
Xco = X[y=='Control',]
Yca = as.numeric(y[y=='Case'])
Yco = as.numeric(y[y=='Control'])-2
# M is nclass by nvars
M = rbind(colMeans(Xco), colMeans(Xca))
Rmc = M - colMeans(M)
S = svd(Rmc)
```

This is working, now it's a matter of permuting to check which components to
take.

```
> S$d
[1] 43.457164  3.306424
```

```r
eig_vals = c()
nperms = 1000
set.seed(42)
for (p in 1:nperms) {
    perm_idx = sample(1:length(y), length(y), replace=F)
    yperm = y[perm_idx]
    Xca_perm = X[yperm=='Case',]
    Xco_perm = X[yperm=='Control',]
    # M is nclass by nvars
    M_perm = rbind(colMeans(Xco_perm), colMeans(Xca_perm))
    Rmc_perm = M_perm - colMeans(M_perm)
    S_perm = svd(Rmc_perm)
    eig_vals = c(eig_vals, max(S_perm$d))
}
sum(eig_vals >= S$d[1])/nperms
```

So, this dataset as is is somewhat significant at p = .03. I'll run the original
data later, and also Caudate. But to compute stability, we need to do:

```r
library(vegan)
eig_vecs = c()
nperms = 1000
set.seed(42)
for (p in 1:nperms) {
    perm_idx = sample(1:length(y), length(y), replace=T)
    yperm = y[perm_idx]
    Xca_perm = X[yperm=='Case',]
    Xco_perm = X[yperm=='Control',]
    # M is nclass by nvars
    M_perm = rbind(colMeans(Xco_perm), colMeans(Xca_perm))
    Rmc_perm = M_perm - colMeans(M_perm)
    S_perm = svd(Rmc_perm)
    pc = procrustes(S$v, S_perm$v)
    eig_vecs = cbind(eig_vecs, pc$Yrot[, 1])
}
v_err = apply(eig_vecs, 1, sd)
vz = S$v/v_err
```

This might work... we have lots of different genes with good Z:

```
> summary(vz[,1])
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
-14.30484  -3.24049  -0.07007  -0.14198   3.05971  12.88531 
```

Sounds good then. Let's run with this.

# TODO
* try
  
  ?
* crank up MBO