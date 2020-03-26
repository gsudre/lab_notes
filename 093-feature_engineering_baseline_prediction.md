# 2020-03-24 11:01:03

Let's do some training/testing in the anatomy dataset first, and then see if
some feature engineering might help the results. For best in family I'll choose
the youngest kid, just as a heuristic.

```r
fname = '~/data/baseline_prediction/prs_start/gf_impute_based_dti_165.csv'
my_sx = 'hi'

nfolds = 10
nreps = 10

library(caret)
data = read.csv(fname)
# so they don't get rescaled
if (grepl(x=fname, pattern='anat')) {
    mymod = 'anat'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3 = as.factor(data$SES_group3)
    var_names = colnames(data)[c(10:17, 18:29, 30:34, 4:6)]
    phen = sprintf("slope_%s_res_trim.x", my_sx)
} else {
    mymod = 'dti'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3_165 = as.factor(data$SES_group3_165)
    var_names = colnames(data)[c(21:28, 29:40, 41:45, 5:6, 95, 9:20)]
    phen = sprintf("slope_%s_res_trim", my_sx)
}

# if imputation is fair game before we split, so is scaling
scale_me = c()
for (v in var_names) {
    if (!is.factor(data[, v])) {
        scale_me = c(scale_me, v)
    } else {
        data[, v] = as.numeric(data[, v])
    }
}
data[, scale_me] = scale(data[, scale_me])

# split traing and test between members of the same family
train_rows = c()
for (fam in unique(data$FAMID)) {
    fam_rows = which(data$FAMID == fam)
    if (length(fam_rows) == 1) {
        train_rows = c(train_rows, fam_rows[1])
    } else {
        # choose the youngest kid in the family for training
        train_rows = c(train_rows,
                       fam_rows[which.min(data[fam_rows, 'base_age'])])
    }
}
X_train <- data[train_rows, var_names]
X_test <- data[-train_rows, var_names]
y_train <- data[train_rows, phen]
y_test <- data[-train_rows, phen]
```

```r
library(caret)
library(caretEnsemble)
library(doParallel)
registerDoParallel(2)
getDoParWorkers()
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('lm', 'svmRadial', 'rf',
                                       'xgbTree', 'xgbLinear'),
                        tuneList = NULL,
                        continue_on_fail = FALSE,
                        preProcess = c('center','scale'))

options(digits = 3)
model_results <- data.frame(
 LM = mean(model_list$lm$results$RMSE),
 SVM = mean(model_list$svmRadial$results$RMSE),
 RF = mean(model_list$rf$results$RMSE),
 XGBT = mean(model_list$xgbTree$results$RMSE),
 XGBL = mean(model_list$xgbLinear$results$RMSE)
 )
print(model_results)

resamples <- resamples(model_list)
dotplot(resamples, metric = 'RMSE')

modelCor(resamples)

set.seed(42)
ensemble_1 <- caretEnsemble(model_list,
                            metric = 'RMSE',
                            trControl = fitControl)
summary(ensemble_1)
plot(ensemble_1)

set.seed(42)
ensemble_2 <- caretStack(model_list, 
                         method = 'glmnet',
                         metric = 'RMSE',
                         trControl = fitControl)
print(ensemble_2)

# PREDICTIONS
pred_lm <- predict.train(model_list$lm, newdata = X_test)
pred_svm <- predict.train(model_list$svmRadial, newdata = X_test)
pred_rf <- predict.train(model_list$rf, newdata = X_test)
pred_xgbT <- predict.train(model_list$xgbTree, newdata = X_test)
pred_xgbL <- predict.train(model_list$xgbLinear, newdata = X_test)
predict_ens1 <- predict(ensemble_1, newdata = X_test)
predict_ens2 <- predict(ensemble_2, newdata = X_test)
predict_dummy = rep(mean(y_train), length(y_test))
# RMSE
pred_RMSE <- data.frame(ensemble_1 = RMSE(predict_ens1, y_test),
                        ensemble_2 = RMSE(predict_ens2, y_test),
                        LM = RMSE(pred_lm, y_test),
                        SVM = RMSE(pred_svm, y_test),
                        RF = RMSE(pred_rf, y_test),
                        XGBT = RMSE(pred_xgbT, y_test),
                        XGBL = RMSE(pred_xgbL, y_test),
                        dummy = RMSE(predict_dummy, y_test))
print(pred_RMSE)

set.seed(42)
best_model <- train(X_train,
                       y_train,
                       trControl = fitControl,
                       method = 'svmRadial',
                       metric = 'RMSE',
                       preProcess = c('center', 'scale'),
                       importance = TRUE)
plot(varImp(best_model))
```

So, without any sort of feature engineering, we're doing better without
ensembles, both in the training set and the test set:

```
> print(model_results)
     LM   SVM    RF XGBT  XGBL
1 0.581 0.491 0.496 0.57 0.577
> print(ensemble_1)
A glm ensemble of 5 base models: lm, svmRadial, rf, xgbTree, xgbLinear

Ensemble results:
Generalized Linear Model 

1250 samples
   5 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1125, 1125, 1125, 1125, 1125, 1125, ... 
Resampling results:

  RMSE   Rsquared  MAE  
  0.495  0.0317    0.345
> print(ensemble_2)
A glmnet ensemble of 5 base models: lm, svmRadial, rf, xgbTree, xgbLinear

Ensemble results:
glmnet 

1250 samples
   5 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1125, 1125, 1125, 1125, 1125, 1125, ... 
Resampling results across tuning parameters:

  alpha  lambda    RMSE   Rsquared  MAE  
  0.10   0.000165  0.495  0.0317    0.345
  0.10   0.001652  0.495  0.0318    0.345
  0.10   0.016518  0.495  0.0322    0.345
  0.55   0.000165  0.495  0.0317    0.345
  0.55   0.001652  0.495  0.0320    0.345
  0.55   0.016518  0.495  0.0322    0.344
  1.00   0.000165  0.495  0.0318    0.345
  1.00   0.001652  0.495  0.0322    0.345
  1.00   0.016518  0.495  0.0329    0.344

RMSE was used to select the optimal model using the smallest value.
The final values used for the model were alpha = 0.1 and lambda = 0.0165.
```

Combining them using weighted averages or glmnet doesn't make a difference in
terms of RMSE in the training set. In fact, SVM does best, then the ensembles,
followed closely by rf.

```
> print(pred_RMSE)
  ensemble_1 ensemble_2    LM   SVM    RF  XGBT  XGBL dummy
1      0.401      0.399 0.608 0.359 0.389 0.386 0.397 0.377
```

The test results loosely followed the training set, but SVM was the only one
that did better than dummy. The consolation is that dummy will be the same
forever provided that training and test sets stay the same. But we can make the
other results better through feature engineering and/or better models. Or even
better tunning the specific models.

So, let's start there, mostly using the models we've had success with:

```r
library(caret)
library(caretEnsemble)
library(doParallel)
registerDoParallel(31)
getDoParWorkers()
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('cforest', 'svmRadial', 'rf', 'svmLinear',
                                       'blackboost', 'blassoAveraged', 'glmboost',
                                       'kernelpls'),
                        tuneList = NULL,
                        continue_on_fail = FALSE,
                        preProcess = c('center','scale'))
null_fit <- train(x=X_train, y=y_train, method = 'null', trControl = fitControl)

options(digits = 3)
model_results <- data.frame(
 cforest = mean(model_list$cforest$results$RMSE),
 SVMr = mean(model_list$svmRadial$results$RMSE),
 RF = mean(model_list$rf$results$RMSE),
 SVMl = mean(model_list$svmLinear$results$RMSE),
 boostedTree = mean(model_list$blackboost$results$RMSE),
 bLassoAveraged = mean(model_list$blassoAveraged$results$RMSE),
 glmBoost = mean(model_list$glmboost$results$RMSE),
 PLS = mean(model_list$kernelpls$results$RMSE),
 null_model = null_fit$results$RMSE
 )
print(model_results)

resamples <- resamples(model_list)
dotplot(resamples, metric = 'RMSE')

modelCor(resamples)

set.seed(42)
ensemble_1 <- caretEnsemble(model_list,
                            metric = 'RMSE',
                            trControl = fitControl)
summary(ensemble_1)
plot(ensemble_1)

set.seed(42)
ensemble_2 <- caretStack(model_list, 
                         method = 'glmnet',
                         metric = 'RMSE',
                         trControl = fitControl)
print(ensemble_2)

# PREDICTIONS
pred_cforest <- predict.train(model_list$cforest, newdata = X_test)
pred_SVMr <- predict.train(model_list$svmRadial, newdata = X_test)
pred_RF <- predict.train(model_list$rf, newdata = X_test)
pred_svmL <- predict.train(model_list$svmLinear, newdata = X_test)
pred_boostedTree <- predict.train(model_list$blackboost, newdata = X_test)
pred_bLassoAveraged <- predict(model_list$blassoAveraged, newdata = X_test)
pred_glmBoost <- predict(model_list$glmboost, newdata = X_test)
pred_pls <- predict(model_list$kernelpls, newdata = X_test)
pred_ens1 <- predict(ensemble_1, newdata = X_test)
pred_ens2 <- predict(ensemble_2, newdata = X_test)
pred_null = predict(null_fit, newdata = X_test)
# RMSE
pred_RMSE <- data.frame(ensemble_1 = RMSE(pred_ens1, y_test),
                        ensemble_2 = RMSE(pred_ens2, y_test),
                        cforest = RMSE(pred_cforest, y_test),
                        SVMr = RMSE(pred_SVMr, y_test),
                        RF = RMSE(pred_RF, y_test),
                        SVMl = RMSE(pred_svmL, y_test),
                        boostedTrees = RMSE(pred_boostedTree, y_test),
                        blassoAveraged = RMSE(pred_bLassoAveraged, y_test),
                        glmBoost = RMSE(pred_glmBoost, y_test),
                        kernelpls = RMSE(pred_pls, y_test),
                        null = RMSE(pred_null, y_test))
print(pred_RMSE)
```

If I use our better classifiers from the past, it's not doing much better:

```
> print(model_results)
  cforest  SVMr    RF  SVMl boostedTree bLassoAveraged glmBoost   PLS null_model
1   0.492 0.491 0.496 0.552       0.485          0.484    0.491 0.503      0.486
```

Not much doing better than the null model...

```
>  print(ensemble_1)
A glm ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
Generalized Linear Model 

1250 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1125, 1125, 1125, 1125, 1125, 1125, ... 
Resampling results:

  RMSE  Rsquared  MAE  
  0.48  0.0896    0.341

>  print(ensemble_2)
A glmnet ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
glmnet 

1250 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1125, 1125, 1125, 1125, 1125, 1125, ... 
Resampling results across tuning parameters:

  alpha  lambda    RMSE   Rsquared  MAE  
  0.10   0.000133  0.480  0.0895    0.340
  0.10   0.001331  0.480  0.0890    0.340
  0.10   0.013315  0.483  0.0799    0.340
  0.55   0.000133  0.480  0.0895    0.340
  0.55   0.001331  0.480  0.0888    0.340
  0.55   0.013315  0.487  0.0686    0.341
  1.00   0.000133  0.480  0.0895    0.340
  1.00   0.001331  0.480  0.0886    0.340
  1.00   0.013315  0.494  0.0389    0.346

RMSE was used to select the optimal model using the smallest value.
The final values used for the model were alpha = 1 and lambda = 0.000133.
```

Although our ensembles do even a bit better than the best predictor (svmRadial),
it fails during testing:

```
> print(pred_RMSE)
  ensemble_1 ensemble_2 cforest  SVMr    RF  SVMl boostedTrees blassoAveraged glmBoost kernelpls  null
1      0.388      0.388   0.381 0.359 0.388 0.491        0.377          0.379    0.396     0.422 0.377
```

I wonder if that's a trend I'll see for the other datasets and sx as well...
let's check.

The result above was for DTI HI. Now, for DTI inatt:

```
> print(model_results)
  cforest  SVMr    RF  SVMl boostedTree bLassoAveraged glmBoost   PLS null_model
1   0.585 0.584 0.586 0.744       0.569           0.57    0.596 0.622      0.572
> print(ensemble_1)
A glm ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
Generalized Linear Model 

1250 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1125, 1125, 1125, 1125, 1125, 1125, ... 
Resampling results:

  RMSE   Rsquared  MAE  
  0.556  0.103     0.424

> print(ensemble_2)
A glmnet ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
glmnet 

1250 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1125, 1125, 1125, 1125, 1125, 1125, ... 
Resampling results across tuning parameters:

  alpha  lambda    RMSE   Rsquared  MAE  
  0.10   0.000168  0.556  0.1024    0.424
  0.10   0.001681  0.556  0.1012    0.423
  0.10   0.016810  0.563  0.0820    0.427
  0.55   0.000168  0.556  0.1024    0.424
  0.55   0.001681  0.556  0.1014    0.423
  0.55   0.016810  0.567  0.0704    0.429
  1.00   0.000168  0.556  0.1025    0.424
  1.00   0.001681  0.556  0.1015    0.423
  1.00   0.016810  0.573  0.0461    0.434

RMSE was used to select the optimal model using the smallest value.
The final values used for the model were alpha = 1 and lambda = 0.000168.

> print(pred_RMSE)
  ensemble_1 ensemble_2 cforest  SVMr    RF  SVMl boostedTrees blassoAveraged glmBoost kernelpls  null
1       0.55      0.549   0.546 0.541 0.571 0.678        0.542          0.544     0.55     0.566 0.542

```

Both ensembles do better, but neither does well enough to predict the test set.

How about the anatomical dataset?

```
> print(model_results)
  cforest  SVMr    RF  SVMl boostedTree bLassoAveraged glmBoost   PLS null_model
1   0.626 0.626 0.631 0.669       0.629          0.624     0.63 0.641      0.624
> print(ensemble_1)
A glm ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
Generalized Linear Model 

1820 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1638, 1638, 1638, 1638, 1638, 1638, ... 
Resampling results:

  RMSE   Rsquared  MAE  
  0.612  0.0889    0.452

> print(ensemble_2)
A glmnet ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
glmnet 

1820 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1638, 1638, 1638, 1638, 1638, 1638, ... 
Resampling results across tuning parameters:

  alpha  lambda    RMSE   Rsquared  MAE  
  0.10   0.000191  0.612  0.0888    0.452
  0.10   0.001914  0.612  0.0878    0.451
  0.10   0.019138  0.619  0.0683    0.453
  0.55   0.000191  0.612  0.0889    0.452
  0.55   0.001914  0.612  0.0883    0.451
  0.55   0.019138  0.626  0.0465    0.456
  1.00   0.000191  0.612  0.0889    0.452
  1.00   0.001914  0.612  0.0886    0.451
  1.00   0.019138  0.628  0.0411    0.457

RMSE was used to select the optimal model using the smallest value.
The final values used for the model were alpha = 1 and lambda = 0.000191.
> print(pred_RMSE)
  ensemble_1 ensemble_2 cforest  SVMr    RF SVMl boostedTrees blassoAveraged glmBoost kernelpls  null
1      0.645      0.645   0.642 0.632 0.633 0.67        0.642          0.639    0.648     0.644 0.642
```

Again, the ensembles do better at training. Just not good enough as the simpler
models for testing.

Finally, we didn't test anatomical dataset HI:

```
> print(model_results)
  cforest SVMr    RF  SVMl boostedTree bLassoAveraged glmBoost   PLS null_model
1   0.563 0.57 0.573 0.572       0.567          0.567    0.562 0.583      0.572
> print(ensemble_1)
A glm ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
Generalized Linear Model 

1820 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1638, 1638, 1638, 1638, 1638, 1638, ... 
Resampling results:

  RMSE   Rsquared  MAE  
  0.552  0.128     0.385

> print(ensemble_2)
A glmnet ensemble of 8 base models: cforest, svmRadial, rf, svmLinear, blackboost, blassoAveraged, glmboost, kernelpls

Ensemble results:
glmnet 

1820 samples
   8 predictor

No pre-processing
Resampling: Cross-Validated (10 fold, repeated 10 times) 
Summary of sample sizes: 1638, 1638, 1638, 1638, 1638, 1638, ... 
Resampling results across tuning parameters:

  alpha  lambda    RMSE   Rsquared  MAE  
  0.10   0.000238  0.552  0.128     0.385
  0.10   0.002384  0.552  0.128     0.384
  0.10   0.023839  0.553  0.127     0.378
  0.55   0.000238  0.552  0.128     0.385
  0.55   0.002384  0.552  0.128     0.384
  0.55   0.023839  0.554  0.126     0.374
  1.00   0.000238  0.552  0.128     0.385
  1.00   0.002384  0.552  0.128     0.383
  1.00   0.023839  0.557  0.122     0.370

RMSE was used to select the optimal model using the smallest value.
The final values used for the model were alpha = 0.55 and lambda = 0.00238.
> print(pred_RMSE)
  ensemble_1 ensemble_2 cforest  SVMr    RF  SVMl boostedTrees blassoAveraged glmBoost kernelpls  null
1      0.454      0.452   0.458 0.439 0.438 0.443        0.436          0.434    0.443     0.459 0.436
```

Same story...

Would this get better if I went for LOOCV, since I don't have much data to begin
with? How about only 3-fold?

## 3-fold CV

caretList is crashing with LOOCV, so I'll check out just 3-fold for now.

... no difference. Same pattern of results :()

In terms of reporting I think we could show the model selection (for best model)
using the confidence intervals from resamples, and show the null model
threshold. Then, we show how well we do in predicting a set of never seen data.

Someting like this, for example:

https://topepo.github.io/caret/model-training-and-tuning.html#plots

But how can we engineer these features for optimal performance? And even the
target?

# 2020-03-25 19:46:46

Let's start with just the usual YeoJohnson transform, without paying too close
attention to anything else... that didn't work.

Then, I tried to TRULYnormalize our target:

```r
fname = '~/data/baseline_prediction/prs_start/gf_impute_based_anatomy_272.csv'
my_sx = 'inatt'

nfolds = 10
nreps = 10

library(caret)
data = read.csv(fname)
# so they don't get rescaled
if (grepl(x=fname, pattern='anat')) {
    mymod = 'anat'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3 = as.factor(data$SES_group3)
    var_names = colnames(data)[c(10:17, 18:29, 30:34, 4:6)]
    phen = sprintf("slope_%s_res_trim.x", my_sx)
} else {
    mymod = 'dti'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3_165 = as.factor(data$SES_group3_165)
    var_names = colnames(data)[c(21:28, 29:40, 41:45, 5:6, 95, 9:20)]
    phen = sprintf("slope_%s_res_trim", my_sx)
}

library(caret)
pp = preProcess(data[, var_names], method = c("YeoJohnson", 'range'))
data2 = predict(pp, data[, var_names])
library(bestNormalize)
bn = bestNormalize(data[, phen])
data2$phen = bn$x.t
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# split traing and test between members of the same family
train_rows = c()
for (fam in unique(data$FAMID)) {
    fam_rows = which(data$FAMID == fam)
    if (length(fam_rows) == 1) {
        train_rows = c(train_rows, fam_rows[1])
    } else {
        # choose the youngest kid in the family for training
        train_rows = c(train_rows,
                       fam_rows[which.min(data[fam_rows, 'base_age'])])
    }
}
# data3 doesn't have the target column!
X_train <- data3[train_rows, ]
X_test <- data3[-train_rows, ]
y_train <- data2[train_rows,]$phen
y_test <- data2[-train_rows,]$phen
```

```r
library(caretEnsemble)
library(doParallel)
registerDoParallel(31)
getDoParWorkers()
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('cforest', 'svmRadial', 'rf', 'svmLinear',
                                       'blackboost', 'blassoAveraged', 'glmboost',
                                       'kernelpls'),
                        tuneList = NULL,
                        continue_on_fail = FALSE)
null_fit <- train(x=X_train, y=y_train, method = 'null', trControl = fitControl)

set.seed(42)
ensemble_1 <- caretEnsemble(model_list,
                            metric = 'RMSE',
                            trControl = fitControl)

set.seed(42)
ensemble_2 <- caretStack(model_list, 
                         method = 'glmnet',
                         metric = 'RMSE',
                         trControl = fitControl)

options(digits = 3)
model_results <- data.frame(
 cforest = min(model_list$cforest$results$RMSE),
 SVMr = min(model_list$svmRadial$results$RMSE),
 RF = min(model_list$rf$results$RMSE),
 SVMl = min(model_list$svmLinear$results$RMSE),
 boostedTree = min(model_list$blackboost$results$RMSE),
 bLassoAveraged = min(model_list$blassoAveraged$results$RMSE),
 glmBoost = min(model_list$glmboost$results$RMSE),
 PLS = min(model_list$kernelpls$results$RMSE),
 ens_1 = ensemble_1$ens_model$results$RMSE,
 ens2 = min(ensemble_2$ens_model$results$RMSE),
 null_model = null_fit$results$RMSE
 )
print(model_results)

resamples <- resamples(model_list)
modelCor(resamples)

# PREDICTIONS
pred_cforest <- predict.train(model_list$cforest, newdata = X_test)
pred_SVMr <- predict.train(model_list$svmRadial, newdata = X_test)
pred_RF <- predict.train(model_list$rf, newdata = X_test)
pred_svmL <- predict.train(model_list$svmLinear, newdata = X_test)
pred_boostedTree <- predict.train(model_list$blackboost, newdata = X_test)
pred_bLassoAveraged <- predict(model_list$blassoAveraged, newdata = X_test)
pred_glmBoost <- predict(model_list$glmboost, newdata = X_test)
pred_pls <- predict(model_list$kernelpls, newdata = X_test)
pred_ens1 <- predict(ensemble_1, newdata = X_test)
pred_ens2 <- predict(ensemble_2, newdata = X_test)
pred_null = predict(null_fit, newdata = X_test)
# RMSE
pred_RMSE <- data.frame(ensemble_1 = RMSE(pred_ens1, y_test),
                        ensemble_2 = RMSE(pred_ens2, y_test),
                        cforest = RMSE(pred_cforest, y_test),
                        SVMr = RMSE(pred_SVMr, y_test),
                        RF = RMSE(pred_RF, y_test),
                        SVMl = RMSE(pred_svmL, y_test),
                        boostedTrees = RMSE(pred_boostedTree, y_test),
                        blassoAveraged = RMSE(pred_bLassoAveraged, y_test),
                        glmBoost = RMSE(pred_glmBoost, y_test),
                        kernelpls = RMSE(pred_pls, y_test),
                        null = RMSE(pred_null, y_test))
print(pred_RMSE)
```

This actually helped a bit! Hard to tell the meaning of the error though... but
at least we're learning better. Let me remove some of the classifiers that are
not contributing much to the model and try it again. Also, I'll scale the output
to be between 0 and 1 just because... and for that matter, I'll normalize every
single numeric predictor the best way we can:

```r
fname = '~/data/baseline_prediction/prs_start/gf_impute_based_anatomy_272.csv'
my_sx = 'hi'

nfolds = 10
nreps = 10

library(caret)
data = read.csv(fname)
# so they don't get rescaled
if (grepl(x=fname, pattern='anat')) {
    mymod = 'anat'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3 = as.factor(data$SES_group3)
    var_names = colnames(data)[c(10:17, 18:29, 30:34, 4:6)]
    phen = sprintf("slope_%s_res_trim.x", my_sx)
} else {
    mymod = 'dti'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3_165 = as.factor(data$SES_group3_165)
    var_names = colnames(data)[c(21:28, 29:40, 41:45, 5:6, 95, 9:20)]
    phen = sprintf("slope_%s_res_trim", my_sx)
}

library(caret)
library(bestNormalize)
# mostly so Box-Cox is an option
pp = preProcess(data[, var_names], method = c('range'), rangeBounds=c(1, 2))
data2 = predict(pp, data[, var_names])
scale_me = c()
for (v in var_names) {
    if (!is.factor(data2[, v])) {
        print(sprintf('Normalizing %s', v))
        bn = bestNormalize(data2[, v])
        data2[, v] = bn$x.t
    }
}
# not so worried about boxcox here, it's almost never needed
bn = bestNormalize(data[, phen])
data2$phen = bn$x.t
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# split traing and test between members of the same family
train_rows = c()
for (fam in unique(data$FAMID)) {
    fam_rows = which(data$FAMID == fam)
    if (length(fam_rows) == 1) {
        train_rows = c(train_rows, fam_rows[1])
    } else {
        # choose the youngest kid in the family for training
        train_rows = c(train_rows,
                       fam_rows[which.min(data[fam_rows, 'base_age'])])
    }
}
# data3 doesn't have the target column!
X_train <- data3[train_rows, ]
X_test <- data3[-train_rows, ]
y_train <- data2[train_rows,]$phen
y_test <- data2[-train_rows,]$phen
```

```r
library(caretEnsemble)
library(doParallel)
registerDoParallel(31)
getDoParWorkers()
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('cforest', 'svmRadial', 'rf', #'svmLinear',
                                       'blackboost', 'blassoAveraged', 'glmboost',
                                       'kernelpls'),
                        tuneList = NULL,
                        continue_on_fail = FALSE)
null_fit <- train(x=X_train, y=y_train, method = 'null', trControl = fitControl)

set.seed(42)
ensemble_1 <- caretEnsemble(model_list,
                            metric = 'RMSE',
                            trControl = fitControl)

set.seed(42)
ensemble_2 <- caretStack(model_list, 
                         method = 'glmnet',
                         metric = 'RMSE',
                         trControl = fitControl)

options(digits = 3)
model_results <- data.frame(
 cforest = min(model_list$cforest$results$RMSE),
 SVMr = min(model_list$svmRadial$results$RMSE),
 RF = min(model_list$rf$results$RMSE),
#  SVMl = min(model_list$svmLinear$results$RMSE),
 boostedTree = min(model_list$blackboost$results$RMSE),
 bLassoAveraged = min(model_list$blassoAveraged$results$RMSE),
 glmBoost = min(model_list$glmboost$results$RMSE),
 PLS = min(model_list$kernelpls$results$RMSE),
 ens_1 = ensemble_1$ens_model$results$RMSE,
 ens2 = min(ensemble_2$ens_model$results$RMSE),
 null_model = null_fit$results$RMSE
 )

resamples <- resamples(model_list)
modelCor(resamples)

# PREDICTIONS
pred_cforest <- predict.train(model_list$cforest, newdata = X_test)
pred_SVMr <- predict.train(model_list$svmRadial, newdata = X_test)
pred_RF <- predict.train(model_list$rf, newdata = X_test)
# pred_svmL <- predict.train(model_list$svmLinear, newdata = X_test)
pred_boostedTree <- predict.train(model_list$blackboost, newdata = X_test)
pred_bLassoAveraged <- predict(model_list$blassoAveraged, newdata = X_test)
pred_glmBoost <- predict(model_list$glmboost, newdata = X_test)
pred_pls <- predict(model_list$kernelpls, newdata = X_test)
pred_ens1 <- predict(ensemble_1, newdata = X_test)
pred_ens2 <- predict(ensemble_2, newdata = X_test)
pred_null = predict(null_fit, newdata = X_test)
# RMSE
pred_RMSE <- data.frame(cforest = RMSE(pred_cforest, y_test),
                        SVMr = RMSE(pred_SVMr, y_test),
                        RF = RMSE(pred_RF, y_test),
                        # SVMl = RMSE(pred_svmL, y_test),
                        boostedTrees = RMSE(pred_boostedTree, y_test),
                        blassoAveraged = RMSE(pred_bLassoAveraged, y_test),
                        glmBoost = RMSE(pred_glmBoost, y_test),
                        kernelpls = RMSE(pred_pls, y_test),
                        ensemble_1 = RMSE(pred_ens1, y_test),
                        ensemble_2 = RMSE(pred_ens2, y_test),
                        null = RMSE(pred_null, y_test))
print(model_results)
print(pred_RMSE)
```

After some engineering, we might have some better results. For example:

```
anatomy; inatt:
> print(model_results)
  cforest SVMr RF boostedTree bLassoAveraged glmBoost   PLS ens_1  ens2 null_model
1   0.994 1.01  1           1          0.996        1 0.991 0.946 0.946      0.998
> print(pred_RMSE)
  cforest  SVMr    RF boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.968 0.974 0.962        0.973          0.967    0.977     0.977      0.989      0.988 0.973
```

Well, the good results ended right after that :(

```
anatomy; hi:
> print(model_results)
  cforest SVMr   RF boostedTree bLassoAveraged glmBoost  PLS ens_1  ens2 null_model
1    1.03 1.04 1.03        1.04           1.04     1.03 1.04 0.998 0.998       1.04
> print(pred_RMSE)
  cforest  SVMr    RF boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.883 0.904 0.888        0.882          0.879    0.878     0.899      0.914      0.913 0.882

dti; inatt:
> print(model_results)
  cforest  SVMr    RF boostedTree bLassoAveraged glmBoost  PLS ens_1  ens2 null_model
1    0.99 0.999 0.996       0.988          0.988     1.01 1.01 0.982 0.982       0.99
> print(pred_RMSE)
  cforest SVMr   RF boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.995 1.01 1.02        0.995          0.995    0.993      1.02       1.02       1.02 0.995

dti; hi:
> print(model_results)
  cforest SVMr   RF boostedTree bLassoAveraged glmBoost  PLS ens_1  ens2 null_model
1    1.02 1.02 1.02        1.02           1.02     1.04 1.03 0.994 0.994       1.03
> print(pred_RMSE)
  cforest  SVMr    RF boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.895 0.938 0.901        0.891          0.896    0.897     0.965      0.953      0.951 0.891
```

I could try PCA or ICA within domain? Worth trying... maybe with the normalized
and the regular slope? No need to normalize the predictors before that... just
the output form PCA or ICA should be enough.

# 2020-03-26 07:26:04

```r
fname = '~/data/baseline_prediction/prs_start/gf_impute_based_anatomy_272.csv'
my_sx = 'hi'

nfolds = 10
nreps = 10

library(caret)
data = read.csv(fname)
# so they don't get rescaled
if (grepl(x=fname, pattern='anat')) {
    mymod = 'anat'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3 = as.factor(data$SES_group3)
    var_names = colnames(data)[c(10:17, 18:29, 30:34, 4:6)]
    phen = sprintf("slope_%s_res_trim.x", my_sx)
} else {
    mymod = 'dti'
    data$sex_numeric = as.factor(data$sex_numeric)
    data$SES_group3_165 = as.factor(data$SES_group3_165)
    var_names = colnames(data)[c(21:28, 29:40, 41:45, 5:6, 95, 9:20)]
    phen = sprintf("slope_%s_res_trim", my_sx)
}

library(caret)
# anat only!
pp = preProcess(data[, var_names[1:8]], method = c('center', 'scale', 'pca'))
tmp_data = predict(pp, data[, var_names[1:8]])
cnames = sapply(colnames(tmp_data), function(x) sprintf('struct_%s', x))
colnames(tmp_data) = cnames
data2 = tmp_data
pp = preProcess(data[, var_names[9:20]], method = c('center', 'scale', 'pca'))
tmp_data = predict(pp, data[, var_names[9:20]])
cnames = sapply(colnames(tmp_data), function(x) sprintf('PRS_%s', x))
colnames(tmp_data) = cnames
data2 = cbind(data2, tmp_data)
pp = preProcess(data[, var_names[21:25]], method = c('center', 'scale', 'pca'))
tmp_data = predict(pp, data[, var_names[21:25]])
cnames = sapply(colnames(tmp_data), function(x) sprintf('cog_%s', x))
colnames(tmp_data) = cnames
data2 = cbind(data2, tmp_data)
data2 = cbind(data2, data[, var_names[26:28]])

library(bestNormalize)
bn = bestNormalize(data[, phen])
data2$phen = bn$x.t
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# split traing and test between members of the same family
train_rows = c()
for (fam in unique(data$FAMID)) {
    fam_rows = which(data$FAMID == fam)
    if (length(fam_rows) == 1) {
        train_rows = c(train_rows, fam_rows[1])
    } else {
        # choose the youngest kid in the family for training
        train_rows = c(train_rows,
                       fam_rows[which.min(data[fam_rows, 'base_age'])])
    }
}
# data3 doesn't have the target column!
X_train <- data3[train_rows, ]
X_test <- data3[-train_rows, ]
y_train <- data2[train_rows,]$phen
y_test <- data2[-train_rows,]$phen
```

Now, let's try everything:

```r
library(caretEnsemble)
library(doParallel)
registerDoParallel(31)
getDoParWorkers()
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('cforest', 'svmRadial', 'rf', 'svmLinear',
                                       'blackboost', 'blassoAveraged', 'glmboost',
                                       'kernelpls'),
                        tuneList = NULL,
                        continue_on_fail = FALSE)
null_fit <- train(x=X_train, y=y_train, method = 'null', trControl = fitControl)

set.seed(42)
ensemble_1 <- caretEnsemble(model_list,
                            metric = 'RMSE',
                            trControl = fitControl)

set.seed(42)
ensemble_2 <- caretStack(model_list, 
                         method = 'glmnet',
                         metric = 'RMSE',
                         trControl = fitControl)

options(digits = 3)
model_results <- data.frame(
 cforest = min(model_list$cforest$results$RMSE),
 SVMr = min(model_list$svmRadial$results$RMSE),
 RF = min(model_list$rf$results$RMSE),
 SVMl = min(model_list$svmLinear$results$RMSE),
 boostedTree = min(model_list$blackboost$results$RMSE),
 bLassoAveraged = min(model_list$blassoAveraged$results$RMSE),
 glmBoost = min(model_list$glmboost$results$RMSE),
 PLS = min(model_list$kernelpls$results$RMSE),
 ens_1 = ensemble_1$ens_model$results$RMSE,
 ens2 = min(ensemble_2$ens_model$results$RMSE),
 null_model = null_fit$results$RMSE
 )

resamples <- resamples(model_list)
modelCor(resamples)

# PREDICTIONS
pred_cforest <- predict.train(model_list$cforest, newdata = X_test)
pred_SVMr <- predict.train(model_list$svmRadial, newdata = X_test)
pred_RF <- predict.train(model_list$rf, newdata = X_test)
pred_svmL <- predict.train(model_list$svmLinear, newdata = X_test)
pred_boostedTree <- predict.train(model_list$blackboost, newdata = X_test)
pred_bLassoAveraged <- predict(model_list$blassoAveraged, newdata = X_test)
pred_glmBoost <- predict(model_list$glmboost, newdata = X_test)
pred_pls <- predict(model_list$kernelpls, newdata = X_test)
pred_ens1 <- predict(ensemble_1, newdata = X_test)
pred_ens2 <- predict(ensemble_2, newdata = X_test)
pred_null = predict(null_fit, newdata = X_test)
# RMSE
pred_RMSE <- data.frame(cforest = RMSE(pred_cforest, y_test),
                        SVMr = RMSE(pred_SVMr, y_test),
                        RF = RMSE(pred_RF, y_test),
                        SVMl = RMSE(pred_svmL, y_test),
                        boostedTrees = RMSE(pred_boostedTree, y_test),
                        blassoAveraged = RMSE(pred_bLassoAveraged, y_test),
                        glmBoost = RMSE(pred_glmBoost, y_test),
                        kernelpls = RMSE(pred_pls, y_test),
                        ensemble_1 = RMSE(pred_ens1, y_test),
                        ensemble_2 = RMSE(pred_ens2, y_test),
                        null = RMSE(pred_null, y_test))
print(model_results)
print(pred_RMSE)
```

```
anat; inatt:
> print(model_results)
  cforest SVMr   RF SVMl boostedTree bLassoAveraged glmBoost   PLS ens_1  ens2 null_model
1       1 1.01 1.01 1.04           1          0.999     0.98 0.991 0.919 0.919      0.998
> print(pred_RMSE)
  cforest  SVMr   RF SVMl boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.971 0.984 0.98 1.05        0.973          0.969    0.974     0.978       1.02       1.02 0.973

anat; hi:
> print(model_results)
  cforest SVMr   RF SVMl boostedTree bLassoAveraged glmBoost  PLS ens_1 ens2 null_model
1    1.04 1.04 1.05 1.13        1.04           1.04     1.03 1.04  1.01 1.01       1.04
> print(pred_RMSE)
  cforest SVMr    RF  SVMl boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.878 0.89 0.901 0.946        0.882          0.879     0.89     0.897      0.932      0.925 0.882
```

I can also try to reduce my PCA variance to decrease the number of variables...
now, using 80%:

```r
# anat only!
pp = preProcess(data[, var_names[1:8]], method = c('center', 'scale', 'ica'), n.comp=4)
tmp_data = predict(pp, data[, var_names[1:8]])
cnames = sapply(colnames(tmp_data), function(x) sprintf('struct_%s', x))
colnames(tmp_data) = cnames
data2 = tmp_data
pp = preProcess(data[, var_names[9:20]], method = c('center', 'scale', 'ica'), n.comp=5)
tmp_data = predict(pp, data[, var_names[9:20]])
cnames = sapply(colnames(tmp_data), function(x) sprintf('PRS_%s', x))
colnames(tmp_data) = cnames
data2 = cbind(data2, tmp_data)
pp = preProcess(data[, var_names[21:25]], method = c('center', 'scale', 'ica'), n.comp=2)
tmp_data = predict(pp, data[, var_names[21:25]])
cnames = sapply(colnames(tmp_data), function(x) sprintf('cog_%s', x))
colnames(tmp_data) = cnames
data2 = cbind(data2, tmp_data)
data2 = cbind(data2, data[, var_names[26:28]])

library(bestNormalize)
bn = bestNormalize(data[, phen])
data2$phen = bn$x.t
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# split traing and test between members of the same family
train_rows = c()
for (fam in unique(data$FAMID)) {
    fam_rows = which(data$FAMID == fam)
    if (length(fam_rows) == 1) {
        train_rows = c(train_rows, fam_rows[1])
    } else {
        # choose the youngest kid in the family for training
        train_rows = c(train_rows,
                       fam_rows[which.min(data[fam_rows, 'base_age'])])
    }
}
# data3 doesn't have the target column!
X_train <- data3[train_rows, ]
X_test <- data3[-train_rows, ]
y_train <- data2[train_rows,]$phen
y_test <- data2[-train_rows,]$phen
```

```
anat; inatt:
> print(model_results)
  cforest SVMr   RF SVMl boostedTree bLassoAveraged glmBoost   PLS ens_1  ens2 null_model
1       1 1.02 1.01 1.07           1              1        1 0.991 0.961 0.961      0.999
> print(pred_RMSE)
  cforest  SVMr    RF  SVMl boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.964 0.965 0.948 0.972        0.973          0.969    0.959     0.975      0.987      0.987 0.973

anat; hi:
> print(model_results)
  cforest SVMr   RF SVMl boostedTree bLassoAveraged glmBoost  PLS ens_1 ens2 null_model
1    1.04 1.04 1.05 1.09        1.04           1.04     1.02 1.04  0.98 0.98       1.04
> print(pred_RMSE)
  cforest  SVMr    RF SVMl boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.876 0.917 0.912 0.97        0.882          0.879    0.891     0.895       1.02       1.02 0.882
```

How about ICA?

```
anat; inatt:
> print(model_results)
  cforest SVMr   RF SVMl boostedTree bLassoAveraged glmBoost  PLS ens_1  ens2 null_model
1       1 1.01 1.01 1.05           1              1     1.01 1.01 0.936 0.936      0.999
> print(pred_RMSE)
  cforest  SVMr    RF SVMl boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.966 0.967 0.962 1.01        0.973          0.969    0.958     0.975       1.01       1.01 0.973

anat; hi:
> print(model_results)
  cforest SVMr   RF SVMl boostedTree bLassoAveraged glmBoost  PLS ens_1 ens2 null_model
1    1.04 1.04 1.05  1.1        1.04           1.04     1.05 1.05  1.02 1.02       1.04
> print(pred_RMSE)
  cforest  SVMr  RF  SVMl boostedTrees blassoAveraged glmBoost kernelpls ensemble_1 ensemble_2  null
1   0.881 0.909 0.9 0.979        0.882          0.878    0.877     0.899      0.903      0.901 0.882
```

Nothing...