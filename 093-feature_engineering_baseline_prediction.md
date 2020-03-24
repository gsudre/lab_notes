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
model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('cforest', 'svmRadial', 'rf', 'svmLinear',
                                       'blackboost', 'blassoAveraged', 'glmboost',
                                       'kernelpls'),
                        tuneList = NULL,
                        continue_on_fail = FALSE,
                        preProcess = c('center','scale'))

options(digits = 3)
model_results <- data.frame(
 cforest = mean(model_list$cforest$results$RMSE),
 SVMr = mean(model_list$svmRadial$results$RMSE),
 RF = mean(model_list$rf$results$RMSE),
 SVMl = mean(model_list$svmLinear$results$RMSE),
 boostedTree = mean(model_list$blackboost$results$RMSE),
 bLassoAveraged = mean(model_list$blassoAveraged$results$RMSE),
 boostedTrees = mean(model_list$glmboost$results$RMSE),
 PLS = mean(model_list$kernelpls$results$RMSE)
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


sibs by end of todat
shuffle slopes
feature engineering for another day
finish derek by end of weekend
dti hi not as important
maybe show different models, like 4 or so as well
laptops older than 4y