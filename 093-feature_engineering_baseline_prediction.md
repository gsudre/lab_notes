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

set.seed(222)
ensemble_1 <- caretEnsemble(model_list,
                            metric = 'RMSE',
                            trControl = fitControl)
summary(ensemble_1)
plot(ensemble_1)

set.seed(222)
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
                        dummy = RMSE(pred_dummy, y_test))
print(pred_RMSE)

set.seed(123)
xgbTree_model <- train(X_train,
                       y_train,
                       trControl = fitControl,
                       method = 'xgbLinear',
                       metric = 'RMSE',
                       preProcess = c('center', 'scale'),
                       importance = TRUE)
plot(varImp(xgbTree_model))
```