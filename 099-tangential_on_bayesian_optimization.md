# 2020-04-08 11:48:27

Let's see if we can narrow down a parameter using Bayesian optimization.

```r
library(caret)
library(xgboost)
library(rBayesianOptimization)

fname = '~/data/baseline_prediction/prs_start/gf_philip_03292020.csv'
phen = 'categ_all.4'
c1 = 'improvers'
c2 = 'stable_symptomatic'

data = read.csv(fname)
data$sex_numeric = as.factor(data$sex_numeric)
data$SES_group3 = as.factor(data$SES_group3)
data$slf_fa = data$slf_all  # just to make it easier to filter out
var_names = c(# PRS
              'ADHD_PRS0.000100.orig', 'ADHD_PRS0.001000.orig',
              'ADHD_PRS0.010000.orig', 'ADHD_PRS0.050000.orig',
              'ADHD_PRS0.100000.orig', 'ADHD_PRS0.200000.orig',
              'ADHD_PRS0.300000.orig', 'ADHD_PRS0.400000.orig',
              'ADHD_PRS0.500000.orig',
              # DTI
              'atr_fa', 'cst_fa', 'cing_cing_fa', 'cing_hipp_fa', 'cc_fa',
              'ilf_fa', 'slf_fa', 'unc_fa',
              # demo
              'sex_numeric', 'SES_group3',
              # cog
              'FSIQ', 'SS_RAW', 'DS_RAW', 'PS_RAW', 'VMI.beery_RAW',
              # anat
              'cerebellum_white', 'cerebellum_grey', 'amygdala',
              'cingulate', 'lateral_PFC', 'OFC', 'striatum', 'thalamus'
              )

# data2 = data[, c(var_names, 'FAMID')]
data2 = data[, var_names]

data2$phen = as.factor(data[, phen])
# data2$FAMID = as.factor(data$FAMID)
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# selecting only kids in the 2 specified groups
keep_me = data2$phen==c1 | data2$phen==c2
X = data3[keep_me, ]
y = factor(data2[keep_me, 'phen'])

test_row = 1
train_rows = setdiff(1:nrow(X), test_row)
X_train <- X[train_rows, ]
X_test <- X[-train_rows, ]
y_train <- y[train_rows]
y_test <- y[-train_rows]

dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = as.numeric(y_train)-1, missing=NA)
dtest <- xgb.DMatrix(data = t(as.matrix(X_test)), label = as.numeric(y_test)-1, missing=NA)
```

```r
set.seed(42)
cv_folds = createFolds(y, k = 10)
xgb_cv_bayes <- function(nround, max.depth, min_child_weight, subsample, eta,
                         gamma,colsample_bytree,max_delta_step) {
    param<-list(booster = "gbtree",
                max_depth = max.depth,
                min_child_weight = min_child_weight,
                eta=eta,gamma=gamma,
                subsample = subsample,
                colsample_bytree = colsample_bytree,
                max_delta_step=max_delta_step,
                lambda = 1, alpha = 0,
                objective = "binary:logistic",
                eval_metric = "auc")
    cv <- xgb.cv(params = param, data = dtrain, folds = cv_folds,
                 nrounds = 1000, early_stopping_rounds = 10, maximize = TRUE,
                 verbose = verbose)
    list(Score = cv$evaluation_log$test_auc_mean[cv$best_iteration],
         Pred=cv$best_iteration)
    # we don't need cross-validation prediction and we need the number of rounds.
    # a workaround is to pass the number of rounds(best_iteration) to the Pred, which is a default parameter in the rbayesianoptimization library.
}

OPT_Res <- BayesianOptimization(xgb_cv_bayes,
                                bounds = list(max.depth =c(3L, 10L),
                                min_child_weight = c(1L, 40L),
                                subsample = c(0.6, 0.9),
                                eta=c(0.01,0.3),
                                gamma = c(0.0, 0.2),
                                colsample_bytree=c(0.5,0.8),
                                max_delta_step=c(1L,10L)),
                                init_grid_dt = NULL, init_points = 10,
                                n_iter = 10, acq = "ucb", kappa = 2.576,
                                eps = 0.0, verbose = verbose)

best_param <- list(
    booster = "gbtree",
    eval.metric = "auc",
    objective = "binary:logistic",
    max_depth = OPT_Res$Best_Par["max.depth"],
    eta = OPT_Res$Best_Par["eta"],
    gamma = OPT_Res$Best_Par["gamma"],
    subsample = OPT_Res$Best_Par["subsample"],
    colsample_bytree = OPT_Res$Best_Par["colsample_bytree"],
    min_child_weight = OPT_Res$Best_Par["min_child_weight"],
    max_delta_step = OPT_Res$Best_Par["max_delta_step"])

nrounds=OPT_Res$Pred[[which.max(OPT_Res$History$Value)]]
xgb_model <- xgb.train (params = best_param, data = dtrain, nrounds = nrounds)
```
