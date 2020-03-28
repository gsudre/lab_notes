# 2020-03-27 20:08:30

We are back to predicting classes. I'll run it through here because it's easier
to do feature engineering if needed. I'll script it out if I eventully need to
test multiple classifiers.

Most of Philip's directions came from a Slack DM on 03/27 around 5:50PM.

For groups, categ_inatt3 and categ_hi3; and the less important one is
categ_inatt2. The three groups are never affected, improvers or stable. So, it
should be easy-ish for the models to distinguish the never affected.
Classification between the other two groups might be too tough...

I think it's still worth trying the pairwise group predictions, eventually!

For the categ_inatt2- it's 4 groups (never affected, improvers, stable- but also
a large emergent/worsening group). I suspect that 4 group classification
might be pushing things too hard given the sample size, so maybe just focus on
the categ_inatt3 and the categ_hi3?

```r
fname = '~/data/baseline_prediction/prs_start/gf_philip_03272020.csv'
phen = 'categ_inatt3'

nfolds = 10
nreps = 10

data = read.csv(fname)
data$sex_numeric = as.factor(data$sex_numeric)
data$SES_group3 = as.factor(data$SES_group3)
data$population_self3 = as.factor(data$population_self3)
var_names = c(# PRS
              'ADHD_PRS0.000100.orig', 'ADHD_PRS0.001000.orig',
              'ADHD_PRS0.010000.orig', 'ADHD_PRS0.050000.orig',
              'ADHD_PRS0.100000.orig', 'ADHD_PRS0.200000.orig',
              'ADHD_PRS0.300000.orig', 'ADHD_PRS0.400000.orig',
              'ADHD_PRS0.500000.orig',
              # DTI
              'atr_fa', 'cst_fa', 'cing_cing_fa', 'cing_hipp_fa', 'cc_fa',
              'ilf_fa', 'slf_all', 'unc_fa',
              # demo
              'sex_numeric', 'SES_group3', 'population_self3',
              # cog
              'FSIQ', 'SS_RAW', 'DS_RAW', 'PS_RAW', 'VMI.beery_RAW',
              # anat
              'cerebellum_white', 'cerebellum_grey', 'amygdala',
              'cingulate', 'EstimatedTotalIntraCranialVol', 'lateral_PFC',
              'OFC', 'striatum', 'thalamus'
              )

covar_names = c(# DTI
                'norm.rot', 'norm.trans', # base_age, gender,
                # cog
                # base_age, gender
                # anat
                'average_qc', # age, gender (but not ICV)
                # PRS
                sapply(1:10, function(x) sprintf('PC%02d', x)) # age, gender
                # demo
                # just base_age (as gender one was of the targets)
                )

library(caret)
# I won't touch the predictors at first because the models that handle missing data don't care about that... I'll only dummify the factors
data2 = data[, c(var_names, covar_names)]
data2$phen = as.factor(data[, phen])
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
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=multiClassSummary)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('C5.0', 'rpart1SE', 'rpart2', 'rpartCost',
                                    #    'C5.0Tree', 'C5.0Rules'),
                                       'null'),
                        tuneList = NULL,
                        continue_on_fail = FALSE,
                        metric='AUC')

options(digits = 3)
train_results = c(model_list$null$results$AUC)
test_results = c()
for (m in names(model_list)) {
    if (m != 'null') {
        train_results = c(mean(model_list[[m]]$results$AUC), train_results)
    }
    preds_class = predict.train(model_list[[m]], newdata=X_test)
    preds_probs = predict.train(model_list[[m]], newdata=X_test, type='prob')
    dat = cbind(data.frame(obs = y_test, pred = preds_class), preds_probs)
    test_results = c(test_results,
                     multiClassSummary(dat, lev=colnames(preds_probs))['AUC'])
}
names(train_results) = names(model_list)
names(test_results) = names(model_list)

print(train_results)
print(test_results)
```
