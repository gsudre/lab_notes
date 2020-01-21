# 2020-01-21 11:10:11

At this point it looks like adding PRS doesn't help at att in predicting
improvement. base_sx is still the best predictor. Let's redo here our best
results, for both HI and inatt:

```r
setwd('~/data/baseline_prediction/prs_start/')
data = read.csv('long_clin_01062020_lt16.csv')

winsorize = function(x, cut = 0.01){
  cut_point_top <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x, cut, na.rm = T)
  i = which(x >= cut_point_top) 
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom) 
  x[j] = cut_point_bottom
  return(x)
}

df = data.frame(MRN=unique(data$MRN))
for (r in 1:nrow(df)) {
    subj_data = data[data$MRN==df$MRN[r], ]
    for (sx in c('inatt', 'hi')) {
        fit = lm(as.formula(sprintf('SX_%s ~ age', sx)), data=subj_data)
        df[r, sprintf('slope_%s', sx)] = fit$coefficients['age']
        base_row = which.min(subj_data$age)
        df[r, sprintf('base_%s', sx)] = subj_data[base_row, sprintf('SX_%s', sx)]
        last_row = which.max(subj_data$age)
        df[r, sprintf('last_%s', sx)] = subj_data[last_row, sprintf('SX_%s', sx)]
        df[r, 'base_age'] = subj_data[base_row, 'age']
        df[r, 'last_age'] = subj_data[last_row, 'age']
        df[r, 'sex'] = subj_data[last_row, 'sex']
    }
}
idx = df$base_inatt>=3 | df$base_hi>=3
for (sx in c('inatt', 'hi')) {
    df[, sprintf('slope_%s_GE3_wp05', sx)] = NA
    junk = winsorize(df[idx, sprintf('slope_%s', sx)], cut=.05)
    df[idx, sprintf('slope_%s_GE3_wp05', sx)] = junk
}
prs = read.csv('/Volumes/NCR/reference/merged_NCR_1KG_PRS_12192019.csv')
data = merge(df, prs, by='MRN', all.x=F, all.y=F)

var_names = colnames(data)[grepl(colnames(data), pattern='ADHD_') |
                           grepl(colnames(data), pattern='PC')]
data$y = NA
data[which(data$slope_hi_GE3_wp05 < -.5),]$y = 'imp'
data[which(data$slope_hi_GE3_wp05 >= -.5),]$y = 'nonimp'
use_me = !is.na(data[,]$y)
X1 = prcomp(data[use_me, var_names[1:12]])$x
colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
X = cbind(X1, as.matrix(data[use_me, var_names[13:22]]),
          as.matrix(data[use_me, c('base_age', 'base_hi')]))
X = scale(X)
colnames(X)[(ncol(X)-1):ncol(X)] = c('base_age', 'base_sx')
X = cbind(X, data[use_me, 'sex'])
colnames(X)[ncol(X)] = 'sex'
Y = data[use_me, ]$y

library(caret)
set.seed(3456)
trainIndex <- createDataPartition(Y, p = .9, list = FALSE, times = 1)
X_train <- X[ trainIndex,]
X_test  <- X[-trainIndex,]
y_train <- Y[trainIndex]
y_test  <- Y[-trainIndex]

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", nfold = 10,
                      type.measure = "auc", alpha = 1, parallel=T)
md_lasso <- glmnet(X_train, y_train, family = "binomial", lambda = cv_lasso$lambda.1se,
                   alpha = 1)
roc(y_test, as.numeric(predict(md_lasso, X_test, type = "response")))
```

To recap, it's .7024 for hi and .6042 for inatt.

Adding raw genomics didn't help much either. Now, let's sprinkle in some
neuroimaging to see if those do any good:

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# m has all scans for everyone who has PRS data
m = merge(m, data, by.x="Medical.Record...MRN...Subjects", by.y='MRN',
          all.x=F, all.y=T)

# restrict based on QC
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
m = m[m$"age_at_scan...Scan...Subjects" < 18, ]
m = m[m$"goodVolumes" <= 61, ]
m = m[m$"numVolumes" < 80, ]
```

So, 107 of the 1035 scans are NaNs (10%). Another way to see it is 107 of the
original 393 subjects don't have scan data. 

We have a few options:

1) use only classifiers that can handle missing data
2) implement a voting scheme within modality
3) reduce the data to only subjects that have all datatypes
4) impute (within training data or across all of it)

But we can decide this later. The actual outlier detection can be done only
within the scans that we have, without imputation:

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# keep only scans for subjects that have PRS
data = merge(df, prs, by='MRN', all.x=F, all.y=F)
keep_me = sapply(1:nrow(m), function(x) m$Medical.Record...MRN...Subjects[x] %in% data$MRN)
m = m[keep_me, ]

# restrict based on QC
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
m = m[m$"age_at_scan...Scan...Subjects" < 18, ]
m = m[m$"goodVolumes" <= 61, ]
m = m[m$"numVolumes" < 80, ]

# down to 928 scans that obey criteria and have PRS
library(solitude)
iso <- isolationForest$new()
iso$fit(m[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]
library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(m[, qc_vars], k = round(.5 * nrow(m)))

qtile=.95
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

tracts = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
# somehow I have two entries for 1418?
x = duplicated(tracts$id)
data = merge(m[idx,], tracts[!x, ], by.x='Mask.ID...Scan', by.y='id')
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])

iso <- isolationForest$new()
iso$fit(data[, tract_names])
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(data[, tract_names], k = round(.5 * nrow(data)))

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

data = data[idx, ]

# down to 801 scans when only scans at .95 in both criteria are used

# selecting earliest scan for each subject, regardless of score, as we're assu ing every scan now is good
keep_me = c()
for (s in unique(data$Medical.Record...MRN...Subjects)) {
    subj_rows = which(data$Medical.Record...MRN...Subjects == s)
    subj_data = data[subj_rows, ]
    min_subj_row = which.min(subj_data$age_at_scan...Scan...Subjects)
    keep_me = c(keep_me, subj_rows[min_subj_row])
}
dprs = merge(df, prs, by='MRN', all.x=F, all.y=F)
ddti = data[keep_me, c(tract_names, qc_vars,
                       'age_at_scan...Scan...Subjects',
                       'Medical.Record...MRN...Subjects')]
```

OK, now let's combine the data and see what's the best way to approach this
missing data issue:

```r
data = merge(dprs, ddti, by.x='MRN', by.y='Medical.Record...MRN...Subjects',
             all.x=T, all.y=F)
```

After the cleaning cut-offs, we have 114 out of the 393 subjects without brain
data. That's right at 29%... let's see how the classification remains if we dump
everyone without data, and then if we impute them. Imputation here is somewhat
unkosher as we'll do it on the entire dataset, but I could also do it within the
training set later if it's an interesting approach.

```r
var_names = colnames(data)[grepl(colnames(data), pattern='ADHD_') |
                           grepl(colnames(data), pattern='PC')]
data$y = NA
data[which(data$slope_hi_GE3_wp05 < -.5),]$y = 'imp'
data[which(data$slope_hi_GE3_wp05 >= -.5),]$y = 'nonimp'
use_me = !is.na(data[,]$y) & !is.na(data$meanX.rot)
X1 = prcomp(data[use_me, var_names[1:12]])$x
colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
X = cbind(X1, as.matrix(data[use_me, var_names[13:22]]),
          as.matrix(data[use_me, c('base_age', 'base_hi')]))
X = scale(X)
colnames(X)[(ncol(X)-1):ncol(X)] = c('base_age', 'base_sx')
X = cbind(X, data[use_me, 'sex'])
colnames(X)[ncol(X)] = 'sex'
Y = data[use_me, ]$y

library(caret)
set.seed(3456)
trainIndex <- createDataPartition(Y, p = .9, list = FALSE, times = 1)
X_train <- X[ trainIndex,]
X_test  <- X[-trainIndex,]
y_train <- Y[trainIndex]
y_test  <- Y[-trainIndex]

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", nfold = 10,
                      type.measure = "auc", alpha = 1, parallel=T)
md_lasso <- glmnet(X_train, y_train, family = "binomial", lambda = cv_lasso$lambda.1se,
                   alpha = 1)
roc(y_test, as.numeric(predict(md_lasso, X_test, type = "response")))
```

We're dealing with 166 ADHD subjects now. LASSO results are actually quite good,
at .83, but as usual only using base_sx. That's all for HI. Does it get better
by adding DTI?

```r
var_names = c(colnames(data)[grepl(colnames(data), pattern='ADHD_') |
                           grepl(colnames(data), pattern='PC')],
              tract_names, qc_vars, 'age_at_scan...Scan...Subjects')
data$y = NA
data[which(data$slope_hi_GE3_wp05 < -.5),]$y = 'imp'
data[which(data$slope_hi_GE3_wp05 >= -.5),]$y = 'nonimp'
use_me = !is.na(data[,]$y) & !is.na(data$meanX.rot)
X1 = prcomp(data[use_me, var_names[1:12]])$x
colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
X = cbind(X1, as.matrix(data[use_me, var_names[13:70]]),
          as.matrix(data[use_me, c('base_age', 'base_hi')]))
X = scale(X)
colnames(X)[(ncol(X)-1):ncol(X)] = c('base_age', 'base_sx')
X = cbind(X, data[use_me, 'sex'])
colnames(X)[ncol(X)] = 'sex'
Y = data[use_me, ]$y

library(caret)
set.seed(3456)
trainIndex <- createDataPartition(Y, p = .9, list = FALSE, times = 1)
X_train <- X[ trainIndex,]
X_test  <- X[-trainIndex,]
y_train <- Y[trainIndex]
y_test  <- Y[-trainIndex]

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", nfold = 10,
                      type.measure = "auc", alpha = 1, parallel=T)
md_lasso <- glmnet(X_train, y_train, family = "binomial", lambda = cv_lasso$lambda.1se,
                   alpha = 1)
roc(y_test, as.numeric(predict(md_lasso, X_test, type = "response")))
```

Same issue, because no DTI variable gets selected... what if we upsample it?

```r
var_names = c(colnames(data)[grepl(colnames(data), pattern='ADHD_') |
                           grepl(colnames(data), pattern='PC')],
              tract_names, qc_vars, 'age_at_scan...Scan...Subjects')
data$y = NA
data[which(data$slope_hi_GE3_wp05 < -.5),]$y = 'imp'
data[which(data$slope_hi_GE3_wp05 >= -.5),]$y = 'nonimp'
use_me = !is.na(data[,]$y) & !is.na(data$meanX.rot)
X1 = prcomp(data[use_me, var_names[1:12]])$x
colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
X = cbind(X1, as.matrix(data[use_me, var_names[13:70]]),
          as.matrix(data[use_me, c('base_age', 'base_hi')]))
X = scale(X)
colnames(X)[(ncol(X)-1):ncol(X)] = c('base_age', 'base_sx')
X = cbind(X, data[use_me, 'sex'])
colnames(X)[ncol(X)] = 'sex'
Y = data[use_me, ]$y

library(caret)
set.seed(3456)
trainIndex <- createDataPartition(Y, p = .9, list = FALSE, times = 1)
X_train <- X[ trainIndex,]
X_test  <- X[-trainIndex,]
y_train <- Y[trainIndex]
y_test  <- Y[-trainIndex]

set.seed(3456)
up_train <- upSample(x = X_train,
                     y = as.factor(y_train))
X_train = up_train[, -ncol(up_train)]
y_train = as.character(up_train$Class)

pkgs <- list("glmnet", "doParallel", "foreach", "pROC")
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)

cv_lasso <- cv.glmnet(X_train, y_train, family = "binomial", nfold = 10,
                      type.measure = "auc", alpha = 1, parallel=T)
md_lasso <- glmnet(X_train, y_train, family = "binomial", lambda = cv_lasso$lambda.1se,
                   alpha = 1)
roc(y_test, as.numeric(predict(md_lasso, X_test, type = "response")))
```

Nope, same result. 

# TODO
* enet?
* inatt?
* other fold sizes?
* voxelwise?

