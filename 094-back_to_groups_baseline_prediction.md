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
                        methodList = c('C5.0', 'rpart1SE', 'rpart2',
                                       'C5.0Tree', 'C5.0Rules',
                                       'AdaBoost.M1', 'AdaBag',
                                       'treebag', 'ada',
                                       'null'),
                        tuneList = NULL,
                        continue_on_fail = TRUE,
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

# 2020-03-28 07:25:18

Using all possible classifiers that can handle missing data crashed overnight,
so I changed the code to continue after fails. In any case, it's wise to get
ready for having to deal with imputations and multiple classifiers. I'll keep
the same framework though, so it's easy to do feature engineering later if
needed. Also, it'll be easy to come up with ensembles later if needed, probably
based on the best models.

Note that caretEnsemble doesn't work for multiclass, so we'd have to do it
manually. Or, do it only for the pairwise comparisons.

So, I created modelList_multiClass.R for that. Now it's just a matter of
scripting it:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_multiClass.R;
out_file=swarm.multiClass
rm $out_file
for clf in `cat multi_clf.txt`; do
    for sx in categ_inatt3 categ_hi3; do
        for imp in anat dti; do
            for cov in T F; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03272020.csv $sx $clf $imp 10 10 8 $cov ${my_dir}/multiClassAUC.csv;" >> $out_file;
            done;
        done;
    done;
done

swarm -g 20 -t 8 --job-name mcAUC --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

That failed because of a few reasons. I'm re-running it now after several fixed,
which include now selecting best in family to be the eldest and a family ID
change.

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_multiClass.R;
out_file=swarm.multiClass
rm $out_file
for clf in `cat multi_clf.txt`; do
    for sx in categ_inatt3 categ_hi3; do
        for imp in anat dti; do
            for cov in T F; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03282020.csv $sx $clf $imp 10 10 8 $cov ${my_dir}/multiClassEldestAUC.csv;" >> $out_file;
            done;
        done;
    done;
done

swarm -g 20 -t 8 --job-name mcAUC --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

And I might as well run the code for 2 classes over night. I'll keep itto the
classifiers that work for more than 2 classes to better evaluate across tests:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_twoClass.R;
out_file=swarm.twoClass
rm $out_file
for clf in `cat multi_clf.txt`; do
    for sx in categ_inatt3 categ_hi3; do
        for imp in anat dti; do
            for cov in T F; do
                for cs in "emerge_stable group_0_0" "emerge_stable improvers" \
                    "group_0_0 improvers"; do
                    echo "Rscript $my_script ${my_dir}/gf_philip_03282020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/twoClassEldestAUC.csv;" >> $out_file;
                done;
            done;
        done;
    done;
done

swarm -g 20 -t 8 --job-name tcAUC --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

While that's running, let's check again what are our best results for the
no-imputation models. They're basically static copies of the two and multi class
scripts I'm running, but forcing the models:

```r
library(caret)
library(caretEnsemble)
library(doParallel)

fname = '~/data/baseline_prediction/prs_start/gf_philip_03282020.csv'
phen = 'categ_inatt3'
nfolds = 10
nreps = 10
ncores = 16
use_covs = FALSE

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

if (use_covs) {
    data2 = data[, c(var_names, covar_names)]
} else {
    data2 = data[, var_names]
}
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
                       fam_rows[which.max(data[fam_rows, 'base_age'])])
    }
}
# data3 doesn't have the target column!
X_train <- data3[train_rows, ]
X_test <- data3[-train_rows, ]
y_train <- data2[train_rows,]$phen
y_test <- data2[-train_rows,]$phen
```

Up to there no much change. Now, we need to define the transformations that are
appropriate without imputations:

```r
# imputation and feature engineering
pp_order = c('zv', 'nzv', 'corr', 'YeoJohnson', 'center', 'scale')
pp = preProcess(X_train, method = pp_order)
X_train = predict(pp, X_train)
X_test = predict(pp, X_test)

registerDoParallel(ncores)
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
                        methodList = c('C5.0', 'rpart1SE', 'rpart2',
                                       'C5.0Tree', 'C5.0Rules',
                                       'AdaBoost.M1', 'AdaBag',
                                       'treebag',
                                       'null'),
                        tuneList = NULL,
                        continue_on_fail = TRUE,
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

print(phen)
print(use_covs)
print(train_results)
print(test_results)

# export fit
out_dir = '~/data/baseline_prediction/prs_start/multiClass/'
fname = sprintf('%s/noImp_modelList_%s_%s_%d_%d.RData',
                out_dir, phen, use_covs, nfolds, nreps)
save(model_list, file=fname)
```

I'm currently burning through allmy interactive cores, but here's the code to run
for the two class classification:

```r
library(caret)
library(caretEnsemble)
library(doParallel)

fname = '~/data/baseline_prediction/prs_start/gf_philip_03282020.csv'
nfolds = 10
nreps = 10
ncores = 5
use_covs = FALSE
phen = 'categ_hi3'
c1 = 'emerge_stable'
c2 = 'group_0_0'
# c1 = 'improvers'

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

if (use_covs) {
    data2 = data[, c(var_names, covar_names)]
} else {
    data2 = data[, var_names]
}
data2$phen = as.factor(data[, phen])
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# selecting only kids in the 2 specified groups
keep_me = data2$phen==c1 | data2$phen==c2
data3 = data3[keep_me, ]
data2 = data2[keep_me, ]
data = data[keep_me, ]

# split traing and test between members of the same family
train_rows = c()
for (fam in unique(data$FAMID)) {
    fam_rows = which(data$FAMID == fam)
    if (length(fam_rows) == 1) {
        train_rows = c(train_rows, fam_rows[1])
    } else {
        # choose the youngest kid in the family for training
        train_rows = c(train_rows,
                       fam_rows[which.max(data[fam_rows, 'base_age'])])
    }
}
# data3 doesn't have the target column!
X_train <- data3[train_rows, ]
X_test <- data3[-train_rows, ]
y_train <- factor(data2[train_rows,]$phen)
y_test <- factor(data2[-train_rows,]$phen)

# imputation and feature engineering
pp_order = c('zv', 'nzv', 'corr', 'YeoJohnson', 'center', 'scale')
pp = preProcess(X_train, method = pp_order)
X_train = predict(pp, X_train)
X_test = predict(pp, X_test)

registerDoParallel(ncores)
getDoParWorkers()
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=twoClassSummary)

model_list <- caretList(X_train,
                        y_train,
                        trControl = fitControl,
                        methodList = c('C5.0', 'rpart1SE', 'rpart2',
                                       'C5.0Tree', 'C5.0Rules',
                                       'AdaBoost.M1', 'AdaBag',
                                       'treebag', 'ada',
                                       'null'),
                        tuneList = NULL,
                        continue_on_fail = TRUE,
                        metric='ROC')

options(digits = 3)
train_results = c(model_list$null$results$ROC)
test_results = c()
for (m in names(model_list)) {
    if (m != 'null') {
        train_results = c(mean(model_list[[m]]$results$ROC), train_results)
    }
    preds_class = predict.train(model_list[[m]], newdata=X_test)
    preds_probs = predict.train(model_list[[m]], newdata=X_test, type='prob')
    dat = cbind(data.frame(obs = y_test, pred = preds_class), preds_probs)
    test_results = c(test_results,
                     twoClassSummary(dat, lev=colnames(preds_probs))['ROC'])
}
names(train_results) = names(model_list)
names(test_results) = names(model_list)

print(c1)
print(c2)
print(phen)
print(use_covs)
print(train_results)
print(test_results)

# export fit
out_dir = '~/data/baseline_prediction/prs_start/twoClass/'
fname = sprintf('%s/noImp_modelList_%s_%s_%s_%s_%d_%d.RData',
                out_dir, phen, c1, c2, use_covs, nfolds, nreps)
save(model_list, file=fname)
```

# 2020-03-29 07:40:57

For the 3-class non-imputation results, I get:

```
> print(phen)
[1] "categ_inatt3"
> print(use_covs)
[1] FALSE
> print(train_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.607       0.612       0.610       0.530       0.550       0.537       0.539       0.562       0.500 
> print(test_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.625       0.595       0.574       0.496       0.502       0.612       0.648       0.589       0.500 
> print(phen)
[1] "categ_inatt3"
> print(use_covs)
[1] TRUE
> print(train_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.587       0.607       0.594       0.512       0.525       0.535       0.545       0.547       0.500 
> print(test_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.623       0.556       0.561       0.566       0.492       0.590       0.608       0.571       0.500 
> print(phen)
[1] "categ_hi3"
> print(use_covs)
[1] FALSE
> print(train_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.590       0.661       0.640       0.573       0.596       0.577       0.568       0.592       0.500 
> print(test_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.516       0.500       0.501       0.519       0.489       0.573       0.583       0.508       0.500 
> print(phen)
[1] "categ_hi3"
> print(use_covs)
[1] TRUE
> print(train_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.595       0.664       0.654       0.579       0.589       0.564       0.551       0.606       0.500 
> print(test_results)
       C5.0    rpart1SE      rpart2    C5.0Tree   C5.0Rules AdaBoost.M1      AdaBag     treebag        null 
      0.546       0.463       0.511       0.528       0.507       0.573       0.578       0.506       0.500 
```

And let's start evaluating the compiled results in twoClassEldestAUC.csv and
multiClassEldestAUC.csv.

Philip also sent out a new gf today with two new categories combining the inatt
and hi:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_multiClass.R;
out_file=swarm.multiClassComb
rm $out_file
for clf in `cat multi_clf.txt`; do
    for sx in categ_all.3 categ_all.4; do
        for imp in anat dti; do
            for cov in T F; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $clf $imp 10 10 8 $cov ${my_dir}/multiClassCombEldestAUC.csv;" >> $out_file;
            done;
        done;
    done;
done

swarm -g 20 -t 8 --job-name mcCombAUC --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

And similar idea for 2 class, but because classes are different I'll need to
split the loop:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_twoClass.R;
out_file=swarm.twoClassComb
rm $out_file
for clf in `cat multi_clf.txt`; do
    for imp in anat dti; do
        for cov in T F; do
            sx="categ_all.3";
            for cs in "improvers never_affected" "improvers symptomatic" \
                "never_affected symptomatic"; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/twoClassCombEldestAUC.csv;" >> $out_file;
            done;
            sx="categ_all.4";
            for cs in "emergent improvers" "emergent never_affected" \
                "emergent stable_symptomatic" "improvers never_affected" \
                "improvers stable_symptomatic" "never_affected stable_symptomatic"; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/twoClassCombEldestAUC.csv;" >> $out_file;
            done;
        done;
    done;
done

swarm -g 20 -t 8 --job-name tcCombAUC --time 4:00:00 -f $out_file \
    -m R --partition norm --logdir trash
```

# 2020-03-29 20:36:00

I forgot to run categ_inatt2, which has 4 classes. Here is goes:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_multiClass.R;
out_file=swarm.multiClassComb
rm $out_file
for clf in `cat multi_clf.txt`; do
    sx='categ_inatt2';
    for imp in anat dti; do
        for cov in T F; do
            echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $clf $imp 10 10 8 $cov ${my_dir}/multiClassInatt2EldestAUC.csv;" >> $out_file;
        done;
    done;
done

swarm -g 20 -t 8 --job-name mc2AUC --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

And similar idea for 2 class, but we have 4 classes:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/modelList_twoClass.R;
out_file=swarm.twoClassComb
rm $out_file
for clf in `cat multi_clf.txt`; do
    for imp in anat dti; do
        for cov in T F; do
            sx="categ_inatt2";
            for cs in "emergent improvers" "emergent group_0_0" \
                "emergent group_2_2" "improvers group_0_0" \
                "improvers group_2_2" "group_0_0 group_2_2"; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/twoClassInatt2EldestAUC.csv;" >> $out_file;
            done;
        done;
    done;
done

swarm -g 20 -t 8 --job-name tc2AUC --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

And while we wait for these to run, let's compile the comb results...

# 2020-03-30 08:01:12

I just noticed that I was reporting the results for the mean over parameters,
instead of the max, which would be the max over parameters. The other option is
to report the mean and other parameters for the resampling, which uses the best
parameter anyways. Let's do that.

The function do to the multiclass is ready, but it'll need to be swarmed:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/resample_multiClass.R;
out_file=swarm.resample_mc
rm $out_file
for clf in `cat multi_clf.txt`; do
    for sx in categ_hi3 categ_inatt3 categ_inatt2 categ_all.3 categ_all.4; do
        for imp in anat dti; do
            for cov in T F; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $clf $imp 10 10 8 $cov ${my_dir}/resamp_multiClassEldestAUC.csv;" >> $out_file;
            done;
        done;
    done;
done

swarm -g 10 -t 1 --job-name resampMC --time 20:00 -f $out_file \
    -m R --partition quick --logdir trash
```

We had lots of racing conditions... let me see if I run parallel it wouldn't be
better.

```bash
cat swarm.resample_mc | parallel --max-args=1 -j 32 {1};
```

Now, time to do the same thing for the 2-class situation:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/resample_twoClass.R;
out_file=swarm.resample_tc
rm $out_file
for clf in `cat multi_clf.txt`; do
    for imp in anat dti; do
        for cov in T F; do
            sx="categ_inatt2";
            for cs in "emergent improvers" "emergent group_0_0" \
                "emergent group_2_2" "improvers group_0_0" \
                "improvers group_2_2" "group_0_0 group_2_2"; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/resamp_twoClassEldestAUC.csv;" >> $out_file;
            done;
            sx="categ_all.3";
            for cs in "improvers never_affected" "improvers symptomatic" \
                "never_affected symptomatic"; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/resamp_twoClassEldestAUC.csv;" >> $out_file;
            done;
            sx="categ_all.4";
            for cs in "emergent improvers" "emergent never_affected" \
                "emergent stable_symptomatic" "improvers never_affected" \
                "improvers stable_symptomatic" "never_affected stable_symptomatic"; do
                echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/resamp_twoClassEldestAUC.csv;" >> $out_file;
            done;
            for sx in categ_inatt3 categ_hi3; do
                for cs in "emerge_stable group_0_0" "emerge_stable improvers" \
                    "group_0_0 improvers"; do
                    echo "Rscript $my_script ${my_dir}/gf_philip_03292020.csv $sx $cs $clf $imp 10 10 8 $cov ${my_dir}/resamp_twoClassEldestAUC.csv;" >> $out_file;
                done;
            done;
        done;
    done;
done
```

Just had a chat with PHilip and we should focus on the all class, and the two
group comparisons. Pick 2-3 models but only one for main text, other for
comparison in supplement. Focus on test set but make sure it's OK for training. 


# TODO
* investigate these two class results:
  categ_hi3 emerge_stable   improvers   bagEarthGCV dti FALSE   10  10  0.371875
  NA  0.5 0.84    0.5
* what if I report the max and median AUC in training?
* need to run non-impute machines in the combined categories and inatt2
* re-run the no-imputation models without transformations to the predictors?
* maybe report MCC and/or F1?