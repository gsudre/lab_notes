# 2020-04-02 20:35:35

Philip said I could try some ML here, just to see what we get. Why not?

I'll start within tissue, and then potentially do something to combine them?
Note that here we have WAY more genes than data, so we'll need to be very
careful with overfitting. Some clever cross-validation might be needed as well. 

At first, our data dimensions are 35,924 genes for 56 subjects,
and that's without taking into consideration the covariates to test.

How about doing LOOCV after leaving out one subject. That will result in nsubjs
models, and we can just average the varImp across all of them? Also, let's focus
on models with built-in cross-validation at first, and add some feature
engineering later if needed (likely).

```r
library(caret)

myregion = 'ACC'

data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data$substance_group = as.factor(data$substance_group)
data$batch = as.factor(data$batch)
# no column names as numbers!
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names

pop_code = read.csv('~/data/rnaseq_derek/file_pop.csv')
data = merge(data, pop_code, by='hbcc_brain_id')
data = data[data$Region==myregion, ]

# variables to be tested/screened
covar_names = c(# brain-related
                "bainbank", 'PMI', 'Manner.of.Death',
                # technical
                'batch', 'RINe',
                # clinical
                'comorbid_group', 'substance_group',
                # others
                'Sex', 'Age', 'POP_CODE')

# only covariates can be binary, and I'm getting stack overflow errors sending
# everything to dummyvars...
data2 = data[, c(covar_names, 'Diagnosis')]
dummies = dummyVars(Diagnosis ~ ., data = data2)
data3 = predict(dummies, newdata = data2)
# remove linear combination variables
comboInfo <- findLinearCombos(data3)
data3 = data3[, -comboInfo$remove]
data4 = cbind(data[, grex_names], data3)

# I'll go ahead and do the pre-processing here because it'll be very costly to
# to it inside LOOCV
set.seed(42)
# data4 doesn't even have Diagnosis, and no NAs
pp_order = c('zv', 'nzv', 'center', 'scale')
pp = preProcess(data4, method = pp_order)
X = predict(pp, data4)
```

OK, at this point I saved the data because the preprocessing took longer than I
wanted. It's not crucial, but it might save a few minutes in the future.

```r
library(caret)
library(doParallel)

X = readRDS('~/data/rnaseq_derek/X_ACCnoPH_zv_nzv_center_scale.rds')
myregion = 'ACC'
ncores = 32
nfolds = 10
nreps = 10
clf_model = 'kernelpls'
just_target = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
y = just_target[just_target$Region==myregion, 'Diagnosis']

my_summary = function(data, lev = NULL, model = NULL) {
    tcs = twoClassSummary(data, lev=lev)
    a = c((tcs['Sens'] + tcs['Spec'])/2, tcs)
    names(a)[1] = 'BalancedAccuracy'
    return(a)
}

registerDoParallel(ncores)
getDoParWorkers()

set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = nfolds,
                           repeats = nreps,
                           savePredictions = 'final',
                           allowParallel = TRUE,
                           classProbs = TRUE,
                           summaryFunction=my_summary)

# let's then do a repeated 10 fold CV within LOOCV. We save the test predictions
# to later compute the overall result.
varimps = X
varimps[,] = NA
test_preds = c()
for (test_rows in 1:length(y)) {
    print(sprintf('Hold out %d / %d', train_rows, length(y)))
    X_train <- X[-test_rows, ]
    X_test <- X[test_rows, ]
    y_train <- factor(y[-test_rows])
    y_test <- factor(y[test_rows])

    set.seed(42)
    fit <- train(X_train, y_train, trControl = fitControl, method = clf_model,
                 metric='BalancedAccuracy')

    preds_class = predict.train(fit, newdata=X_test)
    preds_probs = predict.train(fit, newdata=X_test, type='prob')
    dat = cbind(data.frame(obs = y_test, pred = preds_class), preds_probs)
    test_preds = rbind(test_preds, dat)

    tmp = varImp(fit, useModel=T)$importance
    varimps[test_rows,] = tmp
}


fit = model_list[[clf_model]]
resamps = resamples(list(fit=fit, tmp=fit))
bacc_stats = summary(resamps)$statistics$BalancedAccuracy['fit',]
cnames = sapply(names(bacc_stats), function(x) sprintf('BalancedAccuracy_%s', x))
names(bacc_stats) = cnames
auc_stats = summary(resamps)$statistics$ROC['fit',]
cnames = sapply(names(auc_stats), function(x) sprintf('AUC_%s', x))
names(auc_stats) = cnames
sens_stats = summary(resamps)$statistics$Sens['fit',]
cnames = sapply(names(sens_stats), function(x) sprintf('Sens_%s', x))
names(sens_stats) = cnames
spec_stats = summary(resamps)$statistics$Spec['fit',]
cnames = sapply(names(spec_stats), function(x) sprintf('Spec_%s', x))
names(spec_stats) = cnames

preds_class = predict.train(fit, newdata=X_test)
preds_probs = predict.train(fit, newdata=X_test, type='prob')
dat = cbind(data.frame(obs = y_test, pred = preds_class), preds_probs)
mcs = my_summary(dat, lev=colnames(preds_probs))
test_results = c(mcs['BalancedAccuracy'], mcs['ROC'], mcs['Sens'], mcs['Spec'])
names(test_results) = c('test_BalancedAccuracy', 'test_AUC', 'test_Sens',
                        'test_Spec')

res = c(phen, clf_model, c1, c2, impute, use_covs, nfolds, nreps,
        auc_stats, sens_stats, spec_stats, bacc_stats, test_results)
line_res = paste(res,collapse=',')
write(line_res, file=out_file, append=TRUE)
print(line_res)

# export variable importance
a = varImp(fit, useModel=T)
b = varImp(fit, useModel=F)
out_dir = '~/data/baseline_prediction/prs_start/twoClassBA/'
fname = sprintf('%s/varimp_%s_%s_%s_%s_%s_%s_%d_%d.csv',
                out_dir, clf_model, phen, c1, c2, impute, use_covs, nfolds, nreps)
# careful here because for non-linear models the rows of the importance matrix
# are not aligned!!!
write.csv(cbind(a$importance, b$importance), file=fname)

# export fit
fname = sprintf('%s/fit_%s_%s_%s_%s_%s_%s_%d_%d.RData',
                out_dir, clf_model, phen, c1, c2, impute, use_covs, nfolds, nreps)
save(fit, file=fname)
```

# TODO
* check for linearities again after preprocessing?
* might need to do most of the prepocessing and save the file for later? Also,
  might need some extra visualizations of the data and do it outside of caret.
  Especially for correlations, PCA, etc. But maybe just try some scaling and
  range at first to run some models with built in feature selection, and go from
  there!
* this would be a good task for 2vs2 decoding, since all we want to know is
  which genes to best in differnetiating the par of DX... worth trying 
* I had to remove the pH variable because it had 23 NAs, and it was the only
  variable with NAs. So, probably good to test if our results change at all with
  it later.

