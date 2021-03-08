args = commandArgs(trailingOnly=TRUE)
s = as.integer(args[1])
root_fname = args[2]
njobs = 6
tuneLength=10
myseed = 1234
###########

fname = sprintf('%s_%d.log', root_fname, myseed)
sink(fname, append=FALSE, split=TRUE)

source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')
tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
load('~/data/baseline_prediction/dti/ad_voxelwise.RData')
dti_vdata = cbind(tract_data$maskid, ad_data)
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf_base = gf[gf$BASELINE=='BASELINE', ]
my_ids = intersect(gf_base$MRN, tract_data$MRN)
merged = mergeOnClosestDate(gf_base, tract_data, my_ids)
rm_me = abs(merged$dateX.minus.dateY.months) > 12
merged = merged[!rm_me, ]
dti_base_vdata = merge(merged$maskid, dti_vdata, by.x=1, by.y=1, all.y=F, all.x=T)

get_needed_residuals = function(y, fm_str, cutoff, merged) {
  fm = as.formula(fm_str)
  fit = lm(fm)
  # selecting which covariates to use
  fm = "y ~ "
  for (r in 2:dim(summary(fit)$coefficients)[1]) {
    if (summary(fit)$coefficients[r, 4] < cutoff) {
      cname = rownames(summary(fit)$coefficients)[r]
      cname = gsub("SEXMale", "SEX", cname)
      fm = sprintf('%s + %s', fm, cname)
    }
  }
  # don't do anything if no variables were significant
  if (fm == 'y ~ ') {
    return(y)
  } else {
    opt_fit = lm(as.formula(fm))
    return(opt_fit$residuals)
  }
}

X = dti_base_vdata[, 2:ncol(dti_base_vdata)]
rm_me = colSums(is.na(X)) > 0
X = X[, !rm_me]
keep_me = merged$age <= 12
X = X[keep_me, ]
y = merged$DX_BASELINE
y[y!='NV'] = 'ADHD'
y = factor(y, levels=c('ADHD', 'NV'))
y = y[keep_me]
library(parallel)
cl <- makeCluster(njobs)
X_resid = parSapply(cl, X, get_needed_residuals, 'y ~ merged$age_at_scan + I(merged$age_at_scan^2) + merged$SEX', .1, merged[keep_me, ])
stopCluster(cl)
X_resid = as.data.frame(X_resid)

# we need to do the univariate filter first otherwise we will be doing pca
# on garbage
print(sprintf('LO %d / %d (%s)', s, length(y), y[s]))
Xtrain <- X_resid[ -s, ]
ytrain <- y[ -s ]
Xtest  <- X_resid[s, ]

library(parallel)
cl <- makeCluster(njobs)
pvals = parSapply(cl, Xtrain, function(d, ytrain) t.test(d ~ ytrain)$p.value, ytrain)
stopCluster(cl)
Xtrain = Xtrain[, which(pvals <= .05)]
keep_me = sapply(colnames(Xtrain), function(d) which(colnames(Xtest) == d))
Xtest = Xtest[, keep_me]

pp <- preProcess(rbind(Xtrain, Xtest),
                 method = c('BoxCox', 'center', 'scale', 'pca'), thresh=.6)
filtXtrain <- predict(pp, Xtrain)
filtXtest <- predict(pp, Xtest)

set.seed(myseed)
index <- createMultiFolds(ytrain, k = 10, times = 10)

fullCtrl <- trainControl(method = "repeatedcv",
                       index = index,
                       savePredictions="final",
                       classProbs=TRUE)

require(doMC)
registerDoMC(cores=njobs)
library(caretEnsemble)

set.seed(myseed)
model_list <- caretList(
filtXtrain, ytrain,
tuneLength=10,
trControl=fullCtrl,
metric='Accuracy',
methodList=c('kernelpls', 'bagEarthGCV', 'knn')
)

greedy_ensemble <- caretEnsemble(
  model_list,
  metric='Accuracy',
  trControl=trainControl(
    number=2,
    classProbs=TRUE
    ))
# ROC stats
print(summary(greedy_ensemble))

preds = lapply(model_list, predict, newdata=filtXtest, type='prob')
print(do.call(rbind, preds))

sink()
