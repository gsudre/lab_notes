args = commandArgs(trailingOnly=TRUE)
if (length(args) != 3) {
  stop("Arguments: seed[-1] pct_var root_fname", call.=FALSE)
} else {
  if (args[1] == '-1') {
    myseed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  } else {
    myseed = as.integer(args[1])
  }
  pct_var = as.numeric(args[2])
  root_fname = args[3]
}

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
cl <- makeCluster(8)
X_resid = parSapply(cl, X, get_needed_residuals, 'y ~ merged$age_at_scan + I(merged$age_at_scan^2) + merged$SEX', .1, merged[keep_me, ])
stopCluster(cl)
X_resid = as.data.frame(X_resid)

set.seed(myseed)
split <- createDataPartition(y, p = .7, list = FALSE)
Xtrain <- X_resid[ split, ]
ytrain <- y[ split ]
Xtest  <- X_resid[-split, ]
ytest = y[-split]

library(parallel)
cl <- makeCluster(8)
pvals = parSapply(cl, Xtrain, function(d, ytrain) t.test(d ~ ytrain)$p.value, ytrain)
stopCluster(cl)
Xtrain = Xtrain[, which(pvals <= .05)]
print(dim(Xtrain))

keep_me = sapply(colnames(Xtrain), function(d) which(colnames(Xtest) == d))
Xtest = Xtest[, keep_me]

pp <- preProcess(Xtrain, method = c('BoxCox', 'center', 'scale', 'pca'), thresh=pct_var)
filtXtrain<- predict(pp, Xtrain)
filtXtest <- predict(pp, Xtest)
print(dim(filtXtrain))

tuneLength=10
set.seed(myseed)
index <- createResample(ytrain, times=20)

set.seed(myseed)
fullCtrl <- trainControl(method = "boot",
                         index = index,
                         savePredictions="final",
                         classProbs=TRUE,
                         summaryFunction=twoClassSummary)

require(doMC)
registerDoMC(cores=8)
library(caretEnsemble)
model_list <- caretList(
  filtXtrain, ytrain,
  tuneLength=10,
  trControl=fullCtrl,
  metric='ROC',
  methodList=c('rf', 'kernelpls', 'svmRadial', 'knn', 'rpart', 'bagEarthGCV', 'LogitBoost', 'lda', 'nb')
  )

greedy_ensemble <- caretEnsemble(
  model_list,
  metric='ROC',
  trControl=trainControl(
    number=2,
    summaryFunction=twoClassSummary,
    classProbs=TRUE
    ))
# ROC stats
summary(greedy_ensemble)
model_preds <- lapply(model_list, function(d, x, y) eval_model(d, x, y, c('ADHD', 'NV'))['ROC'], filtXtest, ytest)
model_preds <- data.frame(model_preds)
print(model_preds)
# sensitivity stats
model_preds <- lapply(model_list, function(d, x, y) eval_model(d, x, y, c('ADHD', 'NV'))['Sens'], filtXtest, ytest)
model_preds <- data.frame(model_preds)
print(model_preds)
# specificity stats
model_preds <- lapply(model_list, function(d, x, y) eval_model(d, x, y, c('ADHD', 'NV'))['Spec'], filtXtest, ytest)
model_preds <- data.frame(model_preds)
print(model_preds)
# Accuracy stats
model_preds <- lapply(model_list, predict, newdata=filtXtest)
model_preds <- lapply(model_preds, function(d, obs) postResample(d, obs)[1], obs=ytest)
model_preds <- data.frame(model_preds)
print(model_preds)
print(sprintf('No information rate: Accuracy=%f', max(table(ytrain)/length(ytrain))))
model_preds <- lapply(model_list, predict, newdata=filtXtest)
model_preds <- lapply(model_preds, function(d) confusionMatrix(d, ref=ytest)$overall['AccuracyPValue'])
model_preds <- data.frame(model_preds)
print(model_preds)

sink()

fname = sprintf('%s_%d.RData', root_fname, myseed)
save(model_list, myseed, pp, file=fname)
