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

X = cbind(merged[keep_me, c('age_at_scan', 'SEX')], X)

Xtrain <- X[ -s, ]
ytrain <- y[ -s ]
Xtest  <- X[s, ]

pp <- preProcess(Xtrain, method = c('BoxCox', 'center', 'scale'))
filtXtrain <- predict(pp, Xtrain)
filtXtest <- predict(pp, Xtest)

library(parallel)
cl <- makeCluster(njobs)
pvals = parSapply(cl, filtXtrain,
                  function(d, ytrain) kruskal.test(d ~ ytrain)$p.value, ytrain)
stopCluster(cl)
filtXtrain = filtXtrain[, which(pvals <= .05)]
keep_me = sapply(colnames(filtXtrain), function(d) which(colnames(filtXtest) == d))
filtXtest = filtXtest[, keep_me]

pp <- preProcess(filtXtrain, method = c('pca'))
filtXtrain <- predict(pp, filtXtrain)
filtXtest <- predict(pp, filtXtest)

print(sprintf('LO %d / %d (%s)', s, length(y), y[s]))

set.seed(myseed)
index <- createResample(ytrain, times = 200)

require(doMC)
registerDoMC(cores=njobs)
set.seed(myseed)
fullCtrl <- trainControl(method = "boot",
                       number=200,
                       index = index,
                       savePredictions="all",
                       classProbs=TRUE,
                       returnResamp = 'all',
                       indexFinal=createResample(ytrain, times=1)[[1]])

library(caretEnsemble)
set.seed(myseed)
model_list <- caretList(
    filtXtrain, ytrain,
    tuneLength=tuneLength,
    trControl=fullCtrl,
    metric='Accuracy',
    methodList=c('kernelpls', 'bagEarthGCV', 'knn', 'svmRadial')
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
