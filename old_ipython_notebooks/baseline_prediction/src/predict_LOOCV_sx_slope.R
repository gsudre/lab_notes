source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')

args = commandArgs(trailingOnly=TRUE)
subj=as.numeric(args[1])
model = args[2]
mysx = args[3]
adhdOnly = as.logical(args[4])
out_fname = sprintf('/data/NCR_SBRB/loocv/%s_sx_slope_%s_%s/s%03d.log', mysx, model,
                    adhdOnly, subj)
library(caret)
library(doMC)
ncpus <- detectBatchCPUs()
registerDoMC(ncpus)

# subj=5
# model = 'rf'
# mysx = 'hi'
# adhdOnly = F
# out_fname = sprintf('~/tmp/s%03d.log', subj)

sink(out_fname, append=FALSE, split=TRUE)
ntimes = 50
myseed = 1234
tuneLength = 10
cpuDiff = 0
dsets = c('geospatial', 'prs', 'neuropsych', 'struct_rois', 'dti_tracts')

get_SX_slope = function (df, ids) {
  inatt = c()
  hi = c()
  mrns = c()
  for (s in ids) {
    idx = df$MRN==s
    # proceed if we have more than one observation in the data
    if (sum(idx) >= 2) {
      inatt = c(inatt, lm(SX_inatt ~ age, data=df[idx,])$coefficients[[2]])
      hi = c(hi, lm(SX_HI ~ age, data=df[idx,])$coefficients[[2]])
      mrns = c(mrns, s)
    }
  }
  res = cbind(inatt, hi)
  colnames(res) = c('inatt', 'hi')
  rownames(res) = mrns
  return(res)
}
source('~/ncr_notebooks/baseline_prediction/src/load_raw_voting_data.R')

nsubjs = nrow(gf_base)

if (adhdOnly) {
  idx = which(gf_base$DX_BASELINE!='NV')
  ids = gf_base[idx,]$MRN
} else {
  idx = 1:nrow(gf_base)
  ids = gf_base$MRN
}
slopes = get_SX_slope(gf, ids)
y = slopes[, mysx]

predict_LOOCV_sx_slope = function(X, s) {
  # recoding SEX
  dummies = dummyVars(~SEX, data=X)
  X = cbind(X, predict(dummies, newdata=X))
  X$SEX = NULL
  
  Xtrain <- X[-s, ]
  ytrain <- y[-s]
  Xtest  <- X[s, ]
  ytest  <- y[s]
  
  # remove training examples that are all NaN. Note that this will never happen for testing!
  rm_me = rowSums(is.na(Xtrain)) == ncol(Xtrain)
  Xtrain = Xtrain[!rm_me,]
  ytrain = ytrain[!rm_me]
  
  pp = preProcess(Xtrain, method=c('medianImpute'))
  Xtrain = predict(pp, Xtrain)
  Xtest = predict(pp, Xtest)
  
  set.seed(myseed)
  fullCtrl <- trainControl(method = "boot632",
                           number = ntimes,
                           savePredictions="final",
                           summaryFunction=defaultSummary)
  set.seed(myseed)
  mymod = train(Xtrain, ytrain,
                tuneLength=tuneLength,
                trControl=fullCtrl,
                method=model)
  
  pred = predict(mymod, newdata=Xtest)
}

print(sprintf('LOOCV %d / %d', subj, nsubjs))
preds = c()
for (dset in dsets) {
  eval(parse(text=sprintf('nfeats = ncol(%s) - 2', dset)))  # remove sex and age
  eval(parse(text=sprintf('na_feats = sum(is.na(%s[subj, ]))', dset)))
  # if this subjects doesn't have all features NA
  if (na_feats < nfeats) {
    eval(parse(text=sprintf('X = %s', dset)))
    print(sprintf('evaluating %s', dset))
    res = predict_LOOCV_sx_slope(X[idx, ], subj)
  }
  else {
    res = NA
  }
  preds = c(preds, res)
}
print(sprintf('evaluating %s', 'all data'))
X = cbind(geospatial, prs, neuropsych, struct_rois, dti_tracts)[idx, ]
# everything in voting_data is already < 12. Let's keep the first age and SEX only,
# which will come from geospatial and therefore there are no NAs (everybody has geospatial data)
X <- X[, !duplicated(colnames(X))]
res = predict_LOOCV_sx_slope(X, subj)
preds = c(preds, res)

names(preds) = c(dsets, 'all')

means = colMeans(slopes)
nir = RMSE(means[mysx], slopes[, mysx])
print(sprintf('NIR %s: %.2f', mysx, nir))

dsets
model
myseed
tuneLength
ntimes
adhdOnly
preds
sink()