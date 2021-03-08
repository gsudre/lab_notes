args = commandArgs(trailingOnly=TRUE)
subj=as.numeric(args[1])
ntimes = 50
model = args[2]
myseed = 1234
tuneLength = 10
cpuDiff = 0
dsets = c('geospatial', 'prs', 'neuropsych', 'struct_rois', 'dti_tracts')

source('~/ncr_notebooks/baseline_prediction/src/load_raw_voting_data.R')
# out_fname = sprintf('/data/NCR_SBRB/loocv/hiOutcome2Class_%s/s%03d.log', model, subj)
out_fname = sprintf('~/tmp/s%03d.log', subj)
sink(out_fname, append=FALSE, split=TRUE)

vote_hiOutcome2Class = function(X, y, s) {
  dummies = dummyVars(~SEX, data=X)
  X = cbind(X, predict(dummies, newdata=X))
  X$SEX = NULL
  
  Xtrain <- X[ -s, ]
  ytrain <- y[ -s ]
  Xtest  <- X[s, ]

  # remove training examples that are all NaN. Note that this will never happen for testing!
  rm_me = rowSums(is.na(Xtrain)) == ncol(Xtrain)
  Xtrain = Xtrain[!rm_me,]
  ytrain = ytrain[!rm_me]

  pp = preProcess(Xtrain, method=c('medianImpute'))
  Xtrain = predict(pp, Xtrain)
  Xtest = predict(pp, Xtest)
  
  set.seed(myseed)
  fullCtrl <- trainControl(method = "boot",
                           number = ntimes,
                           savePredictions="final",
                           classProbs=TRUE,
                           summaryFunction=defaultSummary,
                           allowParallel=T)
  set.seed(myseed)
  mymod = train(Xtrain, ytrain,
                tuneLength=tuneLength,
                trControl=fullCtrl,
                metric='Accuracy',
                method=model)
  print(varImp(mymod))
  return(predict(mymod, newdata=Xtest, type='prob'))
}

nsubjs = nrow(gf_base)
library(doMC)
ncpus <- detectBatchCPUs()
njobs <- ncpus - cpuDiff
registerDoMC(njobs)

y = gf_base$HI3_named
keep_me = y!='never_affected'
y = factor(y[keep_me])

# get a vote for each dset if the participant has data in the domain
print(sprintf('LOOCV %d / %d', subj, nsubjs))
preds = c()
for (dset in dsets) {
  eval(parse(text=sprintf('nfeats = ncol(%s)', dset)))
  eval(parse(text=sprintf('na_feats = sum(is.na(%s[subj, ]))', dset)))
  if (na_feats < nfeats) {
    eval(parse(text=sprintf('X = %s', dset)))
    print(sprintf('evaluating %s', dset))
    res = vote_hiOutcome2Class(X[keep_me,], y, subj)['rapid_improvers'][[1]]
  } else {
    res = NA
  }
  preds = c(preds, res)
}
dsets
model
myseed
tuneLength
ntimes
preds
sink()
