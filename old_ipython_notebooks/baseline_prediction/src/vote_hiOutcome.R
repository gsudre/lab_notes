args = commandArgs(trailingOnly=TRUE)
s=as.numeric(args[1])
ntimes = 50
model = args[2]
myseed = 1234
tuneLength = 10
cpuDiff = 0
dsets = c('prs', 'geospatial', 'neuropsych', 'struct_rois', 'dti_tracts')

source('~/ncr_notebooks/baseline_prediction/src/load_voting_data.R')
out_fname = sprintf('/data/NCR_SBRB/loocv/hiOutcome_%s/s%03d.log', model, s)
sink(out_fname, append=FALSE, split=TRUE)

vote_hiOutcome = function(X, s, uni=T, pca=T, do_rfe=F) {
  # need to use it from merged instead of gf_base because unique() use dfor my_ids changes the order of MRNs!
  y = merged$HI3_named
  y = factor(y)

  Xtrain <- X[ -s, ]
  ytrain <- y[ -s ]
  Xtest  <- X[s, ]

  # remove training examples that are all NaN. Note that this will never happen for testing!
  rm_me = rowSums(is.na(Xtrain)) == ncol(Xtrain)
  Xtrain = Xtrain[!rm_me,]
  ytrain = ytrain[!rm_me]

  if (pca) {
    pp <- preProcess(Xtrain, method = c('pca'), thresh=.9)
    Xtrain<- predict(pp, Xtrain)
    Xtest <- predict(pp, Xtest)
  }

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
                method=model,
                preProcess = c('medianImpute'))
  return(predict(mymod, newdata=Xtest))
}

nsubjs = nrow(gf_base)
library(doMC)
ncpus <- detectBatchCPUs()
njobs <- ncpus - cpuDiff
registerDoMC(njobs)

# get a vote for each dset if the participant has data in the domain
print(sprintf('LOOCV %d / %d', s, nsubjs))
ptm <- proc.time()
preds = c()
for (dset in dsets) {
  eval(parse(text=sprintf('nfeats = ncol(%s)', dset)))
  eval(parse(text=sprintf('na_feats = sum(is.na(%s[s, ]))', dset)))
  if (na_feats < nfeats) {
    eval(parse(text=sprintf('X = %s', dset)))
    print(sprintf('evaluating %s', dset))
    res = vote_hiOutcome(X, s, pca=F)[[1]]
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
