args = commandArgs(trailingOnly=TRUE)
s=as.numeric(args[1])
ntimes = 50
model = args[2]
myseed = 1234
tuneLength = 10
cpuDiff = 0
dsets = c('prs', 'geospatial', 'neuropsych', 'struct_rois', 'dti_tracts')

source('~/ncr_notebooks/baseline_prediction/src/load_voting_data.R')
out_fname = sprintf('/data/NCR_SBRB/loocv/hiOutcomeStepwise_%s/s%03d.log', model, s)
sink(out_fname, append=FALSE, split=TRUE)

if (model %in% c('ranger')) {
  allowParallel = FALSE
} else {
  allowParallel = TRUE
}

vote_hiOutcomeStepwise = function(X, s, uni=T, pca=T, do_rfe=F) {
  y = as.character(merged$HI3_named)
  y[y!='never_affected'] = 'affected'
  y = factor(y, levels=c('never_affected', 'affected'))
  
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
                           summaryFunction=twoClassSummary,
                           allowParallel=allowParallel)
  set.seed(myseed)
  mymod = train(Xtrain, ytrain,
                tuneLength=tuneLength,
                trControl=fullCtrl,
                metric='ROC',
                method=model,
                preProcess = c('medianImpute'))
  
  pred = predict(mymod, newdata=Xtest)
  print(pred)
  if (pred == 'never_affected') {
    return(as.character(pred))
  } else {
    y = as.character(merged$HI3_named)
    Xtrain <- X[ -s, ]
    ytrain <- y[ -s ]
    
    keep_me = ytrain!='never_affected'
    ytrain = ytrain[keep_me]
    ytrain = factor(ytrain, levels=c('rapid_improvers', 'severe'))
    Xtrain = Xtrain[keep_me,]
    # Xtest is the same as before
    
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
                             summaryFunction=twoClassSummary,
                             allowParallel=allowParallel)
    set.seed(myseed)
    mymod = train(Xtrain, ytrain,
                  tuneLength=tuneLength,
                  trControl=fullCtrl,
                  metric='ROC',
                  method=model,
                  preProcess = c('medianImpute'))
    
    pred = predict(mymod, newdata=Xtest)
    return(as.character(pred))
  }
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
  res = vote_hiOutcomeStepwise(X, s, pca=F)[[1]]
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
