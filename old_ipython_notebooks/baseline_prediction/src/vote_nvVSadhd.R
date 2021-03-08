args = commandArgs(trailingOnly=TRUE)
s=as.numeric(args[1])
model = args[2]

ntimes = 50
myseed = 1234
tuneLength = 10
cpuDiff = 0
if (model %in% c('ranger')) {
    allowParallel = FALSE
} else {
    allowParallel = TRUE
}
out_fname = sprintf('/data/NCR_SBRB/loocv/nvVSadhd_%s/s%03d.log', model, s)
sink(out_fname, append=FALSE, split=TRUE)
source('~/ncr_notebooks/baseline_prediction/src/load_voting_data.R')
dsets = c('prs', 'geospatial', 'neuropsych', 'struct_rois', 'dti_tracts')#,
          # 'brain_fa', 'brain_ad', 'brain_rd')#,
        #   'brain_thickness', 'brain_volume', 'brain_area')
vote_nvVSadhd = function(X, s, uni=T, pca=T, do_rfe=F) {
  # need to use it from merged instead of gf_base because unique() use dfor my_ids changes the order of MRNs!
  y = merged$DX_BASELINE
  y[y!='NV'] = 'ADHD'
  y = factor(y, levels=c('ADHD', 'NV'))

  Xtrain <- X[ -s, ]
  ytrain <- y[ -s ]
  Xtest  <- X[s, ]

  # remove training examples that are all NaN. Note that this will never happen for testing!
  rm_me = rowSums(is.na(Xtrain)) == ncol(Xtrain)
  Xtrain = Xtrain[!rm_me,]
  ytrain = ytrain[!rm_me]

  # keep only good univariate variables
  if (uni > 0) {
      cl <- makeCluster(ncpus)
      pvals = parSapply(cl, Xtrain, function(d, ytrain) t.test(d ~ ytrain)$p.value, ytrain)
      stopCluster(cl)
    good_vars = which(pvals <= uni)
    if (length(good_vars) > 1) {
      Xtrain = Xtrain[, good_vars]
      keep_me = sapply(colnames(Xtrain), function(d) which(colnames(Xtest) == d))
      Xtest = Xtest[, keep_me]
    } else {
      return(data.frame(ADHD = NA, NV = NA))
    }
  }

  if (pca) {
    pp <- preProcess(Xtrain, method = c('pca'), thresh=.9)
    Xtrain<- predict(pp, Xtrain)
    Xtest <- predict(pp, Xtest)
  }

  if (do_rfe) {
    set.seed(myseed)

    ctrl <- rfeControl(functions = treebagFuncs,
                       method = "boot",
                       repeats = 5,
                       verbose = FALSE
                       )
    subsets = floor(seq(from=1, to=ncol(Xtrain)/2, length.out=min(10, ncol(Xtrain)/2)))
    mymod <- rfe(Xtrain, ytrain,
                     sizes = subsets,
                     rfeControl = ctrl)
  } else{
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
  }
  return(predict(mymod, newdata=Xtest, type='prob'))
}

nsubjs = nrow(gf_base)
library(doMC)
ncpus <- detectBatchCPUs()
njobs <- ncpus - cpuDiff
registerDoMC(njobs)

ptm <- proc.time()
# get a vote for each dset if the participant has data in the domain
print(sprintf('LOOCV %d / %d', s, nsubjs))
preds = c()
for (dset in dsets) {
    eval(parse(text=sprintf('nfeats = ncol(%s)', dset)))
    eval(parse(text=sprintf('na_feats = sum(is.na(%s[s, ]))', dset)))
    if (na_feats < nfeats) {
      eval(parse(text=sprintf('X = %s', dset)))
      if (dset %in% c('brain_fa', 'brain_ad', 'brain_rd',
                      'brain_thickness', 'brain_volume', 'brain_area')) {
        print('Univariate filters + PCA')
        uni = 0.05
        pca = T
        do_rfe = F
      } else {
        uni = 0
        pca = F
        do_rfe = F
      }
      print(sprintf('evaluating %s', dset))
      uni
      pca
      do_rfe
      res = vote_nvVSadhd(X, s, uni=uni, pca=pca, do_rfe=do_rfe)['ADHD'][[1]]
  }
  else {
      res = NA
  }
  preds = c(preds, res)
}
print(proc.time() - ptm)
dsets
model
myseed
tuneLength
ntimes
preds
sink()
