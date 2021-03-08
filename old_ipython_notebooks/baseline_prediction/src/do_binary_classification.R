# save X and y if not already done so
fname = sprintf('%s_Xy.RData', root_fname)
if(!file.exists(fname)){
  save(X, y, file=fname, compress=T)
}

set.seed(myseed)
split <- createDataPartition(y, p = .8, list = FALSE)
Xtrain <- X[ split, ]
ytrain <- y[ split ]
Xtest  <- X[-split, ]
ytest = y[-split]

pp = preProcess(Xtrain, method=c('YeoJohnson', 'center', 'scale'))
filtXtrain = predict(pp, Xtrain)
filtXtest = predict(pp, Xtest)
nzv = nearZeroVar(filtXtrain)
print(nzv)
if (length(nzv) > 0) {
  filtXtrain = filtXtrain[, -nzv]
}
correlations = cor(filtXtrain, use='na.or.complete')

highCorr = findCorrelation(correlations, cutoff=.75)
print(length(highCorr))
if (length(highCorr) > 0) {
  noncorrXtrain = filtXtrain[, -highCorr]
  noncorrXtest = filtXtest[, -highCorr]
} else {
  noncorrXtrain = filtXtrain
  noncorrXtest = filtXtest
}

library(pROC)

set.seed(myseed)
index <- createMultiFolds(ytrain, k = 5, times = 5)

library(doMC)
ncpus <- detectBatchCPUs()
njobs <- ncpus - cpuDiff
registerDoMC(njobs)
selFunc = 'best'  # maybe try oneSE and tolerance as well?

set.seed(myseed)
fullCtrl <- trainControl(method = "repeatedcv",
                         index = index,
                         search='grid',
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE)

ptm <- proc.time()
m1 <- train(noncorrXtrain, ytrain,
            method = mymod,
            trControl = fullCtrl,
            tuneLength = tuneLength,
            metric = 'ROC')
print(proc.time() - ptm)
print(m1)
pred = predict(m1, noncorrXtest)
print(postResample(pred, ytest))
# note that this value below is an approximation. When we run the classSummary
# function it takes the probabilities of each class into consideration, which
# can be more reliable?
print(roc(as.numeric(ytest), as.numeric(pred)))

fname = sprintf('%s_%04d.RData', root_fname, myseed)
save_list = c('m1', 'myseed', 'index', 'split')
save(list=save_list, file=fname)
