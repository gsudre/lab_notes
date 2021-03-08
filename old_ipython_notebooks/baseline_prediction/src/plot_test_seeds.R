get_random_vote = function(ts, nperms) {
  res = matrix(nrow = nperms, ncol = 4)
  colnames(res) = c('accuracy', 'auc', 'sensitivity', 'specificity')
  for (i in 1:nperms) {
    perm_ts = ts
    perm_ts$obs = sample(ts$obs)
    if (length(levels(ts$obs)) > 2) {
      tmp = multiClassSummary(perm_ts, lev=levels(ts$obs))
      res[i, ] = tmp[c('Accuracy', 'Mean_ROC', 'Mean_Sensitivity', 'Mean_Specificity')]
    } else {
      tmp = multiClassSummary(perm_ts, lev=levels(ts$obs))
      res[i, ] = tmp[c('Accuracy', 'ROC', 'Sensitivity', 'Specificity')]
    }
  }
  res = as.data.frame(res)
  return(res)
}

library(caret)
models = c('AdaBag', 'AdaBoost')
fname_path = '~/data/baseline_prediction/results/bw/'
fname_root = 'allNeuropsychCleanGASsplit'
pct = .95
nperms = 1000

model_results = list()
rnd_res = c()
for (m in models) {
  load(sprintf('%s/%s_%s_Xy.RData', fname_path, fname_root, m))
  tfiles = list.files(path=fname_path, pattern=glob2rx(sprintf('%s_%s_???*.RData', fname_root, m)))
  nfiles = length(tfiles)
  print(sprintf('Found %d test result files for %s_%s.', nfiles, fname_root, m))
  res = matrix(nrow = nfiles, ncol = 4)
  colnames(res) = c('accuracy', 'auc', 'sensitivity', 'specificity')
  for (f in 1:nfiles) {
    cat(sprintf('%d ', f))
    load(sprintf('%s/%s', fname_path, tfiles[f]))

    Xtrain <- X[ split, ]
    ytrain <- y[ split ]
    Xtest  <- X[-split, ]
    ytest = y[-split]

    pp = preProcess(Xtrain, method=c('YeoJohnson', 'center', 'scale'))
    filtXtrain = predict(pp, Xtrain)
    filtXtest = predict(pp, Xtest)
    nzv = nearZeroVar(filtXtrain)
    if (length(nzv) > 0) {
      filtXtrain = filtXtrain[, -nzv]
    }
    correlations = cor(filtXtrain, use='na.or.complete')

    highCorr = findCorrelation(correlations, cutoff=.75)
    if (length(highCorr) > 0) {
      noncorrXtrain = filtXtrain[, -highCorr]
      noncorrXtest = filtXtest[, -highCorr]
    } else {
      noncorrXtrain = filtXtrain
      noncorrXtest = filtXtest
    }

    preds = predict(m1, newdata=noncorrXtest)
    probs = predict(m1, newdata=noncorrXtest, type="prob")
    ts = data.frame(obs=ytest, pred=preds, probs)
    if (length(levels(ts$obs)) > 2) {
      tmp = multiClassSummary(ts, lev=levels(ts$obs))
      res[f, ] = tmp[c('Accuracy', 'Mean_ROC', 'Mean_Sensitivity', 'Mean_Specificity')]
    } else {
      # we use multi here because it gives Accuracy, while two doesn't.
      # but the titles change
      tmp = multiClassSummary(ts, lev=levels(ts$obs))
      res[f, ] = tmp[c('Accuracy', 'ROC', 'Sensitivity', 'Specificity')]
    }

    rnd_res = rbind(rnd_res, get_random_vote(ts, ceiling(nperms / nfiles)))
  }
  cat('\n')
  eval(parse(text=sprintf('model_results$%s = res', m)))
}

qtiles = vector(mode = 'numeric')
for (k in 1:ncol(rnd_res)) {
  qtiles = c(qtiles, quantile(rnd_res[, k], pct, na.rm=T))
}
names(qtiles) = colnames(rnd_res)

par(oma=c(0,0,2,0))
par(mfrow=c(2, 2))
boxplot(sapply(model_results, function(d) d[,c('accuracy')]), ylab='accuracy', las=2, ylim=c(.4, 1))
abline(h=qtiles['accuracy'], col='red', lty=2)
boxplot(sapply(model_results, function(d) d[,c('auc')]), ylab='AUC', las=2, ylim=c(.4, 1))
abline(h=qtiles['auc'], col='red', lty=2)
boxplot(sapply(model_results, function(d) d[,c('sensitivity')]), ylab='sensitivity', las=2, ylim=c(.4, 1))
abline(h=qtiles['sensitivity'], col='red', lty=2)
boxplot(sapply(model_results, function(d) d[,c('specificity')]), ylab='specificity', las=2, ylim=c(.4, 1))
abline(h=qtiles['specificity'], col='red', lty=2)
title(main=sprintf('%s_nperms%d_pct%.2f', fname_root, nperms, pct),outer=T)
