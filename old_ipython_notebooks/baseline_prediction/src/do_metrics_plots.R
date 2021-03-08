compile_metrics = function(res, col, hdr) {
  b = vector()
  for (d in res) {
    b = cbind(b, d[,col])
  }
  colnames(b) = hdr
  return(b)
}

par(mfrow=c(2, 3))
boxplot(compile_metrics(mylist, 1, run_models), ylab='logLoss', las=2)
boxplot(compile_metrics(mylist, 2, run_models), ylab='AUC', las=2)
boxplot(compile_metrics(mylist, 3, run_models), ylab='Accuracy', las=2, ylim=c(.4, 1))
library(Hmisc)
res = binconf(x=max(table(groups)), n=length(groups), alpha=.05)
abline(h=res[1], col='red')
abline(h=res[2], col='red', lty=2)
abline(h=res[3], col='red', lty=2)
boxplot(compile_metrics(mylist, 4, run_models), ylab='Kappa', las=2)
boxplot(compile_metrics(mylist, 5, run_models), ylab='Sensitivity', las=2)
boxplot(compile_metrics(mylist, 6, run_models), ylab='Specificity', las=2)