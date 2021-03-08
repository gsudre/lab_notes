args = commandArgs(trailingOnly=TRUE)
if (length(args) != 4) {
  stop("Arguments: seed[-1] mymod target root_fname", call.=FALSE)
} else {
  if (args[1] == '-1') {
    myseed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  } else {
    myseed = as.integer(args[1])
  }
  target = args[3]
  mymod = args[2]
  root_fname = args[4]
}

###########

fname = sprintf('%s_%04d.log', root_fname, myseed)
sink(fname, append=FALSE, split=TRUE)
source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')

tuneLength=10
struct_data = read.csv('~/data/baseline_prediction/stripped/structural.csv')
load('~/data/baseline_prediction/struct_thickness.RData')
vdata = cbind(struct_data$Mask.ID...Scan, lh_thickness, rh_thickness)
rm(lh_thickness)
rm(rh_thickness)
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf_base = gf[gf$BASELINE=='BASELINE', ]
my_ids = intersect(gf_base$MRN, struct_data$MRN)
mstruct = mergeOnClosestDate(gf_base, struct_data, my_ids)
rm_me = abs(mstruct$dateX.minus.dateY.months) > 12
mstruct = mstruct[!rm_me, ]
struct_base_vdata = merge(mstruct$Mask.ID...Scan, vdata, by.x=1, by.y=1, all.y=F, all.x=T)
rm(vdata)

X = struct_base_vdata[, 2:ncol(struct_base_vdata)]
rm_me = colSums(is.na(X)) > 0
X = X[, !rm_me]

keep_me = mstruct$age <= 12
X = X[keep_me, ]

if (target == 'inatt') {
    y = mstruct$inatt3_named
    y = factor(y, levels=c('low', 'medium', 'high'))
} else {
    y = mstruct$HI3_named
    y = factor(y)
}
y = y[keep_me]

set.seed(myseed)
split <- createDataPartition(y, p = .8, list = FALSE)
Xtrain <- X[ split, ]
ytrain <- y[ split ]
Xtest  <- X[-split, ]
ytest = y[-split]

set.seed(myseed)
index <- createMultiFolds(ytrain, k = 5, times = 5)
fullCtrl <- trainControl(method = "repeatedcv",
                         index = index,
                         savePredictions="final")

library(doMC)
ncpus <- detectBatchCPUs()
registerDoMC(ncpus)

library(caretEnsemble)
set.seed(myseed)
model_list <- caretList(
  Xtrain, ytrain,
  tuneLength=tuneLength,
  trControl=fullCtrl,
  methodList=c(mymod)
  )

model_perf = data.frame(lapply(model_list, function(d) getTrainPerf(d)[1]))
names(model_perf) = names(model_list)
model_preds <- lapply(model_list, predict, newdata=Xtest)
model_preds <- lapply(model_preds, function(d, obs) postResample(d, obs)[1], obs=ytest)
model_preds <- data.frame(model_preds)
names(model_preds) = names(model_list)
print('training')
print(model_perf)
print('testing')
print(model_preds)
print(sprintf('No information rate: Accuracy=%f', max(table(ytrain)/length(ytrain))))

sink()
