do_classification <- function(run_models, inTrain, njobs) {
  require(doParallel)
  require(caret)
  # setting up parallelization
  # registerDoParallel(njobs,cores=njobs)
  # getDoParWorkers()
  
  require(doMC)
  print(njobs)
  registerDoMC(cores=njobs)
  
  # spitting out informative stuff
  gf_fname
  data_fname
  group_var
  metric
  default_preproc
  
  rndForestAll = c()
  lrAll = c()
  lsvmAll = c()
  rsvmAll = c()
  xgbAll = c()
  gbmAll = c()
  save_list = c('train_idx')
  for (m in run_models) {
    save_list = c(save_list, sprintf('%sFit', m))
  }
  
  # for each train/test split
  for (i in 1:length(inTrain)) {
    cat(sprintf('\nWorking on split %d of %d', i, length(inTrain)))
    # Xtrain = ldata[ inTrain[[i]], ]
    Xtrain = as.data.frame(ldata[ inTrain[[i]], ])
    colnames(Xtrain) = colnames(ldata) # do this in two steps because it wasn't working to define col.names in the constructor
    ytrain = groups[ inTrain[[i]]]
    # Xtest  = ldata[-inTrain[[i]], ]
    Xtest = as.data.frame(ldata[ -inTrain[[i]], ])
    colnames(Xtest) = colnames(ldata)  # do this in two steps because it wasn't working to define col.names in the constructor
    ytest = groups[-inTrain[[i]]]
    
    if ('rndForest' %in% run_models) {
      ptm <- proc.time()
      cat('\nRunning random forests\n')
      rndForestGrid <- expand.grid(.mtry=c(1:sqrt(ncol(Xtrain))))
      rndForestFit <- train(Xtrain, ytrain,
                            method = "rf",
                            trControl = ctrl_cv,
                            tuneGrid=rndForestGrid,
                            metric = metric,
                            preProcess = default_preproc)
      rndForestRes = eval_model(rndForestFit, Xtest, ytest)
      rndForestAll = rbind(rndForestAll, rndForestRes)
      print(proc.time() - ptm)
    }
    
    if ('lr' %in% run_models) {
      ptm <- proc.time()
      cat('\nRunning logistic regression\n')
      lrGrid <- expand.grid(.nIter=seq(1,50,5))
      lrFit <- train(Xtrain, ytrain,
                     method = "LogitBoost",
                     trControl = ctrl_cv,
                     tuneGrid=lrGrid,
                     metric = metric,
                     preProc = default_preproc)
      lrRes = eval_model(lrFit, Xtest, ytest)
      lrAll = rbind(lrAll, lrRes)
      print(proc.time() - ptm)
    }
    
    if ('lsvm' %in% run_models) {
      ptm <- proc.time()
      cat('\nRunning Linear SVM\n')
      lsvmGrid <- expand.grid(.C=c(.01, .1, 1, 10, 100, 10^3, 10^4))
      lsvmFit <- train(Xtrain, ytrain,
                       method = "svmLinear",
                       trControl = ctrl_cv,
                       tuneGrid=lsvmGrid,
                       metric = metric,
                       preProc = default_preproc)
      lsvmRes = eval_model(lsvmFit, Xtest, ytest)
      lsvmAll = rbind(lsvmAll, lsvmRes)
      print(proc.time() - ptm)
    }
    
    if ('rsvm' %in% run_models) {
      ptm <- proc.time()
      cat('\nRunning RBF SVM\n')
      rsvmGrid <- expand.grid(.C=seq(1e-2, 1e+4, length.out=7),
                              .sigma=seq(1e-5, 1e+1, length.out=7))
      rsvmFit <- train(Xtrain, ytrain,
                       method = "svmRadial",
                       verbose=F,
                       trControl = ctrl_cv,
                       tuneGrid=rsvmGrid,
                       metric = metric,
                       preProc = default_preproc)
      rsvmRes = eval_model(rsvmFit, Xtest, ytest)
      rsvmAll = rbind(rsvmAll, rsvmRes)
      print(proc.time() - ptm)
    }
    
    if ('xgb' %in% run_models) {
      ptm <- proc.time()
      cat('\nRunning XGBoost\n')
      xgbGrid <- expand.grid(.nrounds = 1000,
                             .eta = c(.0001, .001, 0.01, 0.05, 0.1),
                             .max_depth = c(2,4,6,8,10,14),
                             .gamma = 1,
                             .colsample_bytree=.8,
                             .min_child_weight=1)
      xgbFit <- train(Xtrain, ytrain,
                      method = "xgbTree",
                      trControl = ctrl_cv,
                      # we need this for the cluster, where the C++ backend already does multicore
                      nthread = 1,
                      tuneGrid=xgbGrid,
                      metric = metric,
                      preProc = default_preproc)
      xgbRes = eval_model(xgbFit, Xtest, ytest)
      xgbAll = rbind(xgbAll, xgbRes)
      print(proc.time() - ptm)
    }
    
    if ('gbm' %in% run_models) {
      ptm <- proc.time()
      cat('\nRunning GBM\n')
      gbmGrid <-  expand.grid(.interaction.depth = c(1:sqrt(ncol(Xtrain))),
                              .n.trees = seq(1,501,10),
                              .shrinkage = seq(.0005, .05,.0005),
                              .n.minobsinnode = 10) #c(5, 10, 15, 20))
      gbmFit <- train(Xtrain, ytrain,
                      method = "gbm",
                      verbose=F,
                      trControl = ctrl_cv,
                      tuneGrid=gbmGrid,
                      metric = metric,
                      preProc = default_preproc)
      gbmRes = eval_model(gbmFit, Xtest, ytest)
      gbmAll = rbind(gbmAll, gbmRes)
      print(proc.time() - ptm)
    }
    
    # saving fit models
    fname = sprintf('%s_split%02d.RData', root_fname, i)
    train_idx = inTrain[[i]]
    save(list=save_list, file=fname)
  }
  stopImplicitCluster()
}