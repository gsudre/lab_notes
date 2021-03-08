library(caret)

compile_metrics = function(res, col, hdr) {
  b = vector()
  for (d in res) {
    b = cbind(b, d[,col])
  }
  colnames(b) = hdr
  return(b)
}

do_metrics_plots = function(root_fname) {
  load(sprintf('%s.RData', root_fname))
  res = collect_results(root_fname)
  mylist = res[[1]]
  run_models = res[[2]]
  
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
}

merge_func = function(id, df1, df2, my_ids, x.id, x.date, y.id, y.date) {
  eval(parse(text=sprintf('d1 = subset(df1, %s==id)', x.id)))
  eval(parse(text=sprintf('d2 = subset(df2, %s==id)', y.id)))
  # if the id is not in df2
  if (nrow(d2) == 0) {
    d3 = merge(d1,d2,by.x=x.id, by.y=y.id, all.x=T)
  } else {
    d1$indices <- sapply(d1[, eval(x.date)], function(d) which.min(abs(d2[, eval(y.date)] - d)))
    d2$indices <- 1:nrow(d2)
    d3 = merge(d1,d2,by.x=c(x.id,'indices'), by.y=c(y.id, 'indices'))
    d3$indices = NULL
  }
  if (x.date == y.date) {
    x.date = sprintf('%s.x', x.date)
    y.date = sprintf('%s.y', y.date)
  }
  d3$dateX.minus.dateY.months = as.numeric(d3[, eval(x.date)] - d3[, eval(y.date)]) / 30
  return(d3)
}

mergeOnClosestDate = function(df1, df2, my_ids, x.id='MRN', x.date='DOA', y.id='MRN', y.date='DOA') {
  # replace %Y format by %y to make reading it in easier (if needed)
  df1[, eval(x.date)] = gsub("[0-9]{2}([0-9]{2})$", "\\1", df1[, eval(x.date)]) 
  df2[, eval(y.date)] = gsub("[0-9]{2}([0-9]{2})$", "\\1", df2[, eval(y.date)])
  # convert from factors to dates
  df1[, eval(x.date)] = as.Date(df1[, eval(x.date)], format='%m/%d/%y')
  df2[, eval(y.date)] = as.Date(df2[, eval(y.date)], format='%m/%d/%y')
  
  # apply the merge function to all ids
  merged_ids <- lapply(my_ids, merge_func, df1, df2, my_ids, x.id, x.date, y.id, y.date)
  # bind by rows the ID-specific merged dataframes
  res_df = do.call(rbind, merged_ids)
  return(res_df)
}

# figure out the baseline scans. For 9 months, min_time_diff = .75
get_baseline_scans = function(data, min_time_diff=0) {
    keep = vector()
    for (mrn in unique(data$MRN)) {
        idx = which(data$MRN == mrn)
        mrn_ages = data$age_at_scan[idx]
        if (length(mrn_ages) > 1 && diff(mrn_ages)[1] < min_time_diff) {
            cat(sprintf('ERROR: Scans for %d are not more than %.2f apart!',
                        mrn, min_time_diff))
        }
        keep = append(keep, idx[sort(mrn_ages, index.return=T)$ix[1]])
    }
    data = data[keep,]
    return (data)
}

# returns metrics for different models
eval_model = function(fit, data, labels, lev=levels(labels)) {
  preds = predict(fit, newdata=data)
  probs = predict(fit, newdata=data, type="prob")
  ts = data.frame(obs=labels, pred=preds, probs)
  if (length(lev) == 2) {
    res = twoClassSummary(ts, lev=lev)
  } else {
    res = multiClassSummary(ts, lev=lev)
  }
  return(res)
}

detectBatchCPUs <- function() { 
  ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK")) 
  if (is.na(ncores)) { 
    ncores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE")) 
  } 
  if (is.na(ncores)) { 
    return(8) # for laptop (use 4 for helix)
  } 
  return(ncores) 
}

do_classification <- function(root_fname, ldata, groups, run_models, inTrain, ctrl_cv, njobs, default_preproc, metric) {
  require(doParallel)
  require(caret)
  
  require(doMC)
  print(njobs)
  registerDoMC(cores=njobs)
  
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

runInCluster <- function(root_fname, ncv=5, nrepeatcv=2, nsplits=5, train_test_ratio=.85, cpuDiff=0, fixed_seed=T,
                         run_models=c('rndForest', 'lr', 'rsvm', 'xgb'), metric = "ROC",
                         default_preproc = c("center", 'scale'), default_options = c()) {
  # loads ldata and groups
  load(sprintf('%s.RData', root_fname))
  sink(sprintf('%s.log', root_fname), append=FALSE, split=TRUE)
  ncpus <- detectBatchCPUs()
  njobs <- ncpus - cpuDiff
  print(njobs)
  if (fixed_seed) {
    set.seed(107)
  }
  if (length(levels(groups)) > 2) {
    my_func = multiClassSummary
    if (metric == 'ROC') {
      metric = 'Mean_ROC'
    }
  } else {
    my_func = twoClassSummary
  }
  inTrain = createDataPartition(y = groups,
                                times= nsplits,
                                p=train_test_ratio)
  ctrl_cv <- trainControl(method = "repeatedcv",
                          number = ncv,
                          repeats = nrepeatcv,
                          classProbs = TRUE,
                          returnData = FALSE,
                          summaryFunction = my_func,
                          preProcOptions = default_options,
                          search='grid')
  do_classification(root_fname, ldata, groups, run_models, inTrain, ctrl_cv, njobs, default_preproc, metric)
  sink()
}

collect_results <- function(root_fname) {
  load(sprintf('%s.RData', root_fname))
  # open first results file to figure out model names
  var_names = load(sprintf('%s_split%02d.RData', root_fname, 1))
  
  # figure out model names and create appropriate collection files
  run_models = c()
  for (v in var_names) {
    if (grepl('Fit$', v)) {
      run_models = c(run_models, sub('Fit', '', v))
    }
  }
  
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
  
  root_dir = sprintf('%s/', dirname(root_fname))
  files = list.files(path=root_dir, pattern=sprintf('%s_split*', basename(root_fname)))
  # for each results file
  for (f in files) {
    load(sprintf('%s/%s', root_dir, f))
    Xtrain = as.data.frame(ldata[train_idx, ])
    colnames(Xtrain) = colnames(ldata)
    ytrain = groups[train_idx]
    Xtest = as.data.frame(ldata[-train_idx, ])
    colnames(Xtest) = colnames(ldata)
    ytest = groups[-train_idx]
    
    for (m in run_models) {
      eval(parse(text=sprintf('%sRes = eval_model(%sFit, Xtest, ytest)', m, m)))
      eval(parse(text=sprintf('%sAll = rbind(%sAll, %sRes)', m, m, m)))
    }
  }
  mylist = list(rndForestAll, lrAll, lsvmAll, rsvmAll, xgbAll, gbmAll)
  return(list(mylist, run_models))
}

# straight-up copy of tinGraphs, but it doesn't stupidly open new devices for each figure
tinGraphs2 = function (res, DESIGN = NULL, x_axis = NULL, y_axis = NULL, inference.info = NULL, 
               color.by.boots = TRUE, boot.cols = c("plum4", "darkseagreen", 
                                                    "firebrick3"), fi.col = NULL, fi.pch = NULL, fii.col = NULL, 
               fii.pch = NULL, fj.col = NULL, fj.pch = NULL, col.offset = NULL, 
               constraints = NULL, xlab = NULL, ylab = NULL, main = NULL, 
               bootstrapBars = TRUE, correlationPlotter = TRUE, showHulls = 0.95, 
               biplots = FALSE) 
{
  pca.types <- c("tepBADA")
  ca.types <- c("tepDICA")
  if (class(res)[1] == "tinpoOutput") {
    if (length(res) == 2) {
      inference.info <- res$Inference.Data
    }
    res <- res$Fixed.Data
  }
  tepPlotInfo <- NULL
  if (class(res)[1] == "texpoOutput") {
    if (length(res) == 2) {
      tepPlotInfo <- res$Plotting.Data
    }
    res <- res$TExPosition.Data
  }
  if (is.null(inference.info)) {
    stop("inGraphs requires inference.info")
  }
  component.p.order <- order(inference.info$components$p.vals)
  which.axes <- sort(component.p.order[1:2])
  if (is.null(x_axis)) {
    x_axis <- which.axes[1]
  }
  if (is.null(y_axis)) {
    y_axis <- which.axes[2]
  }
  if (!(class(res)[1] %in% c(pca.types, ca.types))) {
    stop("Unknown TExPosition class. Plotting has stopped.")
  }
  else {
    if (is.null(main)) {
      main <- paste("Inferential ", deparse(substitute(res)), 
                    ". Omni p=", inference.info$omni$p.val, sep = "")
    }
    if (length(unlist(strsplit(main, ""))) > 40) {
      main <- paste("Inferential Results. Omni p=", inference.info$omni$p.val, 
                    sep = "")
    }
    if (is.null(xlab)) {
      xlab <- paste("Component ", x_axis, " variance: ", 
                    round(res$t[x_axis], 3), "%, p=", inference.info$components$p.vals[x_axis], 
                    sep = "")
    }
    else {
      xlab <- paste(xlab, "; p=", inference.info$components$p.vals[x_axis], 
                    sep = "")
    }
    if (is.null(ylab)) {
      ylab <- paste("Component ", y_axis, " variance: ", 
                    round(res$t[y_axis], 3), "%, p=", inference.info$components$p.vals[y_axis], 
                    sep = "")
    }
    else {
      ylab <- paste(ylab, "; p=", inference.info$components$p.vals[y_axis], 
                    sep = "")
    }
    if (!is.null(tepPlotInfo)) {
      if (!(nrow(res$fi) == nrow(tepPlotInfo$fi.col)) || 
          !(nrow(res$fj) == nrow(tepPlotInfo$fj.col)) || 
          !(nrow(res$fii) == nrow(tepPlotInfo$fii.col))) {
        print("Dimension mismatch. tepPlotInfo will be reset, no hulls can be shown.")
        tepPlotInfo$fii.col <- NULL
        tepPlotInfo$fi.col <- NULL
        tepPlotInfo$fj.col <- NULL
        tepPlotInfo$constraints <- NULL
      }
    }
    else {
      tepPlotInfo <- list(fii.col = NULL, fi.col = NULL, 
                          fj.col = NULL, constraints = NULL)
    }
    if (is.null(fii.col) || is.null(fi.col) || nrow(fi.col) != 
        nrow(res$fi) || nrow(fii.col) != nrow(res$fi)) {
      if (is.null(tepPlotInfo$fii.col) || is.null(tepPlotInfo$fi.col)) {
        if (is.null(DESIGN)) {
          stop("fii.col and DESIGN are NULL. You must provide one or the other.")
        }
        else {
          DESIGN <- texpoDesignCheck(DATA = NULL, DESIGN = DESIGN, 
                                     make_design_nominal = FALSE)
          obs.cols <- createColorVectorsByDesign(DESIGN, 
                                                 offset = col.offset)
          fii.col <- obs.cols$oc
          fi.col <- obs.cols$gc
        }
      }
      else {
        fii.col <- tepPlotInfo$fii.col
        fi.col <- tepPlotInfo$fi.col
      }
    }
    if (is.null(fj.col) || nrow(fj.col) != nrow(res$fj)) {
      if (is.null(tepPlotInfo$fj.col)) {
        fj.col <- createColorVectorsByDesign(matrix(1, 
                                                    nrow(res$fj), 1), hsv = FALSE, offset = col.offset)$oc
      }
      else {
        fj.col <- tepPlotInfo$fj.col
      }
    }
    if (is.null(fi.pch) || nrow(fi.pch) != nrow(res$fi)) {
      if (is.null(tepPlotInfo$fi.pch)) {
        fi.pch <- as.matrix(rep(21, nrow(res$fi)))
      }
      else {
        fi.pch <- tepPlotInfo$fi.pch
      }
    }
    if (is.null(fii.pch) || nrow(fii.pch) != nrow(res$fii)) {
      if (is.null(tepPlotInfo$fii.pch)) {
        fii.pch <- as.matrix(rep(21, nrow(res$fii)))
      }
      else {
        fii.pch <- tepPlotInfo$fii.pch
      }
    }
    if (is.null(fj.pch) || nrow(fj.pch) != nrow(res$fj)) {
      if (is.null(tepPlotInfo$fj.pch)) {
        fj.pch <- as.matrix(rep(21, nrow(res$fj)))
      }
      else {
        fj.pch <- tepPlotInfo$fj.pch
      }
    }
    if (is.null(constraints)) {
      if (!is.null(tepPlotInfo$constraints)) {
        constraints <- tepPlotInfo$constraints
      }
      constraints <- calculateConstraints(results = res, 
                                          x_axis = x_axis, y_axis = y_axis, constraints = constraints)
    }
    fj.boot.tests <- rowSums(inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[, 
                                                                                         c(x_axis, y_axis)])
    fj.no.boot.axes <- which(fj.boot.tests == 0)
    fj.both.boot.axes <- which(fj.boot.tests == 2)
    fj.x.boot.axis <- which((inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[, 
                                                                                         c(x_axis)] - inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[, 
                                                                                                                                                                  c(y_axis)]) == 1)
    fj.y.boot.axis <- which((inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[, 
                                                                                         c(y_axis)] - inference.info$boot.data$fj.boot.data$tests$sig.boot.ratios[, 
                                                                                                                                                                  c(x_axis)]) == 1)
    fi.boot.tests <- rowSums(inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[, 
                                                                                         c(x_axis, y_axis)])
    fi.no.boot.axes <- which(fi.boot.tests == 0)
    fi.both.boot.axes <- which(fi.boot.tests == 2)
    fi.x.boot.axis <- which((inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[, 
                                                                                         c(x_axis)] - inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[, 
                                                                                                                                                                  c(y_axis)]) == 1)
    fi.y.boot.axis <- which((inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[, 
                                                                                         c(y_axis)] - inference.info$boot.data$fi.boot.data$tests$sig.boot.ratios[, 
                                                                                                                                                                  c(x_axis)]) == 1)
    orig.fi.col <- fi.col
    if (length(fj.no.boot.axes) != 0) {
      fj.col[fj.no.boot.axes, 1] <- "gray"
    }
    if (length(fi.no.boot.axes) != 0) {
      fi.col[fi.no.boot.axes, 1] <- "gray"
    }
    if (color.by.boots) {
      if (is.null(boot.cols) || length(boot.cols) != 3) {
        boot.cols <- c("plum4", "darkseagreen", "firebrick3")
      }
      if (length(fj.x.boot.axis) != 0) {
        fj.col[fj.x.boot.axis, 1] <- boot.cols[1]
      }
      if (length(fi.x.boot.axis) != 0) {
        fi.col[fi.x.boot.axis, 1] <- boot.cols[1]
      }
      if (length(fj.y.boot.axis) != 0) {
        fj.col[fj.y.boot.axis, 1] <- boot.cols[2]
      }
      if (length(fi.y.boot.axis) != 0) {
        fi.col[fi.y.boot.axis, 1] <- boot.cols[2]
      }
      if (length(fj.both.boot.axes) != 0) {
        fj.col[fj.both.boot.axes, 1] <- boot.cols[3]
      }
      if (length(fi.both.boot.axes) != 0) {
        fi.col[fi.both.boot.axes, 1] <- boot.cols[3]
      }
      fj.col.y <- fj.col.x <- fj.col
      fj.col.y[fj.x.boot.axis, 1] <- "gray"
      fj.col.x[fj.y.boot.axis, 1] <- "gray"
      fi.col.y <- fi.col.x <- fi.col
      fi.col.y[fi.x.boot.axis, 1] <- "gray"
      fi.col.x[fi.y.boot.axis, 1] <- "gray"
    }
    fii.plot.info <- prettyPlot(res$fii, x_axis = x_axis, 
                                y_axis = y_axis, col = fii.col, pch = fii.pch, axes = TRUE, 
                                xlab = xlab, ylab = ylab, main = main, constraints = constraints, 
                                contributionCircles = FALSE, dev.new = F)
    fi.plot.info <- prettyPlot(res$fi, x_axis = x_axis, y_axis = y_axis, 
                               col = orig.fi.col, pch = fi.pch, axes = FALSE, contributionCircles = TRUE, 
                               contributions = abs(inference.info$boot.data$fi.boot.data$tests$boot.ratios), 
                               dev.new = FALSE, new.plot = FALSE)
    if (showHulls > 0 && showHulls <= 1) {
      colorDesign <- makeNominalData(fii.col)
      for (i in 1:nrow(res$fi)) {
        peeledHull(res$fii[which(fii.col[, 1] == orig.fi.col[i, 
                                                             1]), ], x_axis = x_axis, y_axis = y_axis, percentage = showHulls, 
                   col = "black", lwd = 3)
        peeledHull(res$fii[which(fii.col[, 1] == orig.fi.col[i, 
                                                             1]), ], x_axis = x_axis, y_axis = y_axis, percentage = showHulls, 
                   col = orig.fi.col[i, ], lwd = 1)
      }
    }
    fi.plot.info <- prettyPlot(res$fi, x_axis = x_axis, y_axis = y_axis, 
                               col = fi.col, pch = fi.pch, axes = TRUE, xlab = xlab, 
                               ylab = ylab, main = main, constraints = constraints, 
                               contributionCircles = TRUE, contributions = abs(inference.info$boot.data$fi.boot.data$tests$boot.ratios), 
                               dev.new = F)
    if (showHulls > 0 && showHulls <= 1) {
      colorDesign <- makeNominalData(fii.col)
      for (i in 1:nrow(res$fi)) {
        boot.items <- t(inference.info$boot.data$fi.boot.data$boots[i, 
                                                                    , ])
        peeledHull(boot.items, x_axis = x_axis, y_axis = y_axis, 
                   percentage = showHulls, col = "black", lwd = 5)
        peeledHull(boot.items, x_axis = x_axis, y_axis = y_axis, 
                   percentage = showHulls, col = fi.col[i, ], 
                   lwd = 3)
      }
    }
    if (biplots) {
      fj.plot.info <- prettyPlot(res$fj, x_axis = x_axis, 
                                 y_axis = y_axis, col = fj.col, pch = fj.pch, 
                                 axes = FALSE, contributionCircles = TRUE, contributions = abs(inference.info$boot.data$fj.boot.data$tests$boot.ratios), 
                                 dev.new = FALSE, new.plot = FALSE)
    }
    else {
      fj.plot.info <- prettyPlot(res$fj, x_axis = x_axis, 
                                 y_axis = y_axis, col = fj.col, pch = fj.pch, 
                                 axes = TRUE, xlab = xlab, ylab = ylab, main = main, 
                                 constraints = constraints, contributionCircles = TRUE, 
                                 contributions = abs(inference.info$boot.data$fj.boot.data$tests$boot.ratios), 
                                 dev.new = F)
    }
    if (bootstrapBars) {
      prettyBars(inference.info$boot.data$fi.boot.data$tests$boot.ratios, 
                 axis = x_axis, fg.col = fi.col.x, dev.new = F, 
                 threshold.line = TRUE, main = paste("Bootstrap Ratios Component: ", 
                                                     x_axis, sep = ""), bg.lims = c(-inference.info$boot.data$fi.boot.data$tests$critical.value, 
                                                                                    inference.info$boot.data$fi.boot.data$tests$critical.value))
      prettyBars(inference.info$boot.data$fi.boot.data$tests$boot.ratios, 
                 axis = y_axis, fg.col = fi.col.y, dev.new = F, 
                 horiz = FALSE, threshold.line = TRUE, main = paste("Bootstrap Ratios Component: ", 
                                                                    y_axis, sep = ""), bg.lims = c(-inference.info$boot.data$fi.boot.data$tests$critical.value, 
                                                                                                   inference.info$boot.data$fi.boot.data$tests$critical.value))
      prettyBars(inference.info$boot.data$fj.boot.data$tests$boot.ratios, 
                 axis = x_axis, fg.col = fj.col.x, dev.new = F, 
                 threshold.line = TRUE, main = paste("Bootstrap Ratios Component: ", 
                                                     x_axis, sep = ""), bg.lims = c(-inference.info$boot.data$fj.boot.data$tests$critical.value, 
                                                                                    inference.info$boot.data$fj.boot.data$tests$critical.value))
      prettyBars(inference.info$boot.data$fj.boot.data$tests$boot.ratios, 
                 axis = y_axis, fg.col = fj.col.y, dev.new = F, 
                 horiz = FALSE, threshold.line = TRUE, main = paste("Bootstrap Ratios Component: ", 
                                                                    y_axis, sep = ""), bg.lims = c(-inference.info$boot.data$fj.boot.data$tests$critical.value, 
                                                                                                   inference.info$boot.data$fj.boot.data$tests$critical.value))
    }
    if (correlationPlotter && class(res)[1] %in% pca.types) {
      correlationPlotter(res$X, res$fi, col = fj.col, x_axis = 1, 
                         y_axis = 2, xlab = xlab, ylab = ylab, main = main, dev.new=F)
    }
  }
}