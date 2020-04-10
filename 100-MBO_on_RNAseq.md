# 2020-04-09 20:20:44

Following our encouraging results with Derek's data using xgbTree, and my recent
developments with MBO, let's pair the two:

```r
myregion = 'ACC'
fname = sprintf('~/data/rnaseq_derek/X_%snoPH_zv_nzv_center_scale.rds', myregion)
ncores = 32
nfolds = 10

X = readRDS(fname)
just_target = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
y = just_target[just_target$Region==myregion, 'Diagnosis']

library(caret)
library(xgboost)

library(mlrMBO)

myseed = 42
nfolds = 10

# bw
ncores = 32
niters = 20 # to run in MBO
init_design = 100 # number of starting parameter combinations
nrounds = 400  # for xgb.cv
nstop = 40  # for xgb.cv

control = makeMBOControl()
control = setMBOControlTermination(control, iters = niters)

y_probs = c()
best_params = c()
for (test_row in 1:nrow(X)) {
    train_rows = setdiff(1:nrow(X), test_row)
    X_train <- X[train_rows, ]
    X_test <- X[-train_rows, ]
    y_train <- y[train_rows]
    y_test <- y[-train_rows]

    print(sprintf('LOOCV %d / %d', test_row, nrow(X)))
    set.seed(myseed)
    cv_folds = createFolds(y_train, k = nfolds)

    dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = as.numeric(y_train)-1, missing=NA)
    dtest <- xgb.DMatrix(data = t(as.matrix(X_test)), label = as.numeric(y_test)-1, missing=NA)

    obj.fun  <- smoof::makeSingleObjectiveFunction(
    name = "xgb_cv_bayes",
    fn =   function(x){
      set.seed(myseed)
      cv <- xgb.cv(params = list(
                                booster          = "gbtree",
                                eta              = x["eta"],
                                max_depth        = x["max_depth"],
                                min_child_weight = x["min_child_weight"],
                                gamma            = x["gamma"],
                                subsample        = x["subsample"],
                                colsample_bytree = x["colsample_bytree"],
                                lambda = x["lambda"],
                                alpha = x["alpha"],
                                objective        = 'binary:logistic',
                                eval_metric     = "auc"),
                                data = dtrain,
                                nround = nrounds,
                                folds=  cv_folds,
                                prediction = FALSE,
                                showsd = TRUE,
                                early_stopping_rounds = nstop,
                                verbose = 0,
                                nthread=ncores)
        cv$evaluation_log[, max(test_auc_mean)]
      },
      par.set = makeParamSet(
        makeNumericParam("eta",              lower = 0.001, upper = 0.5),
        makeNumericParam("gamma",            lower = 0,     upper = 7),
        makeIntegerParam("max_depth",        lower= 1,      upper = 5),
        makeIntegerParam("min_child_weight", lower= 1,      upper = 10),
        makeNumericParam("subsample",        lower = 0.2,   upper = 1),
        makeNumericParam("colsample_bytree", lower = 0.2,   upper = 1),
        makeNumericParam("lambda", lower = 1,   upper = 10),
        makeNumericParam("alpha", lower = 0,   upper = 10)
      ),
      minimize = FALSE
    )

    des = generateDesign(n=init_design,
                        par.set = getParamSet(obj.fun))
    run = mbo(fun = obj.fun,
              design = des,
              control = control,
              show.info = FALSE)

    xg_mod <- xgboost(data = dtrain, params = run$x, nround = nrounds,
                      verbose = F, nthread=ncores)
    y_probs = c(y_probs, predict(xg_mod, dtest))
    tmp = data.frame(run$x, y=run$y)
    best_params = rbind(best_params, tmp)
    print(tmp)
}

y_preds = factor(ifelse (y_probs > 0.5, c2, c1), levels=levels(y))
dat = data.frame(obs = y_test, pred = y_preds, junk=y_probs)
colnames(dat)[3] = c1
twoClassSummary(dat, lev=levels(y_test))
```

Waiting for free nodes to run this...

