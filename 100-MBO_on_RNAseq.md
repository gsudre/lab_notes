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
nfolds = 5

# bw
ncores = 32
niters = 20 # to run in MBO
init_design = 100 # number of starting parameter combinations
nrounds = 400  # for xgb.cv
nstop = 40  # for xgb.cv

control = makeMBOControl()
control = setMBOControlTermination(control, iters = niters)

set.seed(myseed)
my_params = makeParamSet(
        makeNumericParam("eta",              lower = 0.001, upper = 0.5),
        makeNumericParam("gamma",            lower = 0,     upper = 7),
        makeIntegerParam("max_depth",        lower= 1,      upper = 5),
        makeIntegerParam("min_child_weight", lower= 1,      upper = 10),
        makeNumericParam("subsample",        lower = 0.2,   upper = 1),
        makeNumericParam("colsample_bytree", lower = 0.2,   upper = 1),
        makeNumericParam("lambda", lower = 1,   upper = 10),
        makeNumericParam("alpha", lower = 0,   upper = 10)
      )
des = generateDesign(n=init_design, par.set = my_params)

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

    dtrain <- xgb.DMatrix(data = as.matrix(X_train),
                          label = as.numeric(y_train)-1, missing=NA)
    dtest <- xgb.DMatrix(data = as.matrix(X_test),
                         label = as.numeric(y_test)-1, missing=NA)

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
      par.set = my_params,
      minimize = FALSE
    )

    run = mbo(fun = obj.fun,
              design = des,
              control = control,
              show.info = TRUE)

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

Running overnight...

# 2020-04-10 06:38:31

Our ys are looking good, but I still only have 34 out of 56 predictions after a
whole night. Maybe that's how it is going to be, but let's see if we can
optimize it a bit more. I'll run Caudate for now, and then I still need to add a
dummy for region.

```r
myregion = 'Caudate'
fname = sprintf('~/data/rnaseq_derek/X_%snoPH_zv_nzv_center_scale.rds', myregion)

X = readRDS(fname)
just_target = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
y = just_target[just_target$Region==myregion, 'Diagnosis']

library(caret)
library(xgboost)

library(mlrMBO)

myseed = 42
nfolds = 5

# bw
ncores = 32
niters = 20 # to run in MBO
init_design = 50 # number of starting parameter combinations
nrounds = 400  # for xgb.cv
nstop = 10  # for xgb.cv

control = makeMBOControl()
control = setMBOControlTermination(control, iters = niters)

set.seed(myseed)
my_params = makeParamSet(
        makeNumericParam("eta",              lower = 0.001, upper = 0.5),
        makeNumericParam("gamma",            lower = 0,     upper = 7),
        makeIntegerParam("max_depth",        lower= 1,      upper = 5),
        makeIntegerParam("min_child_weight", lower= 1,      upper = 6),
        makeNumericParam("subsample",        lower = 0.4,   upper = 1),
        makeNumericParam("colsample_bytree", lower = 0.4,   upper = 1),
        makeNumericParam("lambda", lower = 1,   upper = 10),
        makeNumericParam("alpha", lower = 0,   upper = 3)
      )
des = generateDesign(n=init_design, par.set = my_params)

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

    dtrain <- xgb.DMatrix(data = as.matrix(X_train),
                          label = as.numeric(y_train)-1, missing=NA)
    dtest <- xgb.DMatrix(data = as.matrix(X_test),
                         label = as.numeric(y_test)-1, missing=NA)

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
      par.set = my_params,
      minimize = FALSE
    )

    run = mbo(fun = obj.fun,
              design = des,
              control = control,
              show.info = TRUE)

    set.seed(myseed)
    xg_mod <- xgboost(data = dtrain, params = run$x, nround = nrounds,
                      verbose = F, nthread=ncores)
    y_probs = c(y_probs, predict(xg_mod, dtest))
    tmp = data.frame(run$x, y=run$y)
    best_params = rbind(best_params, tmp)
    print(best_params)

    # updated LOOCV predictions
    y_preds = factor(ifelse (y_probs > 0.5, levels(y)[2], levels(y)[1]),
                     levels=levels(y))
    dat = data.frame(obs = y[1:length(y_probs)], pred = y_preds, junk=y_probs)
    colnames(dat)[3] = levels(y)[2]
    print(twoClassSummary(dat, lev=levels(y)))
}
```

And we'll need some modifications to get it working for both regions at the same
time:

```r
fname = '~/data/rnaseq_derek/X_noPH_zv_nzv_center_scale.rds'

X = readRDS(fname)
just_target = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
y = just_target[, 'Diagnosis']
```

Maybe I should script that out so I'm not running it in interactive mode all the
time?

```bash
my_dir=~/data/rnaseq_derek/
cd $my_dir
my_script=~/research_code/rnaseq_LOOCV_MBO.R;
out_file=swarm.loocv
ncores=32;
rm $out_file
for r in ACC Caudate both; do
    for m in auc error; do
        res_file=$my_dir/loocv_${r}_${m}.RData;
        echo "Rscript $my_script $r $ncores $m $res_file;" >> $out_file;
    done;
done

swarm -g 30 -t $ncores --job-name loocv --time 24:00:00 -f $out_file \
    -m R --partition norm --logdir trash
```

# 2020-04-11 09:30:58

I noticed I wasn't setting the metric and objective explicitly in the final
model before using test data, and also I wasn't optimizing nrounds. So, let's run
it again:

```bash
my_dir=~/data/rnaseq_derek/
cd $my_dir
my_script=~/research_code/rnaseq_LOOCV_MBO.R;
out_file=swarm.loocv
ncores=32;
rm $out_file
for r in ACC Caudate both; do
    for m in auc error; do
        res_file=$my_dir/loocv_${r}_${m}_v2.RData;
        echo "Rscript $my_script $r $ncores $m $res_file;" >> $out_file;
    done;
done

swarm -g 30 -t $ncores --job-name loocv2 --time 24:00:00 -f $out_file \
    -m R --partition norm --logdir trash
```

The last two are still running, but we can take a look at what we found so far.

```r
> library(caret)
Loading required package: lattice
Loading required package: ggplot2
> load('loocv_ACC_auc_v2.RData')
> dat$Case = 1-dat$Control
> print(twoClassSummary(dat, lev=levels(dat$obs)))
      ROC      Sens      Spec 
0.4787097 0.3600000 0.6129032 
> load('loocv_ACC_error_v2.RData')
> dat$Case = 1-dat$Control
> print(twoClassSummary(dat, lev=levels(dat$obs)))
      ROC      Sens      Spec 
0.4541935 0.3200000 0.6451613 
> load('loocv_Caudate_auc_v2.RData')
> dat$Case = 1-dat$Control
> print(twoClassSummary(dat, lev=levels(dat$obs)))
      ROC      Sens      Spec 
0.4860140 0.3461538 0.5757576 
> load('loocv_Caudate_error_v2.RData')
> dat$Case = 1-dat$Control
> print(twoClassSummary(dat, lev=levels(dat$obs)))
      ROC      Sens      Spec 
0.5093240 0.4615385 0.6060606 
>
```

not doing great, even though our training AUC is getting around .8...

We might need to try simpler classifiers.



# TODO
