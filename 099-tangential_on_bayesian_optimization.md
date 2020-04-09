# 2020-04-08 11:48:27

Let's see if we can narrow down a parameter using Bayesian optimization.

```r
library(caret)
library(xgboost)

fname = '~/data/baseline_prediction/prs_start/gf_philip_03292020.csv'
phen = 'categ_all.4'
c1 = 'improvers'
c2 = 'stable_symptomatic'

data = read.csv(fname)
data$sex_numeric = as.factor(data$sex_numeric)
data$SES_group3 = as.factor(data$SES_group3)
data$slf_fa = data$slf_all  # just to make it easier to filter out
var_names = c(# PRS
              'ADHD_PRS0.000100.orig', 'ADHD_PRS0.001000.orig',
              'ADHD_PRS0.010000.orig', 'ADHD_PRS0.050000.orig',
              'ADHD_PRS0.100000.orig', 'ADHD_PRS0.200000.orig',
              'ADHD_PRS0.300000.orig', 'ADHD_PRS0.400000.orig',
              'ADHD_PRS0.500000.orig',
              # DTI
              'atr_fa', 'cst_fa', 'cing_cing_fa', 'cing_hipp_fa', 'cc_fa',
              'ilf_fa', 'slf_fa', 'unc_fa',
              # demo
              'sex_numeric', 'SES_group3',
              # cog
              'FSIQ', 'SS_RAW', 'DS_RAW', 'PS_RAW', 'VMI.beery_RAW',
              # anat
              'cerebellum_white', 'cerebellum_grey', 'amygdala',
              'cingulate', 'lateral_PFC', 'OFC', 'striatum', 'thalamus'
              )

# data2 = data[, c(var_names, 'FAMID')]
data2 = data[, var_names]

data2$phen = as.factor(data[, phen])
# data2$FAMID = as.factor(data$FAMID)
dummies = dummyVars(phen ~ ., data = data2)
data3 = predict(dummies, newdata = data2)

# selecting only kids in the 2 specified groups
keep_me = data2$phen==c1 | data2$phen==c2
X = data3[keep_me, ]
y = factor(data2[keep_me, 'phen'])

test_row = 1
train_rows = setdiff(1:nrow(X), test_row)
X_train <- X[train_rows, ]
X_test <- X[-train_rows, ]
y_train <- y[train_rows]
y_test <- y[-train_rows]

dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = as.numeric(y_train)-1, missing=NA)
dtest <- xgb.DMatrix(data = t(as.matrix(X_test)), label = as.numeric(y_test)-1, missing=NA)
```

```r
library(rBayesianOptimization)
set.seed(42)
cv_folds = createFolds(y, k = 10)
xgb_cv_bayes <- function(nround, max.depth, min_child_weight, subsample, eta,
                         gamma,colsample_bytree,max_delta_step,verbose) {
    param<-list(booster = "gbtree",
                max_depth = max.depth,
                min_child_weight = min_child_weight,
                eta=eta,gamma=gamma,
                subsample = subsample,
                colsample_bytree = colsample_bytree,
                max_delta_step=max_delta_step,
                lambda = 1, alpha = 0,
                objective = "binary:logistic",
                eval_metric = "auc")
    cv <- xgb.cv(params = param, data = dtrain, folds = cv_folds,
                 nrounds = 1000, early_stopping_rounds = 10, maximize = TRUE,
                 verbose = T)
    list(Score = cv$evaluation_log$test_auc_mean[cv$best_iteration],
         Pred=cv$best_iteration)
    # we don't need cross-validation prediction and we need the number of rounds.
    # a workaround is to pass the number of rounds(best_iteration) to the Pred, which is a default parameter in the rbayesianoptimization library.
}

OPT_Res <- BayesianOptimization(xgb_cv_bayes,
                                bounds = list(max.depth =c(3L, 10L),
                                min_child_weight = c(1L, 40L),
                                subsample = c(0.6, 0.9),
                                eta=c(0.01,0.3),
                                gamma = c(0.0, 0.2),
                                colsample_bytree=c(0.5,0.8),
                                max_delta_step=c(1L,10L)),
                                init_grid_dt = NULL, init_points = 10,
                                n_iter = 10, acq = "ucb", kappa = 2.576,
                                eps = 0.0)

best_param <- list(
    booster = "gbtree",
    eval.metric = "auc",
    objective = "binary:logistic",
    max_depth = OPT_Res$Best_Par["max.depth"],
    eta = OPT_Res$Best_Par["eta"],
    gamma = OPT_Res$Best_Par["gamma"],
    subsample = OPT_Res$Best_Par["subsample"],
    colsample_bytree = OPT_Res$Best_Par["colsample_bytree"],
    min_child_weight = OPT_Res$Best_Par["min_child_weight"],
    max_delta_step = OPT_Res$Best_Par["max_delta_step"])

nrounds=OPT_Res$Pred[[which.max(OPT_Res$History$Value)]]
xgb_model <- xgb.train (params = best_param, data = dtrain, nrounds = nrounds)
```

Nahh... this code is broken. 

Let's try this one instead: https://www.simoncoulombe.com/2019/01/bayesian/

```r
library(mlrMBO)

set.seed(42)
cv_folds = createFolds(y_train, k = 10)
obj.fun  <- smoof::makeSingleObjectiveFunction(
  name = "xgb_cv_bayes",
  fn =   function(x){
    set.seed(12345)
    cv <- xgb.cv(params = list(
      booster          = "gbtree",
      eta              = x["eta"],
      max_depth        = x["max_depth"],
      min_child_weight = x["min_child_weight"],
      gamma            = x["gamma"],
      subsample        = x["subsample"],
      colsample_bytree = x["colsample_bytree"],
      objective        = 'binary:logistic', 
      eval_metric     = "auc"),
      data = dtrain,
      nround = 30,
      folds=  cv_folds,
      prediction = FALSE,
      showsd = TRUE,
      early_stopping_rounds = 10,
      verbose = 0)
    cv$evaluation_log[, max(test_auc_mean)]
  },
  par.set = makeParamSet(
    makeNumericParam("eta",              lower = 0.001, upper = 0.05),
    makeNumericParam("gamma",            lower = 0,     upper = 5),
    makeIntegerParam("max_depth",        lower= 1,      upper = 10),
    makeIntegerParam("min_child_weight", lower= 1,      upper = 10),
    makeNumericParam("subsample",        lower = 0.2,   upper = 1),
    makeNumericParam("colsample_bytree", lower = 0.2,   upper = 1)
  ),
  minimize = FALSE
)

# generate an optimal design with only 10  points
des = generateDesign(n=100,
                     par.set = getParamSet(obj.fun), 
                     fun = lhs::randomLHS)  ## . If no design is given by the user, mlrMBO will generate a maximin Latin Hypercube Design of size 4 times the number of the black-box function’s parameters.

# bayes will have 10 additional iterations
control = makeMBOControl()
control = setMBOControlTermination(control, iters = 10)

# run this!
run = mbo(fun = obj.fun, 
          design = des,  
          control = control, 
          show.info = TRUE)

run <- read_rds( "temp_files/run.rds")
# print a summary with run
#run
# return  best model hyperparameters using run$x
# return best log likelihood using run$y
# return all results using run$opt.path$env$path
run$opt.path$env$path  %>% 
  mutate(Round = row_number()) %>%
  mutate(type = case_when(
    Round==1  ~ "1- hardcoded",
    Round<= 11 ~ "2 -design ",
    TRUE ~ "3 - mlrMBO optimization")) %>%
  ggplot(aes(x= Round, y= y, color= type)) + 
  geom_point() +
  labs(title = "mlrMBO optimization")+
  ylab("loglikelihood")
```

OK, the error was that I was KFolding y, instead of y_train!

In any case, this is working. It's just a matter of giving it a good enough
range of parameters, and figuring out a bit how to best use this package, in
terms of multicore and other parameters.

# 2020-04-09 06:43:49

PLaying with this a bit more, we get:

```r
library(mlrMBO)

myseed = 42
ncores = 2
nfolds = 10
niters = 20 # to run in MBO
init_design = 100 # number of starting parameter combinations

set.seed(myseed)
cv_folds = createFolds(y_train, k = nfolds)
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
      nround = 1000,
      folds=  cv_folds,
      prediction = FALSE,
      showsd = TRUE,
      early_stopping_rounds = 50,
      verbose = 0)
    cv$evaluation_log[, max(test_auc_mean)]
  },
  par.set = makeParamSet(
    makeNumericParam("eta",              lower = 0.001, upper = 0.5),
    makeNumericParam("gamma",            lower = 0,     upper = 5),
    makeIntegerParam("max_depth",        lower= 1,      upper = 5),
    makeIntegerParam("min_child_weight", lower= 1,      upper = 10),
    makeNumericParam("subsample",        lower = 0.2,   upper = 1),
    makeNumericParam("colsample_bytree", lower = 0.2,   upper = 1),
    makeNumericParam("lambda", lower = 1,   upper = 10),
    makeNumericParam("alpha", lower = 0,   upper = 10)
  ),
  minimize = FALSE
)

# generate an optimal design with only 10  points
des = generateDesign(n=init_design,
                     par.set = getParamSet(obj.fun))

library(parallelMap)
parallelStartMulticore(cpus = ncores, show.info = FALSE)
# bayes will have 10 additional iterations
control = makeMBOControl(propose.points = ncores)
control = setMBOControlTermination(control, iters = niters)

# run this!
run = mbo(fun = obj.fun, 
          design = des,  
          control = control, 
          show.info = TRUE)
parallelStop()
```

Now I should probably save the best parameters for each CV so that we can have a
better range later...

```r
a = data.frame(run$x, y=run$y)
```

So, in a CV setting it'd be something like this:

```r
library(mlrMBO)
library(parallelMap)

myseed = 42
nfolds = 10

# laptop
ncores = 2
niters = 5 # to run in MBO
init_design = 10 # number of starting parameter combinations
nrounds = 100  # for xgb.cv
nstop = 5  # for xgb.cv

# bw
ncores = 32
niters = 20 # to run in MBO
init_design = 100 # number of starting parameter combinations
nrounds = 1000  # for xgb.cv
nstop = 10  # for xgb.cv

parallelStartMulticore(cpus = ncores, show.info = FALSE)
control = makeMBOControl(propose.points = ncores)
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
                                verbose = 0)
        cv$evaluation_log[, max(test_auc_mean)]
      },
      par.set = makeParamSet(
        makeNumericParam("eta",              lower = 0.001, upper = 0.5),
        makeNumericParam("gamma",            lower = 0,     upper = 5),
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
                      verbose = F)
    y_probs = c(y_probs, predict(xg_mod, dtest))
    tmp = data.frame(run$x, y=run$y)
    best_params = rbind(best_params, tmp)
    print(tmp)
}
parallelStop()

y_preds = factor(ifelse (y_probs > 0.5, c2, c1), levels=levels(y))
dat = data.frame(obs = y_test, pred = y_preds, junk=y_probs)
colnames(dat)[3] = c1
twoClassSummary(dat, lev=levels(y_test))
```

I think I'm ready to run some interactive sessions here. This might not be the
best setup for parallel processing though, based on how much the CPUs are being
used. What if I hyperthread xgboost instead?




# TODO
* check that probability assignments to the classes is right