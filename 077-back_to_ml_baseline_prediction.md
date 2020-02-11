# 2020-02-05 14:33:09

Let's give it another try for ML in baseline prediction. Now that I have a set
of data I'm working with (see 076), it should be easier to find a nicer model.

I'll try not to residualize or impute anything. At least not at first.

```r
data = readRDS('~/data/baseline_prediction/prs_start/complete_massaged_data_02032020.rds')
data$externalizing = as.factor(data$externalizing)
min_sx = 6
for (sx in c('inatt', 'hi')) {
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data[, phen] = 'notGE6adhd'
    my_nvs = which(is.na(data[, phen_slope]))
    idx = data[my_nvs, 'base_inatt'] <= 2 & data[my_nvs, 'base_hi'] <= 2
    data[my_nvs[idx], phen] = 'nv012'
    data[which(data[, phen_slope] < thresh), phen] = 'imp'
    data[which(data[, phen_slope] >= thresh), phen] = 'nonimp'
    data[, phen] = factor(data[, phen], ordered=F)
    data[, phen] = relevel(data[, phen], ref='nv012')
    ophen = sprintf('ORDthresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data[, ophen] = factor(data[, phen],
                         levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                         ordered=T)
}
```

Just for kicks I'll start with the clinical group first. If anything, it'll take
less time to train:

```r
phen = 'thresh0.00_inatt_GE6_wp05' # NVs don't change across SX
adhd = data[, phen] == 'nonimp' | data[, phen] == 'imp'
data2 = data[adhd, ]
data2[, phen] = factor(data2[, phen], ordered=F)
data2[, phen] = relevel(data2[, phen], ref='imp')
training = data2[data2$bestInFamily, ]
testing = data2[!data2$bestInFamily, ]
```

```r
set.seed(42)
fitControl <- trainControl(## 10-fold CV
                           method = "repeatedcv",
                           number = 10,
                           ## repeated ten times
                           repeats = 10,
                           classProbs=T)

var_names = colnames(data)[grepl(colnames(data), pattern='ADHD')]
svmFit <- train(x = training[, var_names], y=training[, phen], 
                 method = "svmRadial", 
                 trControl = fitControl, 
                 preProc = c("center", "scale"),
                 tuneLength = 8,
                 metric = "ROC")
rdaFit <- train(x = training[, var_names], y=training[, phen], 
                 method = "rda", 
                 trControl = fitControl, 
                 preProc = c("center", "scale"),
                 tuneLength = 8,
                 metric = "ROC")
gmbFit <- train(x = training[, var_names], y=training[, phen], 
                 method = "gbm", 
                 trControl = fitControl, 
                 preProc = c("center", "scale"),
                 tuneLength = 8,
                 metric = "ROC")
resamps <- resamples(list(GBM = gmbFit,
                          SVM = svmFit,
                          RDA = rdaFit))
summary(resamps)

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(resamps, layout = c(3, 1))

trellis.par.set(caretTheme())
dotplot(resamps, metric = "Accuracy")
```

# 2020-02-06 09:09:41

OK, so let's play with methods that can handle NAs, which is one of the main
reasons to do this anyways:

```r
library(caret)

my_cols = c(38:49, 86:123, 125, 128:138)
phen = 'thresh0.00_inatt_GE6_wp05'
var_names = c(colnames(data)[my_cols], 'base_age', 'sex', 'base_inatt')
phen = 'thresh0.50_hi_GE6_wp05'
var_names = c(colnames(data)[my_cols], 'base_age', 'sex', 'base_hi')
models = c('AdaBoost.M1', 'AdaBag', 'ada', 'C5.0', 'rpart', 'rpart1SE',
           'rpart2', 'rpartScore', 'C5.0Cost', 'rpartCost', 'C5.0Rules',
           'C5.0Tree')
models = c('C5.0', 'rpart', 'rpart1SE',
           'rpart2', 'rpartScore', 'C5.0Cost', 'rpartCost', 'C5.0Rules',
           'C5.0Tree')

adhd = data[, phen] == 'nonimp' | data[, phen] == 'imp'
data2 = data[adhd, ]
data2[, phen] = factor(data2[, phen], ordered=F)
data2[, phen] = relevel(data2[, phen], ref='imp')
training = data2[data2$bestInFamily, ]
testing = data2[!data2$bestInFamily, ]
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           classProbs=T)
for (m in models) {
    print(m)
    model_str = sprintf('fit_%s <- train(x = training[, var_names],
                                         y=training[, phen],
                                         method = "%s",
                                         trControl = fitControl,
                                         tuneLength = 10,
                                         metric = "AUC")', m, m)
    eval(parse(text=model_str))
}
```

Note that it might be better to train within-modality classifiers first just so
we can deal with the covariates properly, instead of using PC covariates with
brain predictors, for example. If we go that route, then we'd have to use
stacked/voting classifiers on top of everything, especially one that takes into
consideration NAs.

# 2020-02-07 10:43:25

AdaBoost and AdaBag took a whole night and still weren't finished. I'll not play
with those for now, until I can actually optimize them.

```r
fits_str = sapply(models, function(m) sprintf('%s = fit_%s', m, m))
rs_str = paste('resamps <- resamples(list(',
               paste(fits_str, collapse=','),
               '))', sep="")
eval(parse(text=rs_str))
summary(resamps)

theme1 <- trellis.par.get()
theme1$plot.symbol$col = rgb(.2, .2, .2, .4)
theme1$plot.symbol$pch = 16
theme1$plot.line$col = rgb(1, 0, 0, .7)
theme1$plot.line$lwd <- 2
trellis.par.set(theme1)
bwplot(resamps, layout = c(3, 1))
```

![](images/2020-02-07-11-16-53.png)

The results haven't been overwhelming. But I think it's a better approach to do
a per-domain classifier first, and then promote those classifiers later. So, the
approach will be:

 * split data into train and test
 * within train, remove any entries with NAs and try the best possible
   classifier there. Make sure that when predicting new data is is NAs, we get
   NA probabilities. 
* train a classifier on top of those that combines the probabilities of the
  different modalities. Maybe weighted average, or even anoter stack.
* use the test data

A few questions arise, such as whether we should do it only for the ADHDs, or do
the 4 group classification. Let's do 4 groups for now, but it shouldn't be too
hard to change everything for binary later.

I'm going to do this with randomforests for now because they're quick, but we
could potentially have different classifiers per domain? Maybe not, too hard to
explain why we did it that way. But at least try other models across domains to
see what we get.

```r
phen = 'thresh0.00_inatt_GE6_wp05'
# phen = 'thresh0.50_hi_GE6_wp05'
model = 'rf'
sx = 'inatt'
adhd = data[, phen] == 'nonimp' | data[, phen] == 'imp'
data2 = data[adhd, ]
data2[, phen] = factor(data2[, phen], ordered=F)
data2[, phen] = relevel(data2[, phen], ref='imp')
training = data2[data2$bestInFamily, ]
testing = data2[!data2$bestInFamily, ]

set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           classProbs=T,
                           summaryFunction=twoClassSummary
                           )

# dti
var_names = colnames(data)[107:121]
keep_me = !is.na(training[var_names[1]])
this_data = training[keep_me, ]
dti_fit <- train(x = this_data[, var_names], y=this_data[, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
dti_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(dti_fit, type='prob')
dti_preds[keep_me, ] = preds

# anat
var_names = colnames(data)[96:106]
keep_me = !is.na(training[var_names[1]])
this_data = training[keep_me, ]
anat_fit <- train(x = this_data[, var_names], y=this_data[, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
anat_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(anat_fit, type='prob')
anat_preds[keep_me, ] = preds

# genomics
var_names = c(colnames(data)[38:49], colnames(data)[86:95])
keep_me = !is.na(training[var_names[1]])
this_data = training[keep_me, ]
dna_fit <- train(x = this_data[, var_names], y=this_data[, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
dna_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(dna_fit, type='prob')
dna_preds[keep_me, ] = preds

# neuropsych
var_names = c('FSIQ', "VMI.beery", "DS.wisc", "SSB.wisc", "SSF.wisc", "DS.wj",
              "VM.wj")
# I'll impute within neuropsych, so I don't have to train a classifier for each test
this_data = training[, var_names]
numNAvars = rowSums(is.na(this_data))
# but first remove anyone that only has one or no neuropsych
keep_me = numNAvars < 5
this_data = this_data[keep_me, ]
impute = preProcess(this_data, method = "bagImpute")
this_data <- predict(impute, this_data)
neuro_fit <- train(x = this_data[, var_names], y=training[keep_me, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
neuro_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(neuro_fit, type='prob')
neuro_preds[keep_me, ] = preds

# demo
var_names = c('base_age', 'sex', 'SES')
keep_me = !is.na(training[var_names[1]])
this_data = training[keep_me, ]
demo_fit <- train(x = this_data[, var_names], y=this_data[, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
demo_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(demo_fit, type='prob')
demo_preds[keep_me, ] = preds

# clinics
var_names = c('internalizing', 'externalizing', sprintf('base_%s', sx))
keep_me = !is.na(training[var_names[1]])
this_data = training[keep_me, ]
clin_fit <- train(x = this_data[, var_names], y=this_data[, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
clin_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(clin_fit, type='prob')
clin_preds[keep_me, ] = preds

# demo
var_names = c('base_age')
keep_me = !is.na(training[var_names[1]])
this_data = training[keep_me, ]
demo_fit <- train(x = this_data[, var_names], y=this_data[, phen],
                 method = model, trControl = fitControl, tuneLength = 10,
                 metric='ROC')
demo_preds = data.frame(imp=rep(NA, nrow(training)), nonimp=rep(NA, nrow(training)))
preds = predict(demo_fit, type='prob')
demo_preds[keep_me, ] = preds

# ensemble
prob_data = cbind(dti_preds[, 1], anat_preds[, 1], dna_preds[, 1],
                  neuro_preds[, 1], demo_preds[, 1], clin_preds[, 1])
colnames(prob_data) = c('dti', 'anat', 'dna', 'neuro', 'demo', 'clin')
ens_fit <- train(x = prob_data, y=training[, phen],
                 method = 'rpart2', trControl = fitControl, tuneLength = 10,
                 metric='ROC')
```

Now that I have an ensemble classifier, let's see how well it does on the test
data:

```r
for (dom in colnames(prob_data)) {
    print(dom)
    eval(parse(text=sprintf('keep_me = !is.na(testing[, colnames(%s_fit$trainingData)[1]])', dom)))
    this_data = testing[keep_me, ]
    eval(parse(text=sprintf('%s_test_preds = data.frame(imp=rep(NA, nrow(testing)), nonimp=rep(NA, nrow(testing)))', dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob", newdata=this_data)', dom)))
    eval(parse(text=sprintf('%s_test_preds[keep_me, ] = preds', dom)))
}
```

That almost worked, but neuropsych is still breaking. Also, I don't want to just
impute that one. Let's then modularize this to break up neuropsych...

```r
phen = 'thresh0.00_inatt_GE6_wp05'
# phen = 'thresh0.50_hi_GE6_wp05'
model = 'rf'
sx = 'inatt'
adhd = data[, phen] == 'nonimp' | data[, phen] == 'imp'
data2 = data[adhd, ]
data2[, phen] = factor(data2[, phen], ordered=F)
data2[, phen] = relevel(data2[, phen], ref='imp')
training = data2[data2$bestInFamily, ]
testing = data2[!data2$bestInFamily, ]

set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           classProbs=T,
                           summaryFunction=twoClassSummary
                           )

domains = list(iq_vmi = c('FSIQ', "VMI.beery"),
               wisc = c("SSB.wisc", "SSF.wisc", 'DSF.wisc', 'DSB.wisc'),
               wj = c("DS.wj", "VM.wj"),
               demo = c('base_age', 'sex', 'SES'),
               clin = c('internalizing', 'externalizing', sprintf('base_%s', sx)),
               gen = c(colnames(data)[38:49], colnames(data)[86:95]),
               dti = colnames(data)[107:121],
               anat = colnames(data)[96:106]
               )

for (dom in names(domains)) {
    print(sprintf('Training %s on %s (sx=%s)', dom, phen, sx))
    var_names = domains[[dom]]
    numNAvars = rowSums(is.na(training[, var_names]))
    keep_me = numNAvars == 0
    this_data = training[keep_me, ]
    eval(parse(text=sprintf('%s_fit <- train(x = this_data[, var_names],
                                             y=this_data[, phen],
                                             method = model,
                                             trControl = fitControl,
                                             tuneLength = 10, metric="ROC")',
                            dom)))
    eval(parse(text=sprintf('%s_preds = data.frame(imp=rep(NA, nrow(training)),
                                                   nonimp=rep(NA, nrow(training)))',
                            dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob")', dom)))
    eval(parse(text=sprintf('%s_preds[keep_me, ] = preds', dom)))
}
# ensemble
preds_str = sapply(names(domains), function(d) sprintf('%s_preds[, 1]', d))
cbind_str = paste('prob_data = cbind(', paste(preds_str, collapse=','), ')',
                  sep="")
eval(parse(text=cbind_str))
colnames(prob_data) = names(domains)
ens_fit <- train(x = prob_data, y=training[, phen],
                 method = 'rpart2', trControl = fitControl, tuneLength = 10,
                 metric='ROC')
```

Now we can do the testing again:

```r
library(pROC)
for (dom in names(domains)) {
    print(dom)
    eval(parse(text=sprintf('var_names = colnames(%s_fit$trainingData)', dom)))
    # ignore .outcome
    numNAvars = rowSums(is.na(testing[, var_names[1:(length(var_names)-1)]]))
    keep_me = numNAvars == 0
    this_data = testing[keep_me, ]
    eval(parse(text=sprintf('%s_test_preds = data.frame(imp=rep(NA, nrow(testing)), nonimp=rep(NA, nrow(testing)))', dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob", newdata=this_data)', dom)))
    eval(parse(text=sprintf('%s_test_preds[keep_me, ] = preds', dom)))
}
preds_str = sapply(names(domains), function(d) sprintf('%s_test_preds[, 1]', d))
cbind_str = paste('prob_test_data = cbind(', paste(preds_str, collapse=','), ')',
                  sep="")
eval(parse(text=cbind_str))
colnames(prob_test_data) = names(domains)
preds = predict(ens_fit, newdata=prob_test_data)
confusionMatrix(data = preds, reference = testing[, phen])
```

Before we start tuning this so we're not overfitting the training set, let's
play with the entire dataset to integrate all 4 groups:

```r
phen = 'thresh0.00_inatt_GE6_wp05'
# phen = 'thresh0.50_hi_GE6_wp05'
model = 'rf'
sx = 'inatt'
training = data[data$bestInFamily, ]
testing = data[!data$bestInFamily, ]

set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           classProbs=T,
                           summaryFunction=multiClassSummary
                           )

domains = list(iq_vmi = c('FSIQ', "VMI.beery"),
               wisc = c("SSB.wisc", "SSF.wisc", 'DSF.wisc', 'DSB.wisc'),
               wj = c("DS.wj", "VM.wj"),
               demo = c('base_age', 'sex', 'SES'),
               clin = c('internalizing', 'externalizing', sprintf('base_%s', sx)),
               gen = c(colnames(data)[38:49], colnames(data)[86:95]),
               dti = colnames(data)[107:121],
               anat = colnames(data)[96:106]
               )

for (dom in names(domains)) {
    print(sprintf('Training %s on %s (sx=%s)', dom, phen, sx))
    var_names = domains[[dom]]
    numNAvars = rowSums(is.na(training[, var_names]))
    keep_me = numNAvars == 0
    this_data = training[keep_me, ]
    eval(parse(text=sprintf('%s_fit <- train(x = this_data[, var_names],
                                             y=this_data[, phen],
                                             method = model,
                                             trControl = fitControl,
                                             tuneLength = 10, metric="AUC")',
                            dom)))
    eval(parse(text=sprintf('%s_preds = data.frame(nv012=rep(NA, nrow(training)),
                                                   imp=rep(NA, nrow(training)),
                                                   nonimp=rep(NA, nrow(training)),
                                                   notGE6adhd=rep(NA, nrow(training)))',
                            dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob")', dom)))
    eval(parse(text=sprintf('%s_preds[keep_me, ] = preds', dom)))
}
# ensemble
preds_str = sapply(names(domains), function(d) sprintf('%s_preds', d))
cbind_str = paste('prob_data = cbind(', paste(preds_str, collapse=','), ')',
                  sep="")
eval(parse(text=cbind_str))
prob_header = c()
for (dom in names(domains)) {
    for (g in colnames(preds)) {
        prob_header = c(prob_header, sprintf('%s_%s', dom, g))
    }
}
colnames(prob_data) = prob_header
ens_fit <- train(x = prob_data, y=training[, phen],
                 method = 'rpart2', trControl = fitControl, tuneLength = 10,
                 metric='AUC')
```

Then we modify the code for testing as well to take into account all 4 groups:

```r
for (dom in names(domains)) {
    print(dom)
    eval(parse(text=sprintf('var_names = colnames(%s_fit$trainingData)', dom)))
    # ignore .outcome
    numNAvars = rowSums(is.na(testing[, var_names[1:(length(var_names)-1)]]))
    keep_me = numNAvars == 0
    this_data = testing[keep_me, ]
    eval(parse(text=sprintf('%s_test_preds = data.frame(nv012=rep(NA, nrow(testing)),
                                                   imp=rep(NA, nrow(testing)),
                                                   nonimp=rep(NA, nrow(testing)),
                                                   notGE6adhd=rep(NA, nrow(testing)))', dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob", newdata=this_data)', dom)))
    eval(parse(text=sprintf('%s_test_preds[keep_me, ] = preds', dom)))
}
preds_str = sapply(names(domains), function(d) sprintf('%s_test_preds', d))
cbind_str = paste('prob_test_data = cbind(', paste(preds_str, collapse=','), ')',
                  sep="")
eval(parse(text=cbind_str))
colnames(prob_test_data) = prob_header
preds = predict(ens_fit, newdata=prob_test_data, type='prob')
multiclass.roc(testing[, phen], preds)
```

OK, so this code is now working. I can certainly play with kinds of metrics I'll
be evaluating, such as AUC or other, but at least training and predictions are
happening. Now, it's just a matter of tunning it.

## tunning the 4 group model

```r
library(caret)
library(pROC)
data = readRDS('~/data/baseline_prediction/prs_start/complete_massaged_data_02032020.rds')
data$externalizing = as.factor(data$externalizing)
phen = 'thresh0.00_inatt_GE6_wp05'
# phen = 'thresh0.50_hi_GE6_wp05'
model = 'lda'
sx = 'inatt'

min_sx = 6
for (sx in c('inatt', 'hi')) {
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data[, phen] = 'notGE6adhd'
    my_nvs = which(is.na(data[, phen_slope]))
    idx = data[my_nvs, 'base_inatt'] <= 2 & data[my_nvs, 'base_hi'] <= 2
    data[my_nvs[idx], phen] = 'nv012'
    data[which(data[, phen_slope] < thresh), phen] = 'imp'
    data[which(data[, phen_slope] >= thresh), phen] = 'nonimp'
    data[, phen] = factor(data[, phen], ordered=F)
    data[, phen] = relevel(data[, phen], ref='nv012')
    ophen = sprintf('ORDthresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data[, ophen] = factor(data[, phen],
                         levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                         ordered=T)
}
domains = list(iq_vmi = c('FSIQ', "VMI.beery"),
               wisc = c("SSB.wisc", "SSF.wisc", 'DSF.wisc', 'DSB.wisc'),
               wj = c("DS.wj", "VM.wj"),
               demo = c('base_age', 'sex', 'SES'),
               clin = c('internalizing', 'externalizing', sprintf('base_%s', sx)),
               gen = c(colnames(data)[38:49], colnames(data)[86:95]),
               dti = colnames(data)[107:121],
               anat = colnames(data)[96:106]
               )
set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10,
                           classProbs=T,
                           summaryFunction=multiClassSummary
                           )
training = data[data$bestInFamily, ]
testing = data[!data$bestInFamily, ]
for (dom in names(domains)) {
    print(sprintf('Training %s on %s (sx=%s, model=%s)', dom, phen, sx, model))
    var_names = domains[[dom]]
    numNAvars = rowSums(is.na(training[, var_names]))
    keep_me = numNAvars == 0
    this_data = training[keep_me, var_names]
    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        } else {
            this_data[, v] = as.numeric(this_data[, v])
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])
    eval(parse(text=sprintf('%s_fit <- train(x = this_data,
                                             y=training[keep_me, phen],
                                             method = model,
                                             trControl = fitControl,
                                             tuneLength = 10, metric="AUC")',
                            dom)))
    eval(parse(text=sprintf('%s_preds = data.frame(nv012=rep(NA, nrow(training)),
                                                   imp=rep(NA, nrow(training)),
                                                   nonimp=rep(NA, nrow(training)),
                                                   notGE6adhd=rep(NA, nrow(training)))',
                            dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob")', dom)))
    eval(parse(text=sprintf('%s_preds[keep_me, ] = preds', dom)))
}
# ensemble
preds_str = sapply(names(domains), function(d) sprintf('%s_preds', d))
cbind_str = paste('prob_data = cbind(', paste(preds_str, collapse=','), ')',
                  sep="")
eval(parse(text=cbind_str))
prob_header = c()
for (dom in names(domains)) {
    for (g in colnames(preds)) {
        prob_header = c(prob_header, sprintf('%s_%s', dom, g))
    }
}
colnames(prob_data) = prob_header
ens_fit <- train(x = prob_data, y=training[, phen],
                 method = 'rpart2', trControl = fitControl, tuneLength = 10,
                 metric='AUC')
print(ens_fit)
```

And let's see how we do in test data:

```r
for (dom in names(domains)) {
    print(dom)
    eval(parse(text=sprintf('var_names = colnames(%s_fit$trainingData)', dom)))
    # remove .outcome
    var_names = var_names[1:(length(var_names)-1)]
    numNAvars = rowSums(is.na(testing[, var_names]))
    keep_me = numNAvars == 0
    this_data = testing[keep_me, var_names]
    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        } else {
            this_data[, v] = as.numeric(this_data[, v])
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])
    eval(parse(text=sprintf('%s_test_preds = data.frame(nv012=rep(NA, nrow(testing)),
                                                   imp=rep(NA, nrow(testing)),
                                                   nonimp=rep(NA, nrow(testing)),
                                                   notGE6adhd=rep(NA, nrow(testing)))', dom)))
    eval(parse(text=sprintf('preds = predict(%s_fit, type="prob", newdata=this_data)', dom)))
    eval(parse(text=sprintf('%s_test_preds[keep_me, ] = preds', dom)))
}
preds_str = sapply(names(domains), function(d) sprintf('%s_test_preds', d))
cbind_str = paste('prob_test_data = cbind(', paste(preds_str, collapse=','), ')',
                  sep="")
eval(parse(text=cbind_str))
colnames(prob_test_data) = prob_header
preds = predict(ens_fit, newdata=prob_test_data, type='prob')
multiclass.roc(testing[, phen], preds)
```

This wasn't terrible... .77 ROC for the 4 group classification inattention. But
that's average ROC over all groups, so it can be a bit misleading. And we have
.62 for HI.

Let's add medication status, then we can go back to tuning the best model and
using just the clinical groups.

# TODO
* add medication status
* tune best classifier for 4 group model
* go back to tune best classifier for clinical groups