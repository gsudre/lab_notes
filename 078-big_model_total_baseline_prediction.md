# 2020-02-06 11:30:50

Philip asked me to run the big model for total symptom count too. Let's just
repeat the code from 75:

```r
library(caret)
data0 = readRDS('~/data/baseline_prediction/prs_start/complete_massagedResids_02052020.rds')
data = data0[!is.na(data0$IFO_fa), ]
data$externalizing = as.factor(data$externalizing)

set.seed(42)
base_vars = c(colnames(data)[42:65], colnames(data)[74:81])
# anatomical
imp_vars = colnames(data)[66:73]
test = preProcess(data[, c(base_vars, imp_vars)], method = "bagImpute")
data[, c(base_vars, imp_vars)] <- predict(test, data[, c(base_vars, imp_vars)])
# beery
imp_vars = colnames(data)[89]
test = preProcess(data[, c(base_vars, imp_vars)], method = "bagImpute")
data[, c(base_vars, imp_vars)] <- predict(test, data[, c(base_vars, imp_vars)])
# wj
imp_vars = colnames(data)[96:98]
test = preProcess(data[, c(base_vars, imp_vars)], method = "bagImpute")
data[, c(base_vars, imp_vars)] <- predict(test, data[, c(base_vars, imp_vars)])
# wisc
imp_vars = colnames(data)[90:95]
test = preProcess(data[, c(base_vars, imp_vars)], method = "bagImpute")
data[, c(base_vars, imp_vars)] <- predict(test, data[, c(base_vars, imp_vars)])
```

Now, let's see how the overall AUC goes up or down the mode we add the different
domains. But before we do that, let's switch to using the actual function,
because this way we can specify maxit:

```r
library(nnet)
library(pROC)
hi_vars = c('VMI.beery', 'VM.wj', 'FSIQ', 'externalizing', 'IFO_fa', 'DS.wj',
            'ADHD_PRS0.001000', 'OFC', 'ATR_fa', 'CST_fa', 'cingulate',
            'DSF.wisc')
inatt_vars = c('FSIQ', 'VMI.beery', 'VM.wj', 'externalizing',
               'ADHD_PRS0.000500', 'DSF.wisc', 'IFO_fa', 'DS.wj')
covars = c('base_age', 'sex')
min_sx = 6

for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('this_data = data[, c(phen, %s_vars, covars)]', sx)))

    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])

    eval(parse(text=sprintf('predictors_str=paste(%s_vars, collapse="+")', sx)))
    fm_str = paste(phen, " ~ ", predictors_str, ' + ',
               paste(covars, collapse='+'),
               sep="")
    fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
    preds = predict(fit, type='prob')
    print(sx)
    print(varImp(fit))
    print(multiclass.roc(this_data[, phen], preds))
}
```

The results are quite disproportional:

```
# weights:  48 (33 variable)
initial  value 249.532985 
iter  10 value 200.388194
iter  20 value 194.484596
iter  30 value 193.542152
iter  40 value 193.278960
iter  50 value 193.253013
final  value 193.252983 
converged
[1] "inatt"
                    Overall
FSIQ              0.7862872
VMI.beery         1.6524547
VM.wj             1.1801762
externalizing1   49.3133069
ADHD_PRS0.000500  0.8268296
DSF.wisc          0.4345088
IFO_fa            0.9686923
DS.wj             1.3378331
base_age          1.7978418
sexMale           1.0346242

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.7411
# weights:  64 (45 variable)
initial  value 249.532985 
iter  10 value 198.713302
iter  20 value 192.503757
iter  30 value 191.683566
iter  40 value 191.339960
iter  50 value 191.262836
iter  60 value 191.257237
final  value 191.257157 
converged
[1] "hi"
                    Overall
VMI.beery         1.7175702
VM.wj             0.7886613
FSIQ              0.9107406
externalizing1   47.7411726
IFO_fa            0.5858586
DS.wj             1.5779318
ADHD_PRS0.001000  0.6740715
OFC               0.9341405
ATR_fa            0.4345321
CST_fa            0.7608150
cingulate         0.2842507
DSF.wisc          0.5487167
base_age          1.2516826
sexMale           1.2585796

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.7578
```

How does it look when we add base_sx?

```r
for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('this_data = data[, c(phen, %s_vars, covars, "base_%s")]',
                            sx, sx)))

    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])

    eval(parse(text=sprintf('predictors_str=paste(%s_vars, collapse="+")', sx)))
    fm_str = paste(phen, " ~ ", predictors_str, sprintf(' + base_%s + ', sx),
               paste(covars, collapse='+'),
               sep="")
    fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
    preds = predict(fit, type='prob')
    print(sx)
    print(varImp(fit))
    print(multiclass.roc(this_data[, phen], preds))
}
```

```
[1] "inatt"
                    Overall
FSIQ              0.6806480
VMI.beery         1.3805707
VM.wj             1.2073122
externalizing1   47.3177951
ADHD_PRS0.000500  1.2894672
DSF.wisc          2.7361367
IFO_fa            1.5614350
DS.wj             0.6733051
base_inatt       27.3284253
base_age          0.9617274
sexMale           2.6653584

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.9138
# weights:  68 (48 variable)
initial  value 249.532985 
iter  10 value 114.169643
iter  20 value 105.005824
iter  30 value 101.410331
iter  40 value 101.225247
iter  50 value 101.224064
iter  60 value 101.221782
final  value 101.221646 
converged
[1] "hi"
                    Overall
VMI.beery         1.8009884
VM.wj             3.2346574
FSIQ              1.3900451
externalizing1   23.7156506
IFO_fa            1.6400383
DS.wj             0.4748365
ADHD_PRS0.001000  2.0801047
OFC               1.5064806
ATR_fa            0.5491814
CST_fa            2.7601051
cingulate         1.1449211
DSF.wisc          0.5819857
base_hi          24.2712899
base_age          3.8420027
sexMale           2.6573462

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.9169
```

Big jump, which makes sense as base_sx is used in defining the groups. But
interesting that externalizing is almost as good, maybe even better...

Now, let's go systematically to see how much improvement we get by adding each
new domain...

## Stepwise domain addition

```r
hi_vars = list(demo = c('base_age', 'sex'),
               clin = c('base_hi', 'externalizing'),
               gen = c('ADHD_PRS0.001000'),
               dti = c('ATR_fa', 'CST_fa', 'IFO_fa'),
               anat = c('OFC', 'cingulate'),
               neuropsych = c('VMI.beery', 'VM.wj', 'FSIQ', 'DS.wj', 'DSF.wisc'))
inatt_vars = list(demo = c('base_age', 'sex'),
                  clin = c('base_inatt', 'externalizing'),
                  gen = c('ADHD_PRS0.000500'),
                  dti = c('IFO_fa'),
                  neuropsych = c('FSIQ', 'VMI.beery', 'VM.wj', 'DSF.wisc', 'DS.wj'))
min_sx = 6

for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('my_vars = %s_vars', sx)))
    cur_vars = c()
    scores = c()
    for (dom in 1:length(my_vars)) {
        cur_vars = c(cur_vars, my_vars[[dom]])
        this_data = data[, c(phen, cur_vars)]
        scale_me = c()
        for (v in colnames(this_data)) {
            if (!is.factor(this_data[, v])) {
                scale_me = c(scale_me, v)
            }
        }
        this_data[, scale_me] = scale(this_data[, scale_me])

        predictors_str = paste(cur_vars, collapse="+")
        fm_str = paste(phen, " ~ ", predictors_str, sep="")
        fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
        preds = predict(fit, type='prob')
        scores = c(scores, multiclass.roc(this_data[, phen], preds)$auc)
    }
    print(sx)
    names(scores) = names(my_vars)
    print(scores)
}
```

Here's the stepwise progression:

```
[1] "inatt"
      demo       clin        gen        dti neuropsych 
 0.6134810  0.8804639  0.8919944  0.8922737  0.9137963 
[1] "hi"
      demo       clin        gen        dti       anat neuropsych 
 0.6065227  0.8345617  0.8455017  0.8616758  0.8765961  0.9169098 
```

Out of curiosity, what are our numbers if we use the exact same framework, but
trying to distinguish just between the two clinical groups?

```r
for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('my_vars = %s_vars', sx)))
    cur_vars = c()
    scores = c()
    for (dom in 1:length(my_vars)) {
        cur_vars = c(cur_vars, my_vars[[dom]])
        this_data = data[, c(phen, cur_vars)]
        this_data = this_data[this_data[, phen] == 'nonimp' | this_data[, phen] == 'imp',]
        this_data[, phen] = factor(this_data[, phen], ordered=F)
        this_data[, phen] = relevel(this_data[, phen], ref='imp')

        scale_me = c()
        for (v in colnames(this_data)) {
            if (!is.factor(this_data[, v])) {
                scale_me = c(scale_me, v)
            }
        }
        this_data[, scale_me] = scale(this_data[, scale_me])

        predictors_str = paste(cur_vars, collapse="+")
        fm_str = paste(phen, " ~ ", predictors_str, sep="")
        fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
        preds = predict(fit, type='prob')
        scores = c(scores, multiclass.roc(this_data[, phen], preds)$auc)
    }
    print(sx)
    names(scores) = names(my_vars)
    print(scores)
}
```

Results are also not bad:

```
[1] "inatt"
      demo       clin        gen        dti neuropsych 
 0.6871640  0.8185484  0.8319892  0.8319892  0.8830645 
[1] "hi"
      demo       clin        gen        dti       anat neuropsych 
 0.6702128  0.8071809  0.8091755  0.8138298  0.8244681  0.8590426 
```

Now, the interesting thing is that the externalizing distribution is very
sparse:

```
> table(data$externalizing, data[, phen])

    nv012 imp nonimp notGE6adhd
  0    74  27     42         25
  1     0   5      5          2
```

So, what happens to the weights if I either don't run that variable, or don't
include nv012?

```r
library(nnet)
library(pROC)
hi_vars = c('VMI.beery', 'VM.wj', 'FSIQ', 'IFO_fa', 'DS.wj',
            'ADHD_PRS0.001000', 'OFC', 'ATR_fa', 'CST_fa', 'cingulate',
            'DSF.wisc')
inatt_vars = c('FSIQ', 'VMI.beery', 'VM.wj',
               'ADHD_PRS0.000500', 'DSF.wisc', 'IFO_fa', 'DS.wj')
covars = c('base_age', 'sex')
min_sx = 6

for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('this_data = data[, c(phen, %s_vars, covars)]', sx)))

    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])

    eval(parse(text=sprintf('predictors_str=paste(%s_vars, collapse="+")', sx)))
    fm_str = paste(phen, " ~ ", predictors_str, ' + ',
               paste(covars, collapse='+'),
               sep="")
    fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
    preds = predict(fit, type='prob')
    print(sx)
    print(varImp(fit))
    print(multiclass.roc(this_data[, phen], preds))
}
```

This is not as bad, and we didn't lose much in terms of AUC:

```
[1] "inatt"
                   Overall
FSIQ             0.6713690
VMI.beery        1.8008210
VM.wj            1.1684526
ADHD_PRS0.000500 0.7708455
DSF.wisc         0.4803280
IFO_fa           1.1708453
DS.wj            1.0548323
base_age         1.6850936
sexMale          0.9483320

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.7229

[1] "hi"
                   Overall
VMI.beery        1.8268554
VM.wj            0.8664441
FSIQ             0.9230327
IFO_fa           0.7473989
DS.wj            1.2511845
ADHD_PRS0.001000 0.6451152
OFC              0.8566875
ATR_fa           0.5427812
CST_fa           0.7208566
cingulate        0.4296081
DSF.wisc         0.6122018
base_age         1.1866683
sexMale          1.1602818

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.747
```

Let's see how it looks after we add in base_sx:

```r
for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('this_data = data[, c(phen, %s_vars, covars, "base_%s")]',
                            sx, sx)))

    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])

    eval(parse(text=sprintf('predictors_str=paste(%s_vars, collapse="+")', sx)))
    fm_str = paste(phen, " ~ ", predictors_str, sprintf(' + base_%s + ', sx),
               paste(covars, collapse='+'),
               sep="")
    fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
    preds = predict(fit, type='prob')
    print(sx)
    print(varImp(fit))
    print(multiclass.roc(this_data[, phen], preds))
}
```

Big jump as expected:

```
[1] "inatt"
                    Overall
FSIQ              0.5190794
VMI.beery         1.2832258
VM.wj             1.1350190
ADHD_PRS0.000500  1.1784921
DSF.wisc          3.0009173
IFO_fa            1.6415447
DS.wj             0.7241186
base_inatt       28.0924540
base_age          0.8509839
sexMale           2.2257656

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.9115

[1] "hi"
                    Overall
VMI.beery         1.8090450
VM.wj             3.2007661
FSIQ              1.4060350
IFO_fa            1.6473047
DS.wj             0.5276963
ADHD_PRS0.001000  2.0131734
OFC               1.5049995
ATR_fa            0.5186103
CST_fa            2.7554278
cingulate         1.1505311
DSF.wisc          0.5744980
base_hi          24.1625875
base_age          3.8087047
sexMale           2.6791743

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 4 levels of this_data[, phen]: nv012, imp, nonimp, notGE6adhd.
Multi-class area under the curve: 0.9163
```

What if we reduce it to a 3 group comparison only?

```r
hi_vars = c('VMI.beery', 'VM.wj', 'FSIQ', 'externalizing', 'IFO_fa', 'DS.wj',
            'ADHD_PRS0.001000', 'OFC', 'ATR_fa', 'CST_fa', 'cingulate',
            'DSF.wisc')
inatt_vars = c('FSIQ', 'VMI.beery', 'VM.wj', 'externalizing',
               'ADHD_PRS0.000500', 'DSF.wisc', 'IFO_fa', 'DS.wj')
covars = c('base_age', 'sex')
for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('this_data = data[, c(phen, %s_vars, covars)]', sx)))

    this_data = this_data[this_data[, phen] != 'nv012',]
    this_data[, phen] = factor(this_data[, phen], ordered=F)
    this_data[, phen] = relevel(this_data[, phen], ref='notGE6adhd')

    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])

    eval(parse(text=sprintf('predictors_str=paste(%s_vars, collapse="+")', sx)))
    fm_str = paste(phen, " ~ ", predictors_str, ' + ',
               paste(covars, collapse='+'),
               sep="")
    fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
    preds = predict(fit, type='prob')
    print(sx)
    print(varImp(fit))
    print(multiclass.roc(this_data[, phen], preds))
}
```

It still carries weight, but at least it's not over the top.

```
[1] "inatt"
                   Overall
FSIQ             1.3936495
VMI.beery        0.4464261
VM.wj            0.7960408
externalizing1   2.4373028
ADHD_PRS0.000500 1.3741709
DSF.wisc         0.5311392
IFO_fa           0.4959420
DS.wj            0.4958580
base_age         1.0964769
sexMale          1.1486840

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 3 levels of this_data[, phen]: notGE6adhd, imp, nonimp.
Multi-class area under the curve: 0.7569

[1] "hi"
                   Overall
VMI.beery        0.4864423
VM.wj            0.5471967
FSIQ             1.6412003
externalizing1   2.5931840
IFO_fa           0.6310506
DS.wj            0.4103319
ADHD_PRS0.001000 1.1681734
OFC              1.0591294
ATR_fa           0.2912992
CST_fa           0.9043042
cingulate        0.2246797
DSF.wisc         0.5846238
base_age         1.0156251
sexMale          1.6160933

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 3 levels of this_data[, phen]: notGE6adhd, imp, nonimp.
Multi-class area under the curve: 0.7823
```

And if we include base_sx...

```r
for (sx in c('inatt', 'hi')) {
    set.seed(42)
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    }
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    eval(parse(text=sprintf('this_data = data[, c(phen, %s_vars, covars, "base_%s")]',
                            sx, sx)))

    this_data = this_data[this_data[, phen] != 'nv012',]
    this_data[, phen] = factor(this_data[, phen], ordered=F)
    this_data[, phen] = relevel(this_data[, phen], ref='notGE6adhd')

    scale_me = c()
    for (v in colnames(this_data)) {
        if (!is.factor(this_data[, v])) {
            scale_me = c(scale_me, v)
        }
    }
    this_data[, scale_me] = scale(this_data[, scale_me])

    eval(parse(text=sprintf('predictors_str=paste(%s_vars, collapse="+")', sx)))
    fm_str = paste(phen, " ~ ", predictors_str, sprintf(' + base_%s + ', sx),
               paste(covars, collapse='+'),
               sep="")
    fit = multinom(as.formula(fm_str), data=this_data, maxit=2000)
    preds = predict(fit, type='prob')
    print(sx)
    print(varImp(fit))
    print(multiclass.roc(this_data[, phen], preds))
}
```

Yep, makes sense:

```
[1] "inatt"
                   Overall
FSIQ             0.7149091
VMI.beery        0.2881959
VM.wj            1.1430780
externalizing1   2.7247453
ADHD_PRS0.000500 1.1872987
DSF.wisc         0.1431103
IFO_fa           0.8085379
DS.wj            0.5387884
base_inatt       5.3156872
base_age         0.7218521
sexMale          2.4801733

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 3 levels of this_data[, phen]: notGE6adhd, imp, nonimp.
Multi-class area under the curve: 0.8596

[1] "hi"
                   Overall
VMI.beery        0.7125827
VM.wj            0.6703544
FSIQ             1.7043603
externalizing1   1.1160982
IFO_fa           1.0440388
DS.wj            0.3050194
ADHD_PRS0.001000 1.9647090
OFC              1.2091195
ATR_fa           0.5853600
CST_fa           1.3282894
cingulate        0.2831388
DSF.wisc         1.2255681
base_hi          4.9314779
base_age         3.7619075
sexMale          1.8857750

Call:
multiclass.roc.default(response = this_data[, phen], predictor = preds)

Data: multivariate predictor preds with 3 levels of this_data[, phen]: notGE6adhd, imp, nonimp.
Multi-class area under the curve: 0.8786
```


# TODO
