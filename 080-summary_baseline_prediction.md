# 2020-02-18 14:48:27

Time to summarize all results so far. We start with a sample of 393 subjects
that have PRS. Some details about the sample:

```r
> data0 = readRDS('~/data/baseline_prediction/prs_start/complete_massagedResids_02052020.rds')
> table(data0$sex)

Female   Male 
   130    263 
> mean(data0$base_age)
[1] 8.238728
> sd(data0$base_age)
[1] 2.64162
> mean(data0$last_age)
[1] 13.08527
> sd(data0$last_age)
[1] 2.979007
> table(data0$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       159        115         73         46 
> table(data0$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       159         76        112         46 
```

There are 6 data domains:

* neuropsych: 'FSIQ', "VMI.beery", "SSB.wisc", "SSF.wisc", 'DSF.wisc',
  'DSB.wisc', "DS.wj", "VM.wj"
* demographic: 'base_age', 'sex', 'SES'
* genomics: PRS scores for the entire cohort (i.e. not the European-only PRS)
* DTI: FA values for JHU tracts, collapsed to reduce variables
* anatomy: thickness data for collapsed Freesurfer lobar regions
* clinics: 'internalizing', 'externalizing', base_sx

## Analysis 1: univariate results

Within each data domain, check which variables are significantly predicted in
the model that uses the linear relationship between the 4 groups as the main
predictor. In other words:

```
model = lme(myvar ~ ordered + covariates, ~1|FAMID)
```

And we collected the p-value and betas for the linear model associated with the
variable **ordered**. The order of the groups is always ('nv012',
'notGE6adhd', 'imp', 'nonimp'). The covariates change per data domain:

* neuropsych: base_age
* demographic: sex
* genomics: population PCs, base_age, sex
* DTI: "meanX.trans", "meanY.trans", "meanZ.trans", "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes", age_at_scan, sex
* anatomy: "mprage_score", "ext_avg", "int_avg", age_at_scan, sex
* clinics: sex

Those were the initial covariates in each domain, but I used stepAIC to select
which covariates we should keep for the best model possible.

Also, the number of observations varies per domain, as there was no imputation
in this analysis. So, the final number per domain is as follows:

```
> idx=!is.na(data0[, 'FSIQ'])
> sum(idx)
[1] 390
> table(data0[idx,]$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       157        114         73         46 
> table(data0[idx,]$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       157         75        112         46 
> idx=!is.na(data0[, 'VMI.beery'])
> sum(idx)
[1] 315
> table(data0[idx,]$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       124         91         61         39 
> table(data0[idx,]$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       124         60         92         39 
> idx=!is.na(data0[, 'SSB.wisc'])
> sum(idx)
[1] 245
> table(data0[idx,]$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
        91         74         42         38 
> table(data0[idx,]$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
        91         49         67         38 
> idx=!is.na(data0[, 'DS.wj'])
> sum(idx)
[1] 339
> table(data0[idx,]$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       137        102         57         43 
> table(data0[idx,]$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       137         66         93         43 
> idx=!is.na(data0[, 'CC_fa'])
> sum(idx)
[1] 180
> table(data0[idx,]$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
        74         48         31         27 
> table(data0[idx,]$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
        74         32         47         27 
> idx=!is.na(data0[, 'parietal'])
> sum(idx)
[1] 286
> table(data0[idx,]$thresh0.00_inatt_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       124         78         47         37 
> table(data0[idx,]$thresh0.50_hi_GE6_wp05)

     nv012        imp     nonimp notGE6adhd 
       124         52         73         37 
```

Finally, the data used for anat and DTI are already the result of excluding
outliers using the same methods we used in the heritability of change paper,
after removing any scans outside the 95th percentile. Specifically, all longitudinal
scans for the 393 subjects with PRS were taken into consideration. After removing the
outliers, I picked the first scan and only kept it if it was acquired within 1
year of the baseline clinical assessment.

### Results

So, the complete list of variables tested was:

```
VMI.beery
VM.wj
FSIQ
externalizing
IFO_fa
DS.wj
ADHD_PRS0.001000 (only the best PRS was taken into FDR)
OFC
ATR_fa
CST_fa
base_age
cingulate
DSF.wisc
UNC_fa
SSB.wisc
temporal
CC_fa
SSF.wisc
parietal
ILF_fa
CIN_fa
frontal
insula
sensorimotor
DSB.wisc
internalizing
occipital
SES
SLF_fa
```

Our chosen approach was to do an FDR across all domains. Then it's a matter of
choosing what q level to go with:

```
> hi_p[p2<.05,'var']
[1] VMI.beery        VM.wj            FSIQ             externalizing   
[5] IFO_fa           DS.wj            ADHD_PRS0.001000 OFC             
[9] ATR_fa          
> hi_p[p2<.1,'var']
 [1] VMI.beery        VM.wj            FSIQ             externalizing   
 [5] IFO_fa           DS.wj            ADHD_PRS0.001000 OFC             
 [9] ATR_fa           CST_fa           base_age         cingulate       
[13] DSF.wisc        
> inatt_p[p2<.05,'var']
[1] FSIQ             VMI.beery        VM.wj            base_age        
[5] externalizing    ADHD_PRS0.000500
> inatt_p[p2<.1,'var']
[1] FSIQ             VMI.beery        VM.wj            base_age        
[5] externalizing    ADHD_PRS0.000500 DSF.wisc         IFO_fa          
[9] DS.wj           
```

Or we could choose to not include base_age, because it's not really a predictor
and more of a covariate anyways (like sex):

```
> hi_p[p2<.05,'var']
[1] VMI.beery        VM.wj            FSIQ             externalizing   
[5] IFO_fa           DS.wj            ADHD_PRS0.001000 OFC             
[9] ATR_fa          
> hi_p[p2<.1,'var']
 [1] VMI.beery        VM.wj            FSIQ             externalizing   
 [5] IFO_fa           DS.wj            ADHD_PRS0.001000 OFC             
 [9] ATR_fa           CST_fa           cingulate        DSF.wisc        
> inatt_p[p2<.05,'var']
[1] FSIQ             VMI.beery        VM.wj            externalizing   
[5] ADHD_PRS0.000500
> inatt_p[p2<.1,'var']
[1] FSIQ             VMI.beery        VM.wj            externalizing   
[5] ADHD_PRS0.000500 DSF.wisc         IFO_fa           DS.wj           
```

## Analysis 2: big model

This analysis had the goal of checking how well we can model the groups by
combining all the "good" univariate variables from analysis 1.

As we cannot deal with NAs here, we decided to impute the data using as the base
the 180 kids who have both PRS and DTI. 

The model is a multinomial logistic regression that ignores the family term. It
is also not ordered, because it performed better (i.e. higher ROC AUC) than the
ordered model:

```
group ~ good_vars + covars
```

The good_vars came from analysis 1:

```r
hi_vars = c('VMI.beery', 'VM.wj', 'FSIQ', 'externalizing', 'IFO_fa', 'DS.wj',
            'ADHD_PRS0.001000', 'OFC', 'ATR_fa', 'CST_fa', 'cingulate',
            'DSF.wisc')
inatt_vars = c('FSIQ', 'VMI.beery', 'VM.wj', 'externalizing',
               'ADHD_PRS0.000500', 'DSF.wisc', 'IFO_fa', 'DS.wj')
covars = c('base_age', 'sex')
```

Also, note that every domain has already been residualized within domain (for
example, PRS was residualized based on PCs, etc). We can evaluate the models
based on ROC AUC, and check how important each variable was in the prediction
(just the sum of the absolute value of the coefficients across the different categories). 

The results are quite disproportional, as one would expect since externalizing
is very unbalanced:

```
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
Multi-class area under the curve: 0.7411

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
Multi-class area under the curve: 0.7578
```

How does it look when we add in base_sx as one of the predictors?

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
Multi-class area under the curve: 0.9138

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
Multi-class area under the curve: 0.9169
```

Big jump in AUC, which makes sense as base_sx is used in defining the groups. But
interesting that externalizing is almost as good, maybe even better...

We can also check the increase in AUC as we add the different domains in a
step-wise fashion:

```
[1] "inatt"
      demo       clin        gen        dti neuropsych 
 0.6134810  0.8804639  0.8919944  0.8922737  0.9137963 
[1] "hi"
      demo       clin        gen        dti       anat neuropsych 
 0.6065227  0.8345617  0.8455017  0.8616758  0.8765961  0.9169098 
```

Out of curiosity, what are our numbers if we use the exact same framework, but
trying to distinguish just between the two clinical groups? Results are also not
that bad...

```
[1] "inatt"
      demo       clin        gen        dti neuropsych 
 0.6871640  0.8185484  0.8319892  0.8319892  0.8830645 
[1] "hi"
      demo       clin        gen        dti       anat neuropsych 
 0.6702128  0.8071809  0.8091755  0.8138298  0.8244681  0.8590426 
```

Now, of course the externalizing distribution is very sparse:

```
> table(data$externalizing, data[, phen])

    nv012 imp nonimp notGE6adhd
  0    74  27     42         25
  1     0   5      5          2
```

So, what happens to the weights if I either don't run that variable in the 4
group comparison? Or don't include nv012, and only run 3 groups?

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
Multi-class area under the curve: 0.747
```

Overall, not as bad, and we didn't lose much in terms of AUC. Let's see how it
looks after we add in base_sx: 

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
Multi-class area under the curve: 0.9163
```

We still see the big jump from before, again expected as base_sx is used to
define the groups. What if we run 3 groups only?

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
Multi-class area under the curve: 0.7823
```

It still carries weight, but at least it's not over the top. AUC is still
decent... And if we include base_sx?

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
Multi-class area under the curve: 0.8786
```

Yep, make sense. Lastly, what happens when we add in
medication_status_at_observation? First, the distribution is very skewed, like
externalizing:

```
> table(data$medication_status_at_observation, data$thresh0.50_hi_GE6_wp05)
         
          nv012 imp nonimp notGE6adhd
  none      159  52     78         37
  nonstim     0   3      2          1
  stim        0  21     32          8
```

So, of course the linear fit is extremely significant. 




## Analysis 3: ML


