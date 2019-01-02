# 2019-01-02 14:20:44

Let's take a look at how a simple model performs in predicting OLS, or odd
ratio, using the combined datasets after the descriptives selection. 

Let's start with trying to predict the actual OLS in each modality, trying to
figure out how much better we do by adding different modalities.

Another thing to consider is that we could use mixed models or not... something
to think about.

## OLS-inatt

```r
load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz')
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
idx = df$diag_group != 'new_onset' & df$DX != 'NV'
struct = df[!is.na(df$HI_vol_rh) & idx,]
fit = lm(OLS_inatt_slope ~ inatt_vol_lh, data=struct)
```

Well, it's not good that we're not doing better than using the mean:

```
> summary(fit)

Call:
lm(formula = OLS_inatt_slope ~ inatt_vol_lh, data = struct)

Residuals:
     Min       1Q   Median       3Q      Max 
-3.05975 -0.32130 -0.05971  0.38179  1.87399 

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept)   1.294836   0.267107   4.848 3.28e-06 ***
inatt_vol_lh -0.016050   0.002952  -5.438 2.34e-07 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.6839 on 140 degrees of freedom
Multiple R-squared:  0.1744,    Adjusted R-squared:  0.1685 
F-statistic: 29.57 on 1 and 140 DF,  p-value: 2.337e-07

> p = predict(fit)
> mean(abs(struct$OLS_inatt_slope - p))
[1] 0.4914881
> m = mean(struct$OLS_inatt_slope)
> mean(abs(struct$OLS_inatt_slope - m))
[1] 0.5147604
```

We're only doing slightly better than predicting with the mean. Not good. What
if we use the actual model we used in finding the voxels?

```r
library(nlme)
library(gdata
)
load('/data/NCR_SBRB/baseline_prediction/struct_volume_11142018_260timeDiff12mo.RData.gz')
struct = merge(struct, data, by='MRN')
mprage = read.xls('/data/NCR_SBRB/baseline_prediction/long_scans_08072018.xlsx',
                  sheet='mprage')
struct = merge(struct, mprage, by.x='mask.id', by.y='Mask.ID...Scan')
qc = read.csv('/data/NCR_SBRB/baseline_prediction/master_qc.csv')
struct = merge(struct, qc, by.x='mask.id', by.y='Mask.ID')
fm = as.formula("OLS_inatt_slope ~ inatt_vol_lh + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC + age_at_scan + I(age_at_scan^2)")
fit2 = lme(fm, random=~1|nuclearFamID, data=struct)
```

```
> summary(fit2)
Linear mixed-effects model fit by REML
 Data: struct 
       AIC      BIC    logLik
  338.2944 367.2728 -159.1472

Random effects:
 Formula: ~1 | nuclearFamID
        (Intercept)  Residual
StdDev:   0.1754667 0.6543485

Fixed effects: list(fm) 
                           Value Std.Error  DF   t-value p-value
(Intercept)            1.5888971 1.0566770 128  1.503673  0.1351
inatt_vol_lh          -0.0157901 0.0030138   6 -5.239252  0.0019
Sex...SubjectsMale     0.1755870 0.1314736   6  1.335530  0.2301
ext_avg_freesurfer5.3  0.2061340 0.1785680   6  1.154372  0.2922
int_avg_freesurfer5.3  0.2328711 0.1886450   6  1.234441  0.2632
mprage_QC             -0.1980678 0.1393510   6 -1.421360  0.2050
age_at_scan           -0.1853818 0.1991918   6 -0.930670  0.3879
I(age_at_scan^2)       0.0077219 0.0104256   6  0.740669  0.4869
 Correlation: 
                      (Intr) intt__ S...SM e__5.3 i__5.3 mpr_QC ag_t_s
inatt_vol_lh          -0.224                                          
Sex...SubjectsMale    -0.097 -0.162                                   
ext_avg_freesurfer5.3 -0.252  0.072 -0.049                            
int_avg_freesurfer5.3 -0.224  0.102  0.015 -0.206                     
mprage_QC             -0.144  0.005  0.048 -0.051 -0.318              
age_at_scan           -0.853 -0.072  0.045  0.019 -0.029  0.055       
I(age_at_scan^2)       0.828  0.058 -0.034 -0.042  0.036 -0.021 -0.990

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-4.32674098 -0.39250228 -0.02872832  0.55649900  2.50200844 

Number of Observations: 142
Number of Groups: 129 

> p2 = predict(fit2)
> mean(abs(struct$OLS_inatt_slope - p2))
[1] 0.4406693
> summary(struct$OLS_inatt_slope)
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-4.00000 -0.43088 -0.06869 -0.12368  0.22419  1.91122 
```

Yep, this is a bit better, but still somewhat high and too close to the mean
prediction. I could also winsorize it and see what the damage is.

```r
winsorize = function(x, cut = 0.01){
  cut_point_top <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x, cut, na.rm = T)
  i = which(x >= cut_point_top) 
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom) 
  x[j] = cut_point_bottom
  return(x)
}
struct$wInatt = winsorize(struct$OLS_inatt_slope)
fm = as.formula("wInatt ~ inatt_vol_lh + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC + age_at_scan + I(age_at_scan^2)")
fit3 = lme(fm, random=~1|nuclearFamID, data=struct)
```

```
> p3 = predict(fit3)
> mean(abs(struct$wInatt - p3))
[1] 0.4527394
> m = mean(struct$wInatt)
>  mean(abs(struct$wInatt - m))
[1] 0.4983164
```

The results are actually better without winsorizing...

I could try this for all different brain regions and their combinations, but
let's try to dichotomize the groups first.

## creating groups

Let's use anyone who has at least one imaging data point:

### inatt

```r
load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz')
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
idx = df$diag_group != 'new_onset' & df$DX != 'NV'
idx2 = !is.na(df$inatt_vol_lh) | !is.na(df$inatt_AD_clu1) | !is.na(df$inatt_melodic_DMN)
imaging = df[idx & idx2,]
idx = df$diag_group != 'new_onset' & df$DX != 'NV'
idx2 = !is.na(df$HI_vol_rh) | !is.na(df$HI_RD_clu1) | !is.na(df$inatt_melodic_DMN)
imaging = df[idx & idx2,]
```

So, first thing is to do a quick check on how the latent classes groups Philip
had generated split, when only using non-NVs and people without new onset, who
had at least one imaging modality. 

```
> table(imaging$group_HI3)

 1  2  3 
65 33 64 
> table(imaging$group_INATT3)

 1  2  3 
71 71 20 
```

So, that would be one way to do it. But considering that those groups took into
account NVs, it might be better to not use them. So, we could also explore
different cutoffs in the OLS ranges. If we go for 3 equally split groups:

```
> quantile(imaging$OLS_inatt_slope, c(.33, .66))
        33%         66% 
-0.37288520  0.05786572 
> quantile(imaging$OLS_inatt_slope)
         0%         25%         50%         75%        100% 
-4.00000000 -0.46406224 -0.08101091  0.19347778  1.91122072 
```

So, in innattention it makes sense for improvement to be < .37, deterioration to
be > .06, and anything in the middle means stable. Don't see why the cutoffs
would need to be symmetric, other than for easy of understanding.

For HI it gets a bit more complicated:

```
> quantile(imaging$OLS_HI_slope)
        0%        25%        50%        75%       100% 
-3.4313725 -0.9321712 -0.5134901  0.0000000  4.0000000 
> quantile(imaging$OLS_HI_slope, c(.33, .66))
       33%        66% 
-0.7769565 -0.1767975 
```

In other words,
> dim(imaging)
[1] 162  44
> sum(imaging$OLS_inatt_slope>=0)
[1] 68
> sum(imaging$OLS_HI_slope>=0)
[1] 42
> 

```