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
library(gdata)
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

And we do have 162 people in this pool (SX >= 3 and no new_onset).

# 2019-01-07 09:57:28

Philip suggested I should go with inatt marked improvent at 0.33 (1 symtom every
three years). Deterioration anything over zero. And mild improvement 0 to 0.33.
For HI, I could have marked improvement at 0.5 (1 symtom every 2 years). Mild at
zero to 0.5, and deterioration anything above zero.

Let's see how many participants I get using those cut-offs:

```
> sum(imaging$OLS_HI_slope <= -.5)
[1] 83
> sum(imaging$OLS_HI_slope > -.5 & imaging$OLS_HI_slope <= 0)
[1] 46
> sum(imaging$OLS_HI_slope > 0)
[1] 33
> sum(imaging$OLS_inatt_slope <= -.33)
[1] 60
> sum(imaging$OLS_inatt_slope > -.33 & imaging$OLS_inatt_slope <= 0)
[1] 39
> sum(imaging$OLS_inatt_slope > 0)
[1] 63
```

The splits in HI are not great, and will also change based on different imaging
modalities, but we have to go with what we have.

## struct inattention

```
> imaging[imaging$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
> imaging[imaging$OLS_inatt_slope > -.33 & imaging$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
> imaging[imaging$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
> imaging$OLS_inatt_categ = as.factor(imaging$OLS_inatt_categ)
> table(imaging$OLS_inatt_categ)
 deter marked   mild 
    63     60     39 
```

So, a few things I learned:

* (of course) scaling the independent variable doesn't change the model quality, only the
  interpretation of the coefficients
* changing the reference level for the factor doesn;t change the relationships
  or quality of the model; just the information that gets output.
* like in a lmfit, adding more independent variables changes how goodness of the
  model, but the coefficients and their significance of other variables doesn't
  change as much (of course, it depends on how correlated the new variables are
  to old ones). But a good approach could be starting with all the variables we
  used in selecting the voxels, and then trim it down to just the significant
  bits for the final model.
* it might be better to scale the brain variables because the model
  interpretation is based on a one-unit increase. One unit in a scaled variable
  equals one standard deviation, which is more meaningful than one FA, for
  example. 

```r
library(gdata)
library(nnet)

load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz')
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
idx = df$diag_group != 'new_onset' & df$DX != 'NV'
idx2 = !is.na(df$inatt_vol_lh) | !is.na(df$inatt_AD_clu1) | !is.na(df$inatt_melodic_DMN)
imaging = df[idx & idx2,]
idx = df$diag_group != 'new_onset' & df$DX != 'NV'
idx2 = !is.na(df$HI_vol_rh) | !is.na(df$HI_RD_clu1) | !is.na(df$inatt_melodic_DMN)
imaging = df[idx & idx2,]
imaging[imaging$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
imaging[imaging$OLS_inatt_slope > -.33 & imaging$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
imaging[imaging$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
imaging$OLS_inatt_categ = as.factor(imaging$OLS_inatt_categ)
imaging$OLS_inatt_categ = relevel(imaging$OLS_inatt_categ, ref='mild')

load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz')
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
idx = df$diag_group != 'new_onset' & df$DX != 'NV'
struct = df[!is.na(df$HI_vol_rh) & idx,]
load('/data/NCR_SBRB/baseline_prediction/struct_volume_11142018_260timeDiff12mo.RData.gz')
struct = merge(struct, data, by='MRN') # put mask ids in combined dataset
mprage = read.xls('/data/NCR_SBRB/baseline_prediction/long_scans_08072018.xlsx',
                  sheet='mprage')
struct = merge(struct, mprage, by.x='mask.id', by.y='Mask.ID...Scan') # get demographics
qc = read.csv('/data/NCR_SBRB/baseline_prediction/master_qc.csv')
struct = merge(struct, qc, by.x='mask.id', by.y='Mask.ID') # get QC scores
df = merge(struct, imaging, by='MRN')
fit1 <- multinom(OLS_inatt_categ ~ scale(inatt_vol_lh.x) + age_at_scan + I(age_at_scan^2) + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC, data = df, na.action=na.omit)
z1 <- summary(fit1)$coefficients/summary(fit1)$standard.errors
p1 <- (1 - pnorm(abs(z1), 0, 1)) * 2
rr1 = exp(coef(fit1))
pp1 = fitted(fit1)
fit2 <- multinom(OLS_inatt_categ ~ scale(inatt_vol_lh.x), data = df, na.action=na.omit)
z2 <- summary(fit2)$coefficients/summary(fit2)$standard.errors
p2 <- (1 - pnorm(abs(z2), 0, 1)) * 2
rr2 = exp(coef(fit2))
pp2 = fitted(fit2)
print(p1)
print(p2)
print(fit1)
print(fit2)
print(rr1)
print(rr2)
```

```
> print(p1)
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3 mprage_QC
deter    0.3828907             0.2333736   0.2604044        0.4003971          0.3339381             0.3739630             0.5228205 0.6061446
marked   0.1018484             0.0548929   0.1381211        0.1536405          0.3888484             0.7775766             0.6703084 0.6208475
> print(p2)
       (Intercept) scale(inatt_vol_lh.x)
deter   0.04058375            0.17085696
marked  0.38821970            0.07587068
> print(fit1)
Call:
multinom(formula = OLS_inatt_categ ~ scale(inatt_vol_lh.x) + 
    age_at_scan + I(age_at_scan^2) + Sex...Subjects + ext_avg_freesurfer5.3 + 
    int_avg_freesurfer5.3 + mprage_QC, data = df, na.action = na.omit)

Coefficients:
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3  mprage_QC
deter     3.866763            -0.3087101  -0.9746556       0.03787411          0.5063154             0.6252584             0.4661391 -0.2744336
marked    7.295847             0.4806992  -1.2774435       0.06239427         -0.4437234            -0.2108567            -0.3215006  0.2800994

Residual Deviance: 279.0317 
AIC: 311.0317 
> print(fit2)
Call:
multinom(formula = OLS_inatt_categ ~ scale(inatt_vol_lh.x), data = df, 
    na.action = na.omit)

Coefficients:
       (Intercept) scale(inatt_vol_lh.x)
deter    0.4423328            -0.3222614
marked   0.1966899             0.4105336

Residual Deviance: 294.0507 
AIC: 302.0507 
> print(rr1)
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3 mprage_QC
deter     47.78743             0.7343936   0.3773223         1.038600          1.6591666             1.8687287             1.5938286 0.7600025
marked  1474.16557             1.6172047   0.2787490         1.064382          0.6416429             0.8098901             0.7250602 1.3232614
> print(rr2)
       (Intercept) scale(inatt_vol_lh.x)
deter     1.556334             0.7245088
marked    1.217367             1.5076220

```

I'd have expected the brain variable to be more significant. It's trending, but
not there yet. Likely because the significance we got in the regression didn't
transfer to the categorical model. We should also put in the results of a
similar model, but using the regressions:

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
library(nlme)
df$wInatt = winsorize(df$OLS_inatt_slope.x)
fm = as.formula("wInatt ~ inatt_vol_lh.x + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC + age_at_scan + I(age_at_scan^2)")
fit3 = lme(fm, random=~1|nuclearFamID.x, data=df)
p3 = predict(fit3)
m = mean(df$wInatt)
fit4 = lm(wInatt ~ inatt_vol_lh.x, data=df)
p4 = predict(fit4)
```

```
> mean(abs(df$wInatt - p3))
[1] 0.4527394
> mean(abs(df$wInatt - p4))
[1] 0.4709515
> mean(abs(df$wInatt - m))
[1] 0.4983164
```

So, what does it all mean?

* The models that predict inatt winsorized OLS performs slightly better than
  chance (i.e. using the mean for prediction). Still not great though, with an
  MAE of .45 for the full model, and .47 for the model that uses only brain
  data.
* When using the categories, the model with just brain data performs better. It
  was not significant (alpha = .05) in either one, though.
* Interpreting the coefficients in the best model, a one-unit increase in the
  volume of the left hemisphere cluster variable (i.e. increase by 1 SD) is
  associated with an increase in the log odds of deteriorating vs mild
  improvement in the amount of .17. That one-unit increase is also associated
  with a .08 increase in the log odds of a marked improvement vs mild
  improvement.
* In terms of relative risk ratio, a one-unit increase in the the volume of that
  brain cluster yields a relative risk ratio of .72 for deteriorating vs. mild
  improvement. The relative risk ratio for a one-unit increase is 1.51 for
  marked improvement vs. mild improvement.
* In other words, the odds of marked improvement, compared to just mild
  improvement, is 1.5 times higher for every 1 SD we increase in that brain
  region. However, it's not really significant.
* We can also look at a plot of the predicted probabilities of being in each
  category, given the brain volume:

```r
a = cbind(pp2, df$inatt_vol_lh.x)
colnames(a)[4] = 'brain'
lpp = melt(as.data.frame(a), value.name='probability', id.vars=c('brain'))
ggplot(lpp, aes(x=brain, y=probability, color=variable)) + geom_line()
```

![](2019-01-07-15-21-30.png)

So, it doesn't look great overall. I know these are probabilities, but I'd
expect at some point for it to be more likely to be in the mild group. At least
it increases in the middle and not in the extremes, which I guess it;s somewhat
expected given how we found the brain clusters.

Let's take a look at the HI cluster now.

## struct HI

```r
imaging[imaging$OLS_HI_slope <= -.5, 'OLS_HI_categ'] = 'marked'
imaging[imaging$OLS_HI_slope > -.5 & imaging$OLS_HI_slope <= 0, 'OLS_HI_categ'] = 'mild'
imaging[imaging$OLS_HI_slope > 0, 'OLS_HI_categ'] = 'deter'
imaging$OLS_HI_categ = as.factor(imaging$OLS_HI_categ)
imaging$OLS_HI_categ = relevel(imaging$OLS_HI_categ, ref='mild')

df = merge(struct, imaging, by='MRN')
fit1 <- multinom(OLS_HI_categ ~ scale(HI_vol_rh.x) + age_at_scan + I(age_at_scan^2) + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC, data = df, na.action=na.omit)
z1 <- summary(fit1)$coefficients/summary(fit1)$standard.errors
p1 <- (1 - pnorm(abs(z1), 0, 1)) * 2
rr1 = exp(coef(fit1))
pp1 = fitted(fit1)
fit2 <- multinom(OLS_HI_categ ~ scale(HI_vol_rh.x), data = df, na.action=na.omit)
z2 <- summary(fit2)$coefficients/summary(fit2)$standard.errors
p2 <- (1 - pnorm(abs(z2), 0, 1)) * 2
rr2 = exp(coef(fit2))
pp2 = fitted(fit2)
print(p1)
print(p2)
print(fit1)
print(fit2)
print(rr1)
print(rr2)
```

```
> print(p1)
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3 mprage_QC
deter    0.3828907             0.2333736   0.2604044        0.4003971          0.3339381             0.3739630             0.5228205 0.6061446
marked   0.1018484             0.0548929   0.1381211        0.1536405          0.3888484             0.7775766             0.6703084 0.6208475
> print(p2)
       (Intercept) scale(inatt_vol_lh.x)
deter   0.04058375            0.17085696
marked  0.38821970            0.07587068
> print(fit1)
Call:
multinom(formula = OLS_inatt_categ ~ scale(inatt_vol_lh.x) + 
    age_at_scan + I(age_at_scan^2) + Sex...Subjects + ext_avg_freesurfer5.3 + 
    int_avg_freesurfer5.3 + mprage_QC, data = df, na.action = na.omit)

Coefficients:
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3  mprage_QC
deter     3.866763            -0.3087101  -0.9746556       0.03787411          0.5063154             0.6252584             0.4661391 -0.2744336
marked    7.295847             0.4806992  -1.2774435       0.06239427         -0.4437234            -0.2108567            -0.3215006  0.2800994

Residual Deviance: 279.0317 
AIC: 311.0317 
> print(fit2)
Call:
multinom(formula = OLS_inatt_categ ~ scale(inatt_vol_lh.x), data = df, 
    na.action = na.omit)

Coefficients:
       (Intercept) scale(inatt_vol_lh.x)
deter    0.4423328            -0.3222614
marked   0.1966899             0.4105336

Residual Deviance: 294.0507 
AIC: 302.0507 
> print(rr1)
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3 mprage_QC
deter     47.78743             0.7343936   0.3773223         1.038600          1.6591666             1.8687287             1.5938286 0.7600025
marked  1474.16557             1.6172047   0.2787490         1.064382          0.6416429             0.8098901             0.7250602 1.3232614
> print(rr2)
       (Intercept) scale(inatt_vol_lh.x)
deter     1.556334             0.7245088
marked    1.217367             1.5076220

```

I'd have expected the brain variable to be more significant. It's trending, but
not there yet. Likely because the significance we got in the regression didn't
transfer to the categorical model. We should also put in the results of a
similar model, but using the regressions:

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
library(nlme)
df$wInatt = winsorize(df$OLS_inatt_slope.x)
fm = as.formula("wInatt ~ inatt_vol_lh.x + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC + age_at_scan + I(age_at_scan^2)")
fit3 = lme(fm, random=~1|nuclearFamID.x, data=df)
p3 = predict(fit3)
m = mean(df$wInatt)
fit4 = lm(wInatt ~ inatt_vol_lh.x, data=df)
p4 = predict(fit4)
```

```
> mean(abs(df$wInatt - p3))
[1] 0.4527394
> mean(abs(df$wInatt - p4))
[1] 0.4709515
> mean(abs(df$wInatt - m))
[1] 0.4983164
```

So, what does it all mean?

* The models that predict inatt winsorized OLS performs slightly better than
  chance (i.e. using the mean for prediction). Still not great though, with an
  MAE of .45 for the full model, and .47 for the model that uses only brain
  data.
* When using the categories, the model with just brain data performs better. It
  was not significant (alpha = .05) in either one, though.
* Interpreting the coefficients in the best model, a one-unit increase in the
  volume of the left hemisphere cluster variable (i.e. increase by 1 SD) is
  associated with an increase in the log odds of deteriorating vs mild
  improvement in the amount of .17. That one-unit increase is also associated
  with a .08 increase in the log odds of a marked improvement vs mild
  improvement.
* In terms of relative risk ratio, a one-unit increase in the the volume of that
  brain cluster yields a relative risk ratio of .72 for deteriorating vs. mild
  improvement. The relative risk ratio for a one-unit increase is 1.51 for
  marked improvement vs. mild improvement.
* In other words, the odds of marked improvement, compared to just mild
  improvement, is 1.5 times higher for every 1 SD we increase in that brain
  region. However, it's not really significant.
* We can also look at a plot of the predicted probabilities of being in each
  category, given the brain volume:

```r
a = cbind(pp2, df$inatt_vol_lh.x)
colnames(a)[4] = 'brain'
lpp = melt(as.data.frame(a), value.name='probability', id.vars=c('brain'))
ggplot(lpp, aes(x=brain, y=probability, color=variable)) + geom_line()
```

![](2019-01-07-15-21-30.png)

So, it doesn't look great overall. I know these are probabilities, but I'd
expect at some point for it to be more likely to be in the mild group. At least
it increases in the middle and not in the extremes, which I guess it;s somewhat
expected given how we found the brain clusters.

# 2019-01-09 14:36:29

After chatting with Philip, it makes sense to use the NVs here for contrast. So,
let's redo the whole group contrast analysis and see what we get (we can keep
the regression the same way, with only SX >= 3).

```r
library(gdata)
library(nnet)

load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz')
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
idx = df$diag_group != 'new_onset'
idx2 = !is.na(df$inatt_vol_lh) | !is.na(df$inatt_AD_clu1) | !is.na(df$inatt_melodic_DMN)
imaging = df[idx & idx2,]
imaging[imaging$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
imaging[imaging$OLS_inatt_slope > -.33 & imaging$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
imaging[imaging$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
imaging[imaging$DX == 'NV', 'OLS_inatt_categ'] = 'NV'
imaging$OLS_inatt_categ = as.factor(imaging$OLS_inatt_categ)
imaging$OLS_inatt_categ = relevel(imaging$OLS_inatt_categ, ref='NV')

load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz')
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
idx = df$diag_group != 'new_onset'
struct = df[!is.na(df$HI_vol_rh) & idx,]
load('/data/NCR_SBRB/baseline_prediction/struct_volume_11142018_260timeDiff12mo.RData.gz')
struct = merge(struct, data, by='MRN') # put mask ids in combined dataset
mprage = read.xls('/data/NCR_SBRB/baseline_prediction/long_scans_08072018.xlsx',
                  sheet='mprage')
struct = merge(struct, mprage, by.x='mask.id', by.y='Mask.ID...Scan') # get demographics
qc = read.csv('/data/NCR_SBRB/baseline_prediction/master_qc.csv')
struct = merge(struct, qc, by.x='mask.id', by.y='Mask.ID') # get QC scores
df = merge(struct, imaging, by='MRN')
fit1 <- multinom(OLS_inatt_categ ~ scale(inatt_vol_lh.x) + age_at_scan + I(age_at_scan^2) + Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC, data = df, na.action=na.omit)
z1 <- summary(fit1)$coefficients/summary(fit1)$standard.errors
p1 <- (1 - pnorm(abs(z1), 0, 1)) * 2
rr1 = exp(coef(fit1))
pp1 = fitted(fit1)
fit2 <- multinom(OLS_inatt_categ ~ scale(inatt_vol_lh.x), data = df, na.action=na.omit)
z2 <- summary(fit2)$coefficients/summary(fit2)$standard.errors
p2 <- (1 - pnorm(abs(z2), 0, 1)) * 2
rr2 = exp(coef(fit2))
pp2 = fitted(fit2)
print(p1)
print(p2)
print(fit1)
print(fit2)
print(rr1)
print(rr2)
```

```
> print(p1)
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3 mprage_QC
deter  0.022533183            0.04738121  0.41395270       0.24891252        0.002623146             0.1023966            0.01002034 0.7713876
marked 0.337080723            0.02818612  0.93171828       0.89479522        0.494204426             0.8607885            0.27814406 0.4239574
mild   0.006545162            0.60732068  0.05893088       0.06125219        0.101269117             0.6581011            0.12028105 0.8247601
> print(p2)
        (Intercept) scale(inatt_vol_lh.x)
deter  8.712702e-04            0.09716676
marked 9.768048e-06            0.01820195
mild   2.037821e-07            0.94229659
> print(fit1)
Call:
multinom(formula = OLS_inatt_categ ~ scale(inatt_vol_lh.x) + 
    age_at_scan + I(age_at_scan^2) + Sex...Subjects + ext_avg_freesurfer5.3 + 
    int_avg_freesurfer5.3 + mprage_QC, data = df, na.action = na.omit)

Coefficients:
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3  mprage_QC
deter    -6.899292            -0.4151225  0.46099773      -0.03412653          1.2462371            0.87482698             1.4453837 -0.1236049
marked   -2.875338             0.4059156  0.04827447      -0.00376367          0.2702481            0.09481253             0.6029639  0.3579345
mild    -10.399575            -0.1157071  1.36308372      -0.06779776          0.7354048            0.26550164             0.9610721  0.1064069

Residual Deviance: 583.6821 
AIC: 631.6821 
> print(fit2)
Call:
multinom(formula = OLS_inatt_categ ~ scale(inatt_vol_lh.x), data = df, 
    na.action = na.omit)

Coefficients:
       (Intercept) scale(inatt_vol_lh.x)
deter   -0.5649501           -0.30315866
marked  -0.8210424            0.41183896
mild    -1.0113611            0.01456543

Residual Deviance: 620.0299 
AIC: 632.0299 
> print(rr1)
        (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale ext_avg_freesurfer5.3 int_avg_freesurfer5.3 mprage_QC
deter  1.008499e-03             0.6602594    1.585655        0.9664492           3.477234              2.398460              4.243480  0.883729
marked 5.639710e-02             1.5006759    1.049459        0.9962434           1.310289              1.099453              1.827527  1.430372
mild   3.044543e-05             0.8907361    3.908227        0.9344494           2.086326              1.304085              2.614498  1.112274
> print(rr2)
       (Intercept) scale(inatt_vol_lh.x)
deter    0.5683885             0.7384819
marked   0.4399728             1.5095913
mild     0.3637236             1.0146720
```

With more subjects we start seeing more interesting things. Some of the
covariates start approaching significance, so we can play a bit more with them:

```
> fitx <- multinom(OLS_inatt_categ ~ scale(inatt_vol_lh.x) + age_at_scan + I(age_at_scan^2) + Sex...Subjects + int_avg_freesurfer5.3, data = df, na.action=na.omit)
# weights:  28 (18 variable)
initial  value 334.096941 
iter  10 value 301.149045
iter  20 value 294.358135
iter  30 value 293.727260
final  value 293.727240 
converged
> zx <- summary(fitx)$coefficients/summary(fitx)$standard.errors
> (1 - pnorm(abs(zx), 0, 1)) * 2
       (Intercept) scale(inatt_vol_lh.x) age_at_scan I(age_at_scan^2) Sex...SubjectsMale int_avg_freesurfer5.3
deter  0.051819699            0.03447798  0.44552792       0.29232510        0.001060536           0.004887488
marked 0.404715564            0.03039763  0.96011932       0.91086332        0.484470738           0.136182084
mild   0.007982884            0.57712954  0.06031023       0.06380006        0.094448976           0.076196413
```

That puts the difference between deterioration and normals, as well as marked
improvement and nvs, significant. What are the relative risks, and probability
plots?

```r
library(reshape2)
library(ggplot2)
myfit <- multinom(OLS_inatt_categ ~ scale(inatt_vol_lh.x) + age_at_scan + I(age_at_scan^2) + Sex...Subjects + int_avg_freesurfer5.3, data = df, na.action=na.omit)
myz <- summary(myfit)$coefficients/summary(myfit)$standard.errors
myp = (1 - pnorm(abs(zx), 0, 1)) * 2
myrr = exp(coef(myfit))
mypp = fitted(myfit)
a = cbind(mypp, df$inatt_vol_lh.x)
colnames(a)[4] = 'brain'
mylpp = melt(as.data.frame(a), value.name='probability', id.vars=c('brain'))
ggplot(mylpp, aes(x=brain, y=probability, color=variable)) + geom_line()
print(myfit)
print(myp)
print(myrr)
```

Moving it to RStudio...