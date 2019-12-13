# 2019-12-11 10:29:04

Now we redo all the numbers in the heritability paper.

```r
dti = read.csv('~/data/heritability_change/dti_JHUtracts_ADRDonly_OD0.95_withBaseAge_clean.csv')
fmri = read.csv('~/data/heritability_change/rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean.csv')
both = merge(dti, fmri, by='ID', all.x=T, all.y=T)
write.csv(both, file='both_modalities.csv', row.names=F)
```

I then used the maximum absolute slope in SX to define the overall SX slope for
each subject (max between the two modalities, or just the one measure if only
one present). 

```r
> sx = read.csv('~/data/heritability_change/sx_both.csv')
> dim(sx)
[1] 288   3
> t.test(sx$SX_HI)

        One Sample t-test

data:  sx$SX_HI
t = -4.0865, df = 287, p-value = 5.689e-05
alternative hypothesis: true mean is not equal to 0
95 percent confidence interval:
 -0.4601992 -0.1609980
sample estimates:
 mean of x 
-0.3105986 

> sd(sx$SX_HI)/sqrt(nrow(sx))
[1] 0.07600637
> idx = sx$DX=='ADHD'
> t.test(sx[idx,]$SX_HI)

        One Sample t-test

data:  sx[idx, ]$SX_HI
t = -3.6086, df = 140, p-value = 0.0004276
alternative hypothesis: true mean is not equal to 0
95 percent confidence interval:
 -0.7561073 -0.2208572
sample estimates:
 mean of x 
-0.4884822 

> sd(sx[idx,]$SX_HI)/sqrt(sum(idx))
[1] 0.1353656
> t.test(sx$SX_inatt)

        One Sample t-test

data:  sx$SX_inatt
t = -2.1358, df = 287, p-value = 0.03355
alternative hypothesis: true mean is not equal to 0
95 percent confidence interval:
 -0.36003811 -0.01469342
sample estimates:
 mean of x 
-0.1873658 

> sd(sx$SX_inatt)/sqrt(nrow(sx))
[1] 0.08772824
> t.test(sx[idx,]$SX_inatt)

        One Sample t-test

data:  sx[idx, ]$SX_inatt
t = -3.7976, df = 140, p-value = 0.000217
alternative hypothesis: true mean is not equal to 0
95 percent confidence interval:
 -0.8518049 -0.2685432
sample estimates:
 mean of x 
-0.5601741 

> sd(sx[idx,]$SX_inatt)/sqrt(sum(idx))
[1] 0.1475078
```

Then, grabbed Prop_on_psychostim    Med_binary  comorbid_binary from
subj_meds.csv and filled in the remaining ones that weren't there in the
previous cohort. There were 3 ADHD kids I couldn't tell, so I'm leaving them out
from this analysis. If the results flip because of them, they weren't too strong
to begin with... just the prop variable, though.

Just need to re-run the fMRI results:

```r
both = read.csv('~/data/heritability_change/subj_meds_new.csv')
fmri = read.csv('~/data/heritability_change/rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean.csv')
fmri2 = merge(fmri, both, by='ID', all.x=F, all.y=F)
write.csv(fmri2, file='~/data/heritability_change/rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds.csv', row.names=F, quote=F, na='')
fmri2 = fmri2[fmri2$comorbid_binary==0,]
write.csv(fmri2, file='~/data/heritability_change/rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds_nocomorbid.csv', row.names=F, quote=F, na='')
```

```bash
cd ~/data/heritability_change/;
for p in conn_SalVentAttnTOCont conn_DorsAttnTOSalVentAttn; do
    for c in Prop_on_psychostim Med_binary comorbid_binary; do
        solar run_phen_var_OD_xcp_meds rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds $p $c;
     done;
     solar run_phen_var_OD_xcp_meds_comb3 rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds $p;
done
for p in conn_SalVentAttnTOCont conn_DorsAttnTOSalVentAttn; do
    for c in Prop_on_psychostim Med_binary; do
        solar run_phen_var_OD_xcp_meds rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds_nocomorbid $p $c;
     done;
     solar run_phen_var_OD_xcp rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds_nocomorbid ${p}
     solar run_phen_var_OD_xcp_meds_comb2 rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean_meds_nocomorbid $p;
done
```

We should now check if our heritability values from before are still there.
Hand-copied the results to medication_results.xlsx.

# 2019-12-12 15:56:54

Let's redo the SX analysis using all (possibly 4) time points for each subject:

```r
dti = read.csv('~/philip/MS_h2_developing_connectivity_11242019/new_fmri_pheno/dti_JHUtracts_ADRDonly_OD0.95_twoTimePoints_noOtherDX.csv')
fmri = read.csv('~/philip/MS_h2_developing_connectivity_11242019/new_fmri_pheno/rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_twoTimePoints.csv')
bothd = dti[, c(1, 2, 3, 4, 9, 11, 12, 95, 96, 97, 98)]
bothf = fmri[, c(1, 2, 40, 41, 45, 46, 47, 51, 52, 53, 54)]
colnames(bothd) = c('MRN', 'maskid', 'DOS', 'sex', 'age_at_scan', 'extid', 
'nucid', 'DOA', 'SX_inatt', 'SX_HI', 'source')
colnames(bothf) = c('MRN', 'maskid', 'DOS', 'sex', 'age_at_scan', 'extid', 
'nucid', 'DOA', 'SX_inatt', 'SX_HI', 'source')
both = rbind(bothd, bothf)
both = both[!duplicated(both$maskid), ]
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(both, tmp[, c('ID', 'FAMID')], by.x='MRN', by.y='ID', all.x=T, all.y=F)
write.csv(data, file='~/philip/MS_h2_developing_connectivity_11242019/new_fmri_pheno/both_modalities_twoTimePoints.csv', row.names=F)
```

```
> fit = lme(as.formula('SX_HI ~ age_at_scan + sex'), data, ~1|FAMID/MRN)
> summary(fit)
Linear mixed-effects model fit by REML
 Data: data
       AIC      BIC    logLik
  3242.868 3270.573 -1615.434
Random effects:
 Formula: ~1 | FAMID
        (Intercept)
StdDev:    1.679479
 Formula: ~1 | MRN %in% FAMID
        (Intercept) Residual
StdDev:     1.66506 1.404346
Fixed effects: as.formula("SX_HI ~ age_at_scan + sex")
                Value Std.Error  DF    t-value p-value
(Intercept)  5.503084 0.4060251 462  13.553555  0.0000
age_at_scan -0.304215 0.0286709 462 -10.610583  0.0000
sexMale      0.576026 0.2933010  82   1.963941  0.0529
 Correlation:
            (Intr) ag_t_s
age_at_scan -0.778
sexMale     -0.486  0.010
Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max
-2.47059187 -0.46180421 -0.07843612  0.41657929  2.95927525
Number of Observations: 751
Number of Groups:
         FAMID MRN %in% FAMID
           205            288
> fit = lme(as.formula('SX_inatt ~ age_at_scan + sex'), data, ~1|FAMID/MRN)
> summary(fit)
Linear mixed-effects model fit by REML
 Data: data
       AIC      BIC   logLik
  3367.319 3395.024 -1677.66
Random effects:
 Formula: ~1 | FAMID
        (Intercept)
StdDev:    2.546319
 Formula: ~1 | MRN %in% FAMID
        (Intercept) Residual
StdDev:    1.533805 1.463071
Fixed effects: as.formula("SX_inatt ~ age_at_scan + sex")
                Value Std.Error  DF   t-value p-value
(Intercept)  4.408017 0.4490113 462  9.817162  0.0000
age_at_scan -0.060154 0.0303247 462 -1.983654  0.0479
sexMale      0.425959 0.3187100  82  1.336510  0.1851
 Correlation:
            (Intr) ag_t_s
age_at_scan -0.747
sexMale     -0.483  0.015
Standardized Within-Group Residuals:
       Min         Q1        Med         Q3        Max
-3.1502905 -0.3586211 -0.1147428  0.4465651  3.0862429
Number of Observations: 751
Number of Groups:
         FAMID MRN %in% FAMID
           205            288
```

Now we create the DX2 subgroup:

```r
dx2 = c()
for (s in unique(data$MRN)) {
    idx = which(data$MRN == s)
    # grabbing inatt and HI at baseline
    base_DOA = which.min(data[idx, 'age_at_scan'])
    base_inatt = data[idx[base_DOA], 'SX_inatt']
    base_HI = data[idx[base_DOA], 'SX_HI']
    if ((base_inatt >= 4) || (base_HI >= 4)) {
        dx2 = c(dx2, rep("ADHD", length(idx)))
    } else {
        dx2 = c(dx2, rep("NV", length(idx)))
    }
}
data$DX2 = dx2
```

```
> fit = lme(as.formula('SX_HI ~ age_at_scan + sex'), data[data$DX2=='ADHD',], ~1|FAMID/MRN)
> summary(fit)
Linear mixed-effects model fit by REML
 Data: data[data$DX2 == "ADHD", ]
       AIC      BIC    logLik
  1913.476 1937.887 -950.7382
Random effects:
 Formula: ~1 | FAMID
         (Intercept)
StdDev: 0.0004258675
 Formula: ~1 | MRN %in% FAMID
        (Intercept) Residual
StdDev:    1.884903 1.608414
Fixed effects: as.formula("SX_HI ~ age_at_scan + sex")
                Value Std.Error  DF    t-value p-value
(Intercept)  9.089354 0.5143013 268  17.673209  0.0000
age_at_scan -0.494384 0.0391790 268 -12.618602  0.0000
sexMale      0.378327 0.3603042  24   1.050022  0.3042
 Correlation:
            (Intr) ag_t_s
age_at_scan -0.813       
sexMale     -0.467 -0.022
Standardized Within-Group Residuals:
         Min           Q1          Med           Q3          Max 
-2.177982144 -0.587895031 -0.004107512  0.567162521  2.663937249 
Number of Observations: 435
Number of Groups: 
         FAMID MRN %in% FAMID 
           141            166 
> fit = lme(as.formula('SX_inatt ~ age_at_scan + sex'), data[data$DX2=='ADHD',], ~1|FAMID/MRN)
> summary(fit)
Linear mixed-effects model fit by REML
 Data: data[data$DX2 == "ADHD", ] 
       AIC      BIC    logLik
  1900.221 1924.632 -944.1106
Random effects:
 Formula: ~1 | FAMID
        (Intercept)
StdDev:   0.5901345
 Formula: ~1 | MRN %in% FAMID
        (Intercept) Residual
StdDev:    1.505758   1.6719
Fixed effects: as.formula("SX_inatt ~ age_at_scan + sex")
                Value Std.Error  DF   t-value p-value
(Intercept)  6.989247 0.4987131 268 14.014563  0.0000
age_at_scan -0.089064 0.0390930 268 -2.278254  0.0235
sexMale     -0.097122 0.3243734  24 -0.299414  0.7672
 Correlation:
            (Intr) ag_t_s
age_at_scan -0.838
sexMale     -0.431 -0.024
Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max
-2.40664796 -0.51373159  0.07207311  0.62022010  2.44177407
Number of Observations: 435
Number of Groups:
         FAMID MRN %in% FAMID
           141            166
```

Now, let's see if the heritability p-values survive FDR across modalities:

```r
dtips = read.csv('~/philip/MS_h2_developing_connectivity_11242019/new_fmri_pheno/polygen_results_dti_JHUtracts_ADRDonly_OD0.95.csv')
fmrips = read.csv('~/philip/MS_h2_developing_connectivity_11242019/new_fmri_pheno/fmri_pheno_h2r.csv')
ps = c(dtips$h_pval, fmrips$h_pval)
p2 = p.adjust(ps, method='fdr')
res = cbind(c(as.character(dtips$phen), as.character(fmrips$phen)), ps, p2)
colnames(res) = c('phen', 'h_pval', 'fdr')
write.csv(res, file='~/tmp/both_ps.csv', row.names=F, quote=F)
```

## Meff

```r
phen = 'rsfmri_7by7from100_4nets_p05SigSum_OD0.95_12052019_clean'
var_names = c("conn_DorsAttnTODorsAttn", "conn_DorsAttnTOSalVentAttn",
              "conn_DorsAttnTOCont", "conn_DorsAttnTODefault", "conn_SalVentAttnTOSalVentAttn", "conn_SalVentAttnTOCont",
              "conn_SalVentAttnTODefault", "conn_ContTOCont",
              "conn_ContTODefault", "conn_DefaultTODefault")
fname = sprintf('~/data/heritability_change/%s.csv', phen)
data = read.csv(fname)
cc = cor(data[, var_names], use='na.or.complete')
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

## Plotting h2 matrix

```r
cc = read.csv('~/tmp/fmri_h2r.csv')
rownames(cc) = cc$net
cc$net=NULL
corrplot(as.matrix(cc), method='color', type='upper', is.corr=F)
```

# TODO
* redo figure 1
* redo figure 3
* supplemental table 1
* add comorbidities and med info to both file
