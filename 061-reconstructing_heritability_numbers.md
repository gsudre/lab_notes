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



# TODO
* redo figure 1
* redo figure 3
* supplemental table 1
* add comorbidities and med info to both file
