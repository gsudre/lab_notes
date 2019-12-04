# 2019-12-02 11:05:01

I just noticed a bug in how the fmri connectivity was being calculated. So,
let's redo all the heritability estimates:

```bash
Rscript ~/research_code/fmri/make_outlier_detection_slopes.R .9
```

I ran that for the different thresholds, but only median and positive for now.

```bash
# locally
for OD in 80 85 90 95; do
	suf='median';
	p='';
	cd ~/data/heritability_change
	phen=rsfmri_7by7from100_5nets_OD0.${OD}${p}_${suf}_12022019;
	for t in "conn_DorsAttnTODorsAttn" \
		"conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
		"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
		"conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
		"conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
		"conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
		"conn_DefaultTODefault"; do
			solar run_phen_var_OD_xcp ${phen} ${t};
	done;
	mv ${phen} ~/data/tmp/
	cd ~/data/tmp/${phen}
	for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
	python ~/research_code/compile_solar_multivar_results.py ${phen}
done
```

Let's glue everything together to see if there is a pattern

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_5nets.csv;
for f in `ls polygen_results_*5nets*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_5nets.csv;
	done
done
```

Not much there. Let me add mean and all connections variants to see if anything
comes up.

```bash
# locally
for OD in 80 85 90 95; do
	for suf in 'median' 'mean'; do
		for p in '' 'All'; do
			cd ~/data/heritability_change
			phen=rsfmri_7by7from100_5nets_OD0.${OD}_${suf}${p}_12022019;
			for t in "conn_DorsAttnTODorsAttn" \
				"conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
				"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
				"conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
				"conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
				"conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
				"conn_DefaultTODefault"; do
					solar run_phen_var_OD_xcp ${phen} ${t};
			done;
			mv ${phen} ~/data/tmp/
			cd ~/data/tmp/${phen}
			for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
			python ~/research_code/compile_solar_multivar_results.py ${phen}
		done;
	done;
done
```

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_5nets.csv;
for f in `ls polygen_results_*5nets*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_5nets.csv;
	done
done
```

![](images/2019-12-02-11-46-32.png)

It looks like if we go with
polygen_results_rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019.csv we have
some interesting results. The picture shows everything < .05, then I just sorted
it on phen name. I then copied those polygen files to the main results folder
just in case. I also removed any phenotypes that had h2r=1.

I'm also going to try a version with only 4 networks too, just in case 5nets doesn't
work and I'm still waiting for SOLAR to run above as well. 

```bash
# locally
for OD in 80 85 90 95; do
	for suf in 'median' 'mean'; do
		for p in '' 'All'; do
			cd ~/data/heritability_change
			phen=rsfmri_7by7from100_4nets_OD0.${OD}_${suf}${p}_12022019;
			for t in "conn_DorsAttnTODorsAttn" \
				"conn_DorsAttnTOSalVentAttn"  "conn_DorsAttnTOCont" \
				"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
				 "conn_SalVentAttnTOCont" "conn_SalVentAttnTODefault"  \
				"conn_ContTOCont" "conn_ContTODefault" "conn_DefaultTODefault"; do
					solar run_phen_var_OD_xcp ${phen} ${t};
			done;
			mv ${phen} ~/data/tmp/
			cd ~/data/tmp/${phen}
			for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
			python ~/research_code/compile_solar_multivar_results.py ${phen}
		done;
	done;
done
```

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_4nets.csv;
for f in `ls polygen_results_*4nets*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_4nets.csv;
	done
done
```

![](images/2019-12-02-11-52-10.png)

Here we do slightly better using
polygen_results_rsfmri_7by7from100_4nets_OD0.85_mean_12022019.csv. I'll stick
with 5nets for now, as it gives us a bigger number of subjects. Let's just make
sure that it survives FDR, and then that there is still association.

![](images/2019-12-02-11-55-38.png)

The only one that survives FDR q < .05 (or .1) is conn_DorsAttnTOSalVentAttn.
But nothing in the 4nets situation. Let's then check whether that connection is
also associated with SX.

```
HG-02113362-DM4:rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep "(Significant)" conn_DorsAttnTOSalVentAttn_polygenic.out
						 H2r is 0.8070358  p = 0.0007376  (Significant)
							  pctSpikesDV  p = 0.0070319  (Significant)
						 motionDVCorrInit  p = 0.0626968  (Significant)
```

```r
library(nlme)
data = read.csv('~/data/heritability_change/rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(data, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)

i = 'conn_DorsAttnTOSalVentAttn'
fm_root = '%s ~ %s + pctSpikesDV + motionDVCorrInit'
i = 'conn_DorsAttnTODefault'
fm_root = '%s ~ %s + pctSpikesDV'
# i = 'conn_SalVentAttnTOLimbic'
# fm_root = '%s ~ %s + pctSpikesDV'
# i = 'conn_DorsAttnTODorsAttn'
# fm_root = '%s ~ %s + pctSpikesDV + normCoverage + relMeanRMSMotion'
# i = 'conn_SalVentAttnTOSalVentAttn'
# fm_root = '%s ~ %s + pctSpikesDV + motionDVCorrInit + normCoverage + motionDVCorrFinal + pctSpikesRMS'

out_fname = sprintf('~/data/heritability_change/assoc_%s.csv', i)
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline', 'DX', 'DX2')
hold=NULL
for (j in predictors) {
	fm_str = sprintf(fm_root, i, j)
	model1<-try(lme(as.formula(fm_str), data, ~1|FAMID, na.action=na.omit))
	if (length(model1) > 1) {
		temp<-summary(model1)$tTable
		a<-as.data.frame(temp)
		a$formula<-fm_str
		a$target = i
		a$predictor = j
		a$term = rownames(temp)
		hold=rbind(hold,a)
	} else {
		hold=rbind(hold, NA)
	}
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX=='ADHD', ]
out_fname = gsub(x=out_fname, pattern='.csv', '_dx1.csv')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
hold=NULL
for (j in predictors) {
	fm_str = sprintf(fm_root, i, j)
	model1<-try(lme(as.formula(fm_str), data2, ~1|FAMID, na.action=na.omit))
	if (length(model1) > 1) {
		temp<-summary(model1)$tTable
		a<-as.data.frame(temp)
		a$formula<-fm_str
		a$target = i
		a$predictor = j
		a$term = rownames(temp)
		hold=rbind(hold,a)
	} else {
		hold=rbind(hold, NA)
	}
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX2=='ADHD', ]
out_fname = gsub(x=out_fname, pattern='dx1', 'dx2')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
hold=NULL
for (j in predictors) {
	fm_str = sprintf(fm_root, i, j)
	model1<-try(lme(as.formula(fm_str), data2, ~1|FAMID, na.action=na.omit))
	if (length(model1) > 1) {
		temp<-summary(model1)$tTable
		a<-as.data.frame(temp)
		a$formula<-fm_str
		a$target = i
		a$predictor = j
		a$term = rownames(temp)
		hold=rbind(hold,a)
	} else {
		hold=rbind(hold, NA)
	}
}
write.csv(hold, out_fname, row.names=F)
```

So, conn_DorsAttnTOSalVentAttn is not even close to significance, in none of the
3 types of regression I'm running. How about the other 4 connections that are
nominally heritable?

```
HG-02113362-DM4:rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep "(Significant)" conn_DorsAttnTODefault_polygenic.out
						 H2r is 0.5520737  p = 0.0149298  (Significant)
							  pctSpikesDV  p = 0.0349292  (Significant)
HG-02113362-DM4:rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep "(Significant)" conn_SalVentAttnTOLimbic_polygenic.out 
						 H2r is 0.5319386  p = 0.0163941  (Significant)
							  pctSpikesDV  p = 0.0031629  (Significant)
HG-02113362-DM4:rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep "(Significant)" conn_DorsAttnTODorsAttn_polygenic.out 
						 H2r is 0.5761180  p = 0.0214219  (Significant)
							 normCoverage  p = 0.0006030  (Significant)
							  pctSpikesDV  p = 0.0000139  (Significant)
						 relMeanRMSMotion  p = 0.0492652  (Significant)
HG-02113362-DM4:rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep "(Significant)" conn_SalVentAttnTOSalVentAttn_polygenic.out 
						 H2r is 0.5336712  p = 0.0343229  (Significant)
							 normCoverage  p = 0.0845580  (Significant)
							  pctSpikesDV  p = 0.0951218  (Significant)
						 motionDVCorrInit  p = 0.0001882  (Significant)
						motionDVCorrFinal  p = 0.0420654  (Significant)
							 pctSpikesRMS  p = 0.0098802  (Significant)
```

The only one is conn_SalVentAttnTOSalVentAttn, which is associated with SX_HI at
.025, at DX2 type.

![](images/2019-12-02-12-21-34.png)

I didn't check the other types, though. But that's the only connection at DX2.

# Genetic correlation

Let's start within fmri:

The last step before robustness is the bivariate heritability analysis.

```bash
phen=rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019;
cd ~/data/heritability_change
for t in "conn_DorsAttnTODorsAttn" \
				"conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
				"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
				"conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
				"conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
				"conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
				"conn_DefaultTODefault"; do
		solar rsfmri_xcp_base_slope_correlation $phen ${t};
done
cd gencor_$phen
grep -r RhoG */polygenic.out > rhog.txt
grep zero rhog.txt
grep "\-1" rhog.txt
```

Not much there for results. They're either 0 or -1, or at least not
significantly different than those.

```
HG-02113362-DM4:gencor_rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep zero rhog.txt
conn_ContTOCont_AND_baseconn_ContTOCont/polygenic.out:         RhoG different from zero  p = 0.0005589
conn_ContTODefault_AND_baseconn_ContTODefault/polygenic.out:           RhoG different from zero  p = 3.325347e-06
conn_DefaultTODefault_AND_baseconn_DefaultTODefault/polygenic.out:             RhoG different from zero  p = 0.9266462
conn_DorsAttnTOCont_AND_baseconn_DorsAttnTOCont/polygenic.out:         RhoG different from zero  p = 0.0693849
conn_DorsAttnTODefault_AND_baseconn_DorsAttnTODefault/polygenic.out:           RhoG different from zero  p = 0.0005720
conn_DorsAttnTODorsAttn_AND_baseconn_DorsAttnTODorsAttn/polygenic.out:         RhoG different from zero  p = 0.0923875
conn_DorsAttnTOLimbic_AND_baseconn_DorsAttnTOLimbic/polygenic.out:             RhoG different from zero  p = 0.1924130
conn_DorsAttnTOSalVentAttn_AND_baseconn_DorsAttnTOSalVentAttn/polygenic.out:           RhoG different from zero  p = 0.0025333
conn_LimbicTOCont_AND_baseconn_LimbicTOCont/polygenic.out:             RhoG different from zero  p = 0.7883320
conn_LimbicTODefault_AND_baseconn_LimbicTODefault/polygenic.out:               RhoG different from zero  p = 0.2135416
conn_LimbicTOLimbic_AND_baseconn_LimbicTOLimbic/polygenic.out:         RhoG different from zero  p = 1.0000000
conn_SalVentAttnTOCont_AND_baseconn_SalVentAttnTOCont/polygenic.out:           RhoG different from zero  p = 0.0320647
conn_SalVentAttnTODefault_AND_baseconn_SalVentAttnTODefault/polygenic.out:             RhoG different from zero  p = 0.4222949
conn_SalVentAttnTOLimbic_AND_baseconn_SalVentAttnTOLimbic/polygenic.out:               RhoG different from zero  p = 0.0940737
conn_SalVentAttnTOSalVentAttn_AND_baseconn_SalVentAttnTOSalVentAttn/polygenic.out:             RhoG different from zero  p = 0.2564076

HG-02113362-DM4:gencor_rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019 sudregp$ grep "\-1" rhog.txt
conn_ContTOCont_AND_baseconn_ContTOCont/polygenic.out:                   RhoG is -1.0000000
conn_ContTODefault_AND_baseconn_ContTODefault/polygenic.out:                     RhoG is -1.0000000
conn_DefaultTODefault_AND_baseconn_DefaultTODefault/polygenic.out:             RhoG different from -1.0  p = 0.2370629
conn_DorsAttnTOCont_AND_baseconn_DorsAttnTOCont/polygenic.out:                   RhoG is -1.0000000
conn_DorsAttnTODefault_AND_baseconn_DorsAttnTODefault/polygenic.out:                     RhoG is -1.0000000
conn_DorsAttnTODorsAttn_AND_baseconn_DorsAttnTODorsAttn/polygenic.out:                   RhoG is -1.0000000
conn_DorsAttnTOLimbic_AND_baseconn_DorsAttnTOLimbic/polygenic.out:                       RhoG is -1.0000000
conn_DorsAttnTOSalVentAttn_AND_baseconn_DorsAttnTOSalVentAttn/polygenic.out:           RhoG different from -1.0  p = 0.1085292
conn_LimbicTOCont_AND_baseconn_LimbicTOCont/polygenic.out:                       RhoG is -1.0000000
conn_LimbicTODefault_AND_baseconn_LimbicTODefault/polygenic.out:                         RhoG is -1.0000000
conn_LimbicTOLimbic_AND_baseconn_LimbicTOLimbic/polygenic.out:         RhoG different from -1.0  p = 0.5000000
conn_SalVentAttnTOCont_AND_baseconn_SalVentAttnTOCont/polygenic.out:                     RhoG is -1.0000000
conn_SalVentAttnTODefault_AND_baseconn_SalVentAttnTODefault/polygenic.out:                       RhoG is -1.0000000
conn_SalVentAttnTOLimbic_AND_baseconn_SalVentAttnTOLimbic/polygenic.out:               RhoG different from -1.0  p = 0.0258793
conn_SalVentAttnTOSalVentAttn_AND_baseconn_SalVentAttnTOSalVentAttn/polygenic.out:             RhoG different from -1.0  p = 0.0905380
```

Now we look at inter-modality rhoG. Start by merging the data across modalities. Then, it should just be a matter of
creating the right procedure in SOLAR.

```r
fmri = read.csv('~/data/heritability_change/rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019.csv')
dti = read.csv('~/data/heritability_change/dti_JHUtracts_ADRDonly_OD0.95.csv')
both = merge(fmri, dti, by='ID', all.x=F, all.y=F)
# 155 participants with both datasets... not bad
write.csv(both, file='~/data/heritability_change/both_fmri_dti_12022019.csv', row.names=F, quote=F)
```

```bash
# locally
cd ~/data/heritability_change/
solar dti_rsfmri_slope_correlation both_fmri_dti_12022019
cd gencor_both_fmri_dti_12022019
grep -r RhoG */polygenic.out > rhog.txt
grep zero rhog.txt | cut -d "/" -f 1 > nets0.txt
grep zero rhog.txt | cut -d "=" -f 2 > ps0.txt
paste nets.txt ps.txt > cross_modal_rhog_pvals.txt;
cat cross_modal_rhog_pvals.txt | awk '{if ($2<.05) print $1 }' > rhogs_diff_0.txt;
grep -f rhogs_diff_0.txt rhog.txt | grep "different from \-1."
```

```
ad_2_AND_conn_ContTODefault/polygenic.out:             RhoG different from -1.0  p = 0.4457974
ad_2_AND_conn_SalVentAttnTOCont/polygenic.out:         RhoG different from -1.0  p = 0.3635204
ad_7_AND_conn_DorsAttnTOCont/polygenic.out:            RhoG different from -1.0  p = 0.4162783
rd_16_AND_conn_SalVentAttnTOCont/polygenic.out:        RhoG different from -1.0  p = 0.4258607
rd_2_AND_conn_ContTOCont/polygenic.out:        RhoG different from -1.0  p = 0.4814827
rd_2_AND_conn_ContTODefault/polygenic.out:             RhoG different from -1.0  p = 0.3946157
rd_2_AND_conn_SalVentAttnTOCont/polygenic.out:         RhoG different from -1.0  p = 0.2888947
rd_6_AND_conn_SalVentAttnTOCont/polygenic.out:         RhoG different from -1.0  p = 0.4879540
rd_8_AND_conn_DorsAttnTODorsAttn/polygenic.out:        RhoG different from -1.0  p = 0.4513781
```

So, those are all the combinations that have rhoG significantly differnent than
0, and which are not -1. None of them are significantly different than -1
though... I wonder if this is simply a sample issue? Just note that even within
modalities we didn't have anything significantly different than 0 and -1.

Let's look at rhoP as well:

```bash
# locally
cd ~/data/heritability_change/gencor_both_fmri_dti_12022019
grep -r RhoP */polygenic.out > rhop.txt
sed " s/ is /=/g" rhop.txt | cut -d "=" -f 2 > ests.txt
paste nets.txt ests.txt > cross_modal_rhop_estimates.txt
```

But SOLAR didn't compute the p-values, so I'll have to do that manually in R,
taking into consideration that we have 156 individuals for all pairs.

```r
n=155
df = read.table('~/data/heritability_change/gencor_both_fmri_dti_12022019/cross_modal_rhop_estimates.txt')
df$t = df$V2/sqrt((1-df$V2**2)/(n-2))
df$p = 2*pt(-abs(df$t),df=n-1)
colnames(df)[1:2] = c('traits', 'rhoP')
write.csv(df, file='~/data/heritability_change/cross_rhop_155.csv', row.names=F)
```

As I recompute everything, here are the numbers Philip asked for, which led me
to find this bug:

```
> dim(b)
[1] 69  7
> table(b[,1])

	   Cont     Default    DorsAttn      Limbic SalVentAttn      SomMot 
		 13          24          15           5          12           0 
		Vis 
		  0 
```

And I checked it now, and the code does what's supposed to do:

```
> net_names
[1] "Vis"         "SomMot"      "DorsAttn"    "SalVentAttn" "Limbic"
[6] "Cont"        "Default"
> i=6
> j=7
> idx = (conn_map[,2]==net_names[i] & conn_map[,3]==net_names[j]) |
+             (conn_map[,3]==net_names[i] & conn_map[,2]==net_names[j])
> length(var_names[idx])
[1] 312
> 13*24
[1] 312
> i=3
> j=4
> idx = (conn_map[,2]==net_names[i] & conn_map[,3]==net_names[j]) |
+             (conn_map[,3]==net_names[i] & conn_map[,2]==net_names[j])
> length(var_names[idx])
[1] 180
> 15*12
[1] 180
```

## Meff

Let's calculate Meff to see if it's a bit more forgiving then fMRI:

```r
fname = '~/data/heritability_change/rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019.csv'
data = read.csv(fname)
var_names = c("conn_DorsAttnTODorsAttn", "conn_DorsAttnTOSalVentAttn",
			  "conn_DorsAttnTOCont", "conn_DorsAttnTODefault", "conn_SalVentAttnTOSalVentAttn", "conn_SalVentAttnTOCont",
			  "conn_SalVentAttnTODefault", "conn_ContTOCont",
			  "conn_ContTODefault", "conn_DefaultTODefault",
			  "conn_DorsAttnTOLimbic", "conn_SalVentAttnTOLimbic",
			  "conn_LimbicTOLimbic", "conn_LimbicTOCont", "conn_LimbicTODefault")
cc = cor(data[, var_names])
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

So, we get 10.01, which gives us .0049 threshold. So, no difference than just
using FDR.

## 3nets?

What if I just look at the interplay between the attention networks and DMN?

```bash
# locally
for OD in 80 85 90 95; do
	for suf in 'median' 'mean'; do
		for p in '' 'All'; do
			cd ~/data/heritability_change
			phen=rsfmri_7by7from100_3nets_OD0.${OD}_${suf}${p}_12022019;
			for t in "conn_DorsAttnTODorsAttn" "conn_DorsAttnTOSalVentAttn"  \
				"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
				"conn_SalVentAttnTODefault" "conn_DefaultTODefault"; do
					solar run_phen_var_OD_xcp ${phen} ${t};
			done;
			mv ${phen} ~/data/tmp/
			cd ~/data/tmp/${phen}
			for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
			python ~/research_code/compile_solar_multivar_results.py ${phen}
		done;
	done;
done
```

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_3nets.csv;
for f in `ls polygen_results_*3nets*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_3nets.csv;
	done
done
```

![](images/2019-12-02-14-28-46.png)

Well, it actually looks quite nice if we look at
polygen_results_rsfmri_7by7from100_3nets_OD0.90_median_12022019.csv. Basically,
4 out of the 6 results are p < .05. Let's see how it goes with FDR and
association:

![](images/2019-12-02-14-31-14.png)

For starters, now none of them survive Bonferroni. But all 4 there are nominal
also survive FDR:

```
> df = read.csv('~/data/heritability_change/polygen_results_rsfmri_7by7from100_3nets_OD0.90_median_12022019.csv')
> p2 = p.adjust(df$h_pval, method='fdr')
> idx = which(p2<.05)
> df[idx,'phen']
[1] conn_DefaultTODefault         conn_DorsAttnTOSalVentAttn   
[3] conn_SalVentAttnTODefault     conn_SalVentAttnTOSalVentAttn
```

Now we can check just for association:

```
HG-02113362-DM4:rsfmri_7by7from100_3nets_OD0.90_median_12022019 sudregp$ grep "(Significant)" conn_DefaultTODefault_polygenic.out 
						 H2r is 0.6659067  p = 0.0120018  (Significant)
									  sex  p = 0.0861024  (Significant)
							  pctSpikesDV  p = 0.0024897  (Significant)
						 motionDVCorrInit  p = 0.0033996  (Significant)
HG-02113362-DM4:rsfmri_7by7from100_3nets_OD0.90_median_12022019 sudregp$ grep "(Significant)" conn_DorsAttnTOSalVentAttn_polygenic.out 
						 H2r is 0.4900904  p = 0.0330507  (Significant)
							  pctSpikesDV  p = 0.0994364  (Significant)
						motionDVCorrFinal  p = 0.0500346  (Significant)
HG-02113362-DM4:rsfmri_7by7from100_3nets_OD0.90_median_12022019 sudregp$ grep "(Significant)" conn_SalVentAttnTODefault_polygenic.out 
						 H2r is 0.7058733  p = 0.0090638  (Significant)
									  sex  p = 0.0566155  (Significant)
							  pctSpikesDV  p = 0.0104104  (Significant)
						 motionDVCorrInit  p = 0.0002635  (Significant)
HG-02113362-DM4:rsfmri_7by7from100_3nets_OD0.90_median_12022019 sudregp$ grep "(Significant)" conn_SalVentAttnTOSalVentAttn_polygenic.out 
						 H2r is 0.4515570  p = 0.0282135  (Significant)
									  sex  p = 0.0358756  (Significant)
							 normCoverage  p = 0.0681788  (Significant)
							  pctSpikesDV  p = 0.0100749  (Significant)
							 pctSpikesRMS  p = 0.0915601  (Significant)
						 relMeanRMSMotion  p = 0.0260734  (Significant)
```

```r
library(nlme)
data = read.csv('~/data/heritability_change/rsfmri_7by7from100_5nets_OD0.90_medianAll_12022019.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(data, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)

i = 'conn_DorsAttnTODefault'
fm_root = '%s ~ %s + sex + pctSpikesDV + motionDVCorrInit'
# i = 'conn_DorsAttnTOSalVentAttn'
# fm_root = '%s ~ %s + pctSpikesDV + motionDVCorrFinal'
# i = 'conn_SalVentAttnTODefault'
# fm_root = '%s ~ %s + sex + pctSpikesDV + motionDVCorrInit'
# i = 'conn_SalVentAttnTOSalVentAttn'
# fm_root = '%s ~ %s + sex + pctSpikesDV + normCoverage + pctSpikesRMS + relMeanRMSMotion'

out_fname = sprintf('~/data/heritability_change/assoc_%s.csv', i)
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline', 'DX', 'DX2')
hold=NULL
for (j in predictors) {
	fm_str = sprintf(fm_root, i, j)
	model1<-try(lme(as.formula(fm_str), data, ~1|FAMID, na.action=na.omit))
	if (length(model1) > 1) {
		temp<-summary(model1)$tTable
		a<-as.data.frame(temp)
		a$formula<-fm_str
		a$target = i
		a$predictor = j
		a$term = rownames(temp)
		hold=rbind(hold,a)
	} else {
		hold=rbind(hold, NA)
	}
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX=='ADHD', ]
out_fname = gsub(x=out_fname, pattern='.csv', '_dx1.csv')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
hold=NULL
for (j in predictors) {
	fm_str = sprintf(fm_root, i, j)
	model1<-try(lme(as.formula(fm_str), data2, ~1|FAMID, na.action=na.omit))
	if (length(model1) > 1) {
		temp<-summary(model1)$tTable
		a<-as.data.frame(temp)
		a$formula<-fm_str
		a$target = i
		a$predictor = j
		a$term = rownames(temp)
		hold=rbind(hold,a)
	} else {
		hold=rbind(hold, NA)
	}
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX2=='ADHD', ]
out_fname = gsub(x=out_fname, pattern='dx1', 'dx2')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
hold=NULL
for (j in predictors) {
	fm_str = sprintf(fm_root, i, j)
	model1<-try(lme(as.formula(fm_str), data2, ~1|FAMID, na.action=na.omit))
	if (length(model1) > 1) {
		temp<-summary(model1)$tTable
		a<-as.data.frame(temp)
		a$formula<-fm_str
		a$target = i
		a$predictor = j
		a$term = rownames(temp)
		hold=rbind(hold,a)
	} else {
		hold=rbind(hold, NA)
	}
}
write.csv(hold, out_fname, row.names=F)
```

conn_SalVentAttnTOSalVentAttn is significant for DX2 SX_HI, and so is
conn_DorsAttnTOSalVentAttn, but that's it (only looked at DX2, though).

## looking back at net4

Let's take a second look at net4. There, our best result would be with
polygen_results_rsfmri_7by7from100_4nets_OD0.85_mean_12022019.csv, which has 3
connections below .05. The universe looks like:

![](images/2019-12-02-16-49-19.png)

```
> df = read.csv('~/data/heritability_change/polygen_results_rsfmri_7by7from100_3nets_OD0.90_median_12022019.csv')
> p2 = p.adjust(df$h_pval, method='fdr')
```

Nothing... not even at q=.1. Is there at least some association? Things might
look good if we just go nominal on this one.

# 2019-12-03 21:01:38

Let's play with the idea of all networks, and then trimming the 100by100 matrix
either by positive only or by within networks. We can then try other groups of
collapsed networks.

```bash
#locally
Rscript ~/research_code/fmri/make_outlier_detection_slopes_raw.R .8
# and the other thresholds
```

The issue is that isolation forest doesn't run with NAs, so let's skip posOnly
for now. We'll just run all connections, and then we could trim it by network if
needed, instead of creating network-dependent SOLAR runs.

Then, it takes a while to run solar, but we do it like this:

```bash
# locally
for OD in 80 85 90 95; do
	cd ~/data/heritability_change
	phen=rsfmri_100x100_OD0.${OD}_12032019;
	for t in {1..4950}; do
		solar run_phen_var_fixed_OD_xcp ${phen} conn${t};
	done;
	mv ${phen} ~/data/tmp/
	cd ~/data/tmp/${phen}
	for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
	python ~/research_code/compile_solar_multivar_results.py ${phen}
done
```

But that will take all night to run. So, let's explore the 7 network model
again, but now with the fixed code. And this SOLAR run shouldn't take as long:

```bash
# sinteractive so we don't screwp up the runs above
for OD in 80 85 90 95; do
	for suf in 'median' 'mean'; do
		for p in '' 'All'; do
			cd ~/data/heritability_change
			phen=rsfmri_7by7from100_7nets_OD0.${OD}_${suf}${p}_12032019;
			for t in "conn_DorsAttnTODorsAttn" "conn_DorsAttnTOSalVentAttn" \                       "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
				"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
				"conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
				"conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
				"conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
				"conn_DefaultTODefault"; do
					solar run_phen_var_OD_xcp ${phen} ${t};
			done;
			mv ${phen} ~/data/tmp/
			cd ~/data/tmp/${phen}
			for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
			python ~/research_code/compile_solar_multivar_results.py ${phen}
		done;
	done;
done
```

# 2019-12-04 08:57:36

Let's start analyzing the 7nets results, because they're easier to parse out:

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_7nets.csv;
for f in `ls polygen_results_*7nets*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_7nets.csv;
	done
done
```

Nothing there... very little consistency in significant results. Will need to go
to the 3 or 4 network cases, if going this route. Let's parse the 100x100 cases
first.

Right off the bat I can see that the number of p=.05 will swamp any sort of FDR
correction, and Bonferroni/MEff don't look too promising either. I'll need to do
some trimming in the connections to see if anything else survives.

```r
phen = 'rsfmri_100x100_OD0.80_12032019'

df = read.csv(sprintf('~/data/tmp/polygen_results_%s.csv', phen))
nrois=100
fname = sprintf('~/research_code/fmri/Schaefer2018_%dParcels_7Networks_order.txt',
				nrois)
nets = read.table(fname)
all_net_names = sapply(as.character(unique(nets[,2])),
					   function(y) strsplit(x=y, split='_')[[1]][3])
net_names = unique(all_net_names)
nnets = length(net_names)
cat('Creating connection map...\n')
nverts = nrow(nets)
cnt = 1
conn_map = c()
for (i in 1:(nverts-1)) {
	for (j in (i+1):nverts) {
		conn = sprintf('conn%d', cnt)
		conn_map = rbind(conn_map, c(conn, all_net_names[i], all_net_names[j]))
		cnt = cnt + 1
	}
}

# use_only = c('Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default')
# use_only = c('DorsAttn', 'SalVentAttn', 'Cont', 'Default')
use_only = c('DorsAttn', 'SalVentAttn', 'Cont', 'Default', 'Limbic')

keep_me = c()
for (i in 1:length(use_only)) {
	for (j in i:length(use_only)) {
		idx = (conn_map[,2]==use_only[i] & conn_map[,3]==use_only[j]) |
			(conn_map[,3]==use_only[i] & conn_map[,2]==use_only[j])
		keep_me = c(keep_me, conn_map[idx, 1])
	}
}

df2 = df[df$phen %in% keep_me, ]
p2 = p.adjust(df2$h_pval, method='fdr')

fname = sprintf('~/data/heritability_change/%s.csv', phen)
data = read.csv(fname)
cc = cor(data[, keep_me])
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)

print(sprintf('Surviving at FDR at .05 = %d', length(which(p2<.05))))
print(sprintf('Surviving at FDR at .1 = %d', length(which(p2<.1))))
cat(sprintf('Galwey Meff = %.2f\n', meff))
print(sprintf('Surviving at Meff = %d', length(which(df2$h_pval<(.05/meff)))))
```

Using 4 nets, I have 1 connection surviving Meff at .80, 3 at .85, 2 at .9 (the
same 2 survive at FDR .1 too), and 3 at .95. 

Using 5 nets (including limbic), I still get 3 for Meff at .95, 2 for .9, 3 for
.85, and 1 for .8.

Before I analyze those connections, let's try to make it run with positive only.
We can impute or just drop the variables. Let's try dropping first. They will
take a while to run in SOLAr, so we can get it going on the background. There's
also a chance that screening would help this, even though it would take longer
to run... let's try it out too.

```bash
# sinteractive
OD=80;
mkdir /lscratch/$SLURM_JOBID/${OD}s
cd /lscratch/$SLURM_JOBID/${OD}s
phen=rsfmri_100x100_posOnly_OD0.${OD}_12042019;
cp ~/data/heritability_change/pedigree.csv .
cp ~/data/heritability_change/procs.tcl .
cp ~/data/heritability_change/${phen}.csv .
for t in {1..4950}; do
	solar run_phen_var_OD_xcp ${phen} conn${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

```bash
# sinteractive
OD=80;
mkdir /lscratch/$SLURM_JOBID/${OD}s_all
cd /lscratch/$SLURM_JOBID/${OD}s_all
phen=rsfmri_100x100_OD0.${OD}_12032019;
cp ~/data/heritability_change/pedigree.csv .
cp ~/data/heritability_change/procs.tcl .
cp ~/data/heritability_change/${phen}.csv .
for t in {1..4950}; do
	solar run_phen_var_OD_xcp ${phen} conn${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

**Just remember to ignore any connections where we had to impute more than 15% of
the data!!!** For now I didn't do that in the script, but when computing the
results I'll need to do it.

## No ventral

Same as the 3 nets variation I did before, but now Philip suggested I should try
DAN, DMN, and Cog. Then we can add the other ones, keeping the same subject set,
just for robustness. Well, let's see if these results are good to begin with:

```bash
# sinteractive so we don't screwp up the runs above
for OD in 80 85 90 95; do
	for suf in 'median' 'mean'; do
		for p in '' 'All'; do
			cd ~/data/heritability_change
			phen=rsfmri_7by7from100_3netsCog_OD0.${OD}_${suf}${p}_12042019;
			for t in "conn_DorsAttnTODorsAttn" "conn_DorsAttnTOCont" \
					 "conn_DorsAttnTODefault" "conn_ContTOCont" \
					 "conn_ContTODefault" "conn_DefaultTODefault"; do
					solar run_phen_var_OD_xcp ${phen} ${t};
			done;
			mv ${phen} ~/data/tmp/
			cd ~/data/tmp/${phen}
			for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
			python ~/research_code/compile_solar_multivar_results.py ${phen}
		done;
	done;
done
```

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_3netsCog.csv;
for f in `ls polygen_results_*3netsCog*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_3netsCog.csv;
	done
done
```

![](images/2019-12-04-11-21-39.png)

Nothing earth-shattering there, but it's also just 6 connections, so let's check
them individually.

meanAll is fine for Bonferroni and FDR, and so is median. medianAll as well. But
the mean results don't survive FDR .05 (they do at .1 though). FDR on 6
comparison would be overkill? Not sure if the method has limitations... It is
somewhat funky that our top 5 variations all came up with different network
connections being significant... not very robust :(

So, our initial 3 net result is still better...

## QC on 100x100

How about doing QC on the 100x100 matrix, running 7x7 matrices, and selecting
  which ones to look at later?

```bash
Rscript ~/research_code/fmri/make_outlier_detection_slopes_rawQC.R .9
```

```bash
# sinteractive so we don't screwp up the runs above
for OD in 80 85 90 95; do
	for suf in 'median' 'mean'; do
		for p in '' 'All'; do
			cd ~/data/heritability_change
			phen=rsfmri_7by7from100_rawQC_OD0.${OD}_${suf}${p}_12042019;
			for t in "conn_DorsAttnTODorsAttn" "conn_DorsAttnTOSalVentAttn" \                       "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
				"conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
				"conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
				"conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
				"conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
				"conn_DefaultTODefault"; do
					solar run_phen_var_OD_xcp ${phen} ${t};
			done;
			mv ${phen} ~/data/tmp/
			cd ~/data/tmp/${phen}
			for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
			python ~/research_code/compile_solar_multivar_results.py ${phen}
		done;
	done;
done
```

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_rawQC.csv;
for f in `ls polygen_results_*rawQC*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_rawQC.csv;
	done
done
```

## evaluating negative imputation

Let's evaluate our results after imputation, but first removing any connections
for which we had to impute more than 15% of the data.

```r
phen = 'rsfmri_100x100_posOnly_OD0.80_12042019'

df = read.csv(sprintf('~/data/tmp/polygen_results_%s.csv', phen))
nrois=100
fname = sprintf('~/research_code/fmri/Schaefer2018_%dParcels_7Networks_order.txt',
				nrois)
nets = read.table(fname)
all_net_names = sapply(as.character(unique(nets[,2])),
					   function(y) strsplit(x=y, split='_')[[1]][3])
net_names = unique(all_net_names)
nnets = length(net_names)
cat('Creating connection map...\n')
nverts = nrow(nets)
cnt = 1
conn_map = c()
for (i in 1:(nverts-1)) {
	for (j in (i+1):nverts) {
		conn = sprintf('conn%d', cnt)
		conn_map = rbind(conn_map, c(conn, all_net_names[i], all_net_names[j]))
		cnt = cnt + 1
	}
}

fname = sprintf('~/data/heritability_change/%s.csv', phen)
data = read.csv(fname)

# use_only = c('Vis', 'SomMot', 'DorsAttn', 'SalVentAttn', 'Limbic', 'Cont', 'Default')
# use_only = c('DorsAttn', 'SalVentAttn', 'Cont', 'Default')
use_only = c('DorsAttn', 'SalVentAttn', 'Cont', 'Default', 'Limbic')

keep_me = c()
for (i in 1:length(use_only)) {
	for (j in i:length(use_only)) {
		idx = (conn_map[,2]==use_only[i] & conn_map[,3]==use_only[j]) |
			(conn_map[,3]==use_only[i] & conn_map[,2]==use_only[j])
		keep_me = c(keep_me, conn_map[idx, 1])
	}
}

df2 = df[df$phen %in% keep_me, ]
p2 = p.adjust(df2$h_pval, method='fdr')

cc = cor(data[, keep_me])
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)

print(sprintf('Surviving at FDR at .05 = %d', length(which(p2<.05))))
print(sprintf('Surviving at FDR at .1 = %d', length(which(p2<.1))))
cat(sprintf('Galwey Meff = %.2f\n', meff))
print(sprintf('Surviving at Meff = %d', length(which(df2$h_pval<(.05/meff)))))
```

Hum... it turns out that most of the brain has lots of negative connections...
that doesn't look good. Maybe I should go back to degree-centrality analysis?
But that's also based on Pearson correlation, I'd think. Maybe we could try ALF
and reho again?

Or how about number of significant absolute connections?

```bash
Rscript ~/research_code/fmri/make_outlier_detection_slopes_rawSigQC.R .9
```

```bash
# sinteractive so we don't screwp up the runs above
for OD in 80 85 90 95; do
	cd ~/data/heritability_change
	phen=rsfmri_7by7from100_rawQC_OD0.${OD}_${suf}${p}_12042019;
	for t in "conn_DorsAttnTODorsAttn" "conn_DorsAttnTOSalVentAttn" \                       "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
  "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
  "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
  "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
  "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
  "conn_DefaultTODefault"; do
		solar run_phen_var_OD_xcp ${phen} ${t};
	done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
done;
```

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output_rawQC.csv;
for f in `ls polygen_results_*rawQC*.csv`; do
	# skip header
	for line in `tail -n +2 $f`; do
		echo $f,$line >> output_rawQC.csv;
	done
done
```



# TODO
* try removing ventral but not cognitive and then just add the other networks
  without changing subjects for robustness
* how about doing QC on the 100x100 matrix, running 7x7 matrices, and selecting
  which ones to look at later?
* review all numbers in the paper (and figures!)
* robustness analysis (and everything else in note 044)
* remember that fmri is .9 and dti is potentially some different OD threshold!
* recompute table that has everyone