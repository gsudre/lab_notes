# 2019-10-01 11:00:23

Let's conduct some extra analysis for the heritability paper.

## Slopes and age

Similar to the Jalmar Teeuw Neuroimage 2019 paper, can we show that FC within
networks increases with age whereas FC between networks decreases with age?

```r
data = read.csv('~/data/heritability_change/rsfmri_7by7from100_5nets_OD0.90_posOnly_median.csv')
idx1 = grepl(colnames(data), pattern='^conn')
idx2 = grepl(colnames(data), pattern='_baseline$')
idx3 = grepl(colnames(data), pattern='SomMot')
idx4 = grepl(colnames(data), pattern='Vis')
mynets = colnames(data[idx1 & !idx2 & !idx3 & !idx4])
res = c()
for (net in mynets) {
    tmp = t.test(data[, net])
    res = rbind(res, c(net, tmp$statistic, tmp$p.value, tmp$estimate))
}
colnames(res) = c('net', 'tstat', 'pval', 'mean')
write.csv(res, file='~/data/heritability_change/inter_intra.csv', quote=F, row.names=F)
```

## Cross-modal genetic correlation

Start by merging the data across modalities. Then, it should just be a matter of
creating the right procedure in SOLAR.

```r
fmri = read.csv('~/data/heritability_change/rsfmri_7by7from100_5nets_OD0.90_posOnly_median.csv')
dti = read.csv('~/data/heritability_change/dti_JHUtracts_ADRDonly_OD0.95.csv')
both = merge(fmri, dti, by='ID', all.x=F, all.y=F)
# 156 participants with both datasets... not bad
write.csv(both, file='~/data/heritability_change/both_fmri_dti.csv', row.names=F, quote=F)
```

```bash
module load solar
cd ~/data/heritability_change/
solar dti_rsfmri_slope_correlation both_fmri_dti
cd gencor_both_fmri_dti
grep -r RhoG */polygenic.out > rhog.txt
grep zero rhog.txt | cut -d "/" -f 1 > nets.txt
grep zero rhog.txt | cut -d "=" -f 2 > ps.txt
paste nets.txt ps.txt > cross_modal_rhog_pvals.txt
```

## Degree centrality

First step is to check if Luke processed all our candidates:

```r
demo = read.csv('~/data/heritability_change/resting_demo_07032019.csv')

# keeping it to kids only to make sure everyone has been processed
demo = demo[demo$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(demo)))

# let's grab QC metrics on everyone
# note that this only works for non-censoring pipelines!
mydir = '/Volumes/Shaw/rsfmri_36P/xcpengine_output_fc-36p_despike/'
qc_data = c()
for (s in demo$Mask.ID) {
    subj = sprintf('sub-%04d', s)
    # if it processed all the way
    std_fname = sprintf('%s/%s/norm/%s_std.nii.gz', mydir, subj, subj)
    if (file.exists(std_fname)) {
        subj_data = read.csv(sprintf('%s/%s/%s_quality.csv', mydir, subj, subj))
        qc_data = rbind(qc_data, subj_data)
    }
}
```

So, we start with 764 scans in qc_data.

```r
missing = c()
for (s in qc_data$id0) {
    subj = gsub(x=s, pattern='-', replacement='')
    dc_fname = sprintf('/Volumes/Shaw/Functional_derivatives_RS36_despiked/smoothed/sdc_weighted/sDegreeCentrality_PositiveWeightedSumBrainMap_%s.nii', subj)
    if (!file.exists(dc_fname)) {
        missing = c(missing, subj)
    }
}
```

There are 359 scans missing, so I asked Luke for them. They're likely scans with
bad motion or not longitudinal, but I need them in my analysis to figure out
outliers. Luke's working on it.

Actually, it'd take a long time to calculate it using the software he used. It's
quite finicky. I'd do it myself for one mask id, and then compare to his
results. Basically, for every voxel in the gray mask, calculate Pearson
correlation to everything else, threshold at .25 (positive), compute degree centrality
on the surviving edges. Spit it out to a .nii and then smooth out the results.

# 2019-10-04 16:59:22

It's getting too time consuming to re-code this. Let me see if removing some of
the bad scans we can reduce the list we'd need Luke to run. I'll go with 95th
pctile just to get as broad as a brush we can:

```r
demo = read.csv('~/data/heritability_change/resting_demo_07032019.csv')

# keeping it to kids only to make sure everyone has been processed
demo = demo[demo$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(demo)))

# let's grab QC metrics on everyone
# note that this only works for non-censoring pipelines!
mydir = '/Volumes/Shaw/rsfmri_36P/xcpengine_output_fc-36p_despike/'
qc_data = c()
for (s in demo$Mask.ID) {
    subj = sprintf('sub-%04d', s)
    # if it processed all the way
    std_fname = sprintf('%s/%s/norm/%s_std.nii.gz', mydir, subj, subj)
    if (file.exists(std_fname)) {
        subj_data = read.csv(sprintf('%s/%s/%s_quality.csv', mydir, subj, subj))
        qc_data = rbind(qc_data, subj_data)
    }
}
# have some higly correlated qc variables, so let's remove the worse offenders (anything above abs(.8))
qc_vars = c('normCoverage', 'meanDV', 'pctSpikesDV',
            'motionDVCorrInit',
            'motionDVCorrFinal', "pctSpikesRMS", "relMeanRMSMotion")

qtile=.95
library(solitude)
iso <- isolationForest$new()
iso$fit(qc_data[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]

library(dbscan)
scores_lof = lof(qc_data[, qc_vars], k = round(.5 * nrow(qc_data)))

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)

idx = scores_lof < thresh_lof & scores_if < thresh_if
```

So, we start with 696 scans in qc_data.

```r
missing = c()
for (s in qc_data[idx,]$id0) {
    dc_fname = sprintf('/Volumes/Shaw/Functional_derivatives_RS36_despiked/smoothed/sdc_weighted/sDegreeCentrality_PositiveWeightedSumBrainMap_%s.nii', s)
    if (!file.exists(dc_fname)) {
        missing = c(missing, s)
    }
}
```

Still missing 298... Luke will run them for me. Let's get the code ready to run.
First, extract the values:

```bash
# desktop
outdir=~/data/heritability_change/sdc_weighted
mkdir -p $outdir;
cd /Volumes/Shaw/Functional_derivatives_RS36_despiked/smoothed/sdc_weighted/;
for f in `ls *nii`; do
    fout=`basename $f .nii`.txt;
    if [ ! -e $outdir/$fout ]; then
        3dmaskdump -mask ~/data/heritability_change/xcp-36p_despike/gray_matter_mask.nii \
            -o ${outdir}/${fout} $f;
    fi;
done;
```

I did the above for sdc_weighted and szdc_weighted, even though Luke says he's
using the latter. Then, I created fmri/make_outlier_detection_slopes_DC.R to
deal with that data. Then, when Luke's done processing those scans, the R
function is ready to go.





# TODO
 * Run Luke's metrics on degree centrality
 * Cross-modal genetic correlation on slopes