# 2018-09-18 11:31:23

The idea here is to start from the list of everybody who has longitudinal
clinical data, and select the earliest good scan for each person. There was a
lot of manual work there. As usual, I split it into different sets of minutes. 

If we don't take when the clinical data was taken, we have 331 (of the 380
possible) with at least one rsFMRI, 289 with at least 3 min, 263 with at least
4min, and 225 with at least 5min. Those numbers always go down when we add in
the clinical data. Let's see:

```r
> library(gdata)
> df = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
> df_base = df[df$baseline_any,]
> dim(df_base)
[1] 331  12
> clin = read.csv('~/data/baseline_prediction/long_clin_0913.csv')
> m = merge(df_base,clin, by.x='Medical.Record...MRN...Subjects', by.y='MRN')
> dim(m)
[1] 331  24
> scan_date = as.Date(as.character(m$record.date.collected...Scan), format='%m/%d/%Y')
> doa_date = as.Date(as.character(m$last_DOA), format='%m/%d/%y')
> d = doa_date - scan_date
> sum(d > 30*12)
[1] 273
> df_base = df[df$baseline_3min,]
> dim(df_base)
[1] 289  12
> m = merge(df_base,clin, by.x='Medical.Record...MRN...Subjects', by.y='MRN')
> scan_date = as.Date(as.character(m$record.date.collected...Scan), format='%m/%d/%Y')
> doa_date = as.Date(as.character(m$last_DOA), format='%m/%d/%y')
> sum((doa_date - scan_date) > 30*12)
[1] 218
> df_base = df[df$baseline_4min,]
> dim(df_base)
[1] 263  12
> m = merge(df_base,clin, by.x='Medical.Record...MRN...Subjects', by.y='MRN')
> scan_date = as.Date(as.character(m$record.date.collected...Scan), format='%m/%d/%Y')
> doa_date = as.Date(as.character(m$last_DOA), format='%m/%d/%y')
> sum((doa_date - scan_date) > 30*12)
[1] 172
> df_base = df[df$baseline_min,]
> df_base = df[df$baseline_5min,]
> m = merge(df_base,clin, by.x='Medical.Record...MRN...Subjects', by.y='MRN')
> scan_date = as.Date(as.character(m$record.date.collected...Scan), format='%m/%d/%Y')
> doa_date = as.Date(as.character(m$last_DOA), format='%m/%d/%y')
> sum((doa_date - scan_date) > 30*12)
[1] 119
```

The threshold to best match structural and DTI would be having any data, but
that seems a bit too lenient. Let's do with at least 3min of data. If we end up
re-processing everything using a new motion threshold based on Marine's
experiments, we can try something else.

```r
> ids = sapply(m[(doa_date - scan_date) > 30*12,]$Mask.ID...Scan, function(x) sprintf("%04d", x))
> length(ids)
[1] 218
> write.table(ids, file='~/data/baseline_prediction/rsfmri_3minWithClinical.tsv', row.names=F, col.names=F, quote=F)
```

Also, let's make sure they all have proper alignment... so, from that list, I
had to replace 1008 by 1824, as it had QC of 4. Then there were 3 other IDs with
QC of 3 (787, 891, 1941) because of missing parts of the brain or bad alignment.
So, I removed them manually for a final number of 215 scans.

As usual, let's create trimmed and non_trimmed versions of the files, and we
should also create correlation files for both Freesurfer parcellations, just to
be safe we're testing everything we can.

```bash
# caterpie
for m in `cat ~/tmp/rsfmri_3minWithClinical.tsv`; do 
    echo $m;
    scp -q /mnt/shaw/freesurfer5.3_subjects/${m}/SUMA/aparc+aseg_REN_all.nii.gz helix:/scratch/sudregp/rsfmri/${m}_aparc.nii.gz;
    scp -q /mnt/shaw/freesurfer5.3_subjects/${m}/SUMA/aparc.a2009s+aseg_REN_all.nii.gz helix:/scratch/sudregp/rsfmri/${m}_aparc.a2009s.nii.gz;
    scp -q /mnt/shaw/data_by_maskID/${m}/afni/${m}.rest.subjectSpace.results/errts.${m}.fanaticor+orig.* helix:/scratch/sudregp/rsfmri/;
done
```

Then, in biowulf (note that I need different job names!): 

```bash
for s in `cat ~/tmp/rsfmri_3minWithClinical.tsv`; do
sbatch -J rois --cpus-per-task 2 \
    --partition quick --gres=lscratch:10 <<EOF
#!/bin/bash
module load afni; bash ~/research_code/fmri/extract_roi_swarm.sh $s
EOF
done
```

Now we collect it:

```bash
python ~/research_code/fmri/collect_atlas_connectivity_subjectSpace.py
mv ~/tmp/aparc*csv ~/data/baseline_prediction/ 
```

# 2018-09-19 09:52:56

And finally change the variable names to conform to the algorithms:

```r
library(gdata)
df = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
df = df[,1:2]
colnames(df) = c('MRN', 'mask.id')
for (i in c('', '.a2009s')) {
    for (j in c('', '_trimmed')) {
        print(sprintf('%s, %s', i, j))
        fname = sprintf('~/data/baseline_prediction/aparc%s%s.csv', i, j)
        data = read.csv(fname)
        data = merge(df, data, by.x='mask.id', by.y='maskid')
        data$mask.id = data$MRN
        var_names = grepl('TO', colnames(data))
        cnames = sapply(colnames(data)[var_names], function(x) sprintf('v_%s', x))
        colnames(data)[var_names] = cnames
        fname = sprintf('~/data/baseline_prediction/aparc%s%s_n215_09182018.RData.gz', i, j)
        save(data, file=fname, compress=T)
    }
}

clin = read.csv('~/data/baseline_prediction/long_clin_0913.csv')
colnames(clin)[1]='mask.id'
write.csv(clin, file='~/data/baseline_prediction/rsfmri_gf_09182018.csv', row.names=F)
```

# 2018-09-28 11:25:35

Added movement, sex, and age variables to the rsfmri file
(rsfmri_09282018.xlsx).

```bash
# caterpie
while read s; do
    mydir=/mnt/shaw/data_by_maskID/${s}/afni/${s}.rest.subjectSpace.results;
    if [ -d $mydir ]; then
        echo $s;
        used=`1d_tool.py -infile ${mydir}//X.xmat.1D -show_rows_cols -verb 0 | cut -d " " -f 1 -`;
        if [ -e ${mydir}//X.nocensor.xmat.1D ]; then
            total=`1d_tool.py -infile ${mydir}//X.nocensor.xmat.1D -show_rows_cols -verb 0 | cut -d " " -f 1 -`;
        else
            total=$used;
        fi;
        1deval -a ${mydir}/motion_${s}_enorm.1D -b ${mydir}/censor_${s}_combined_2.1D -expr 'a*b' > ~/tmp/rm.ec.1D;
        mvmt=`3dTstat -prefix - -nzmean ~/tmp/rm.ec.1D\'`;
        echo ${s},${total},${used},${mvmt} >> ~/tmp/good_trs.csv;
    fi;
done < ~/tmp/rsfmri.txt
```

# 2018-10-25 16:12:30

I had to redo the fMRI data because Aman found an issue that censored TRs were
still in the errts file, so we need to first remove any TRs with time series 0
before computing correlations. I changed the compute file today, and now should
have more accurate correlation matrices.

```bash
python ~/research_code/fmri/collect_atlas_connectivity_subjectSpace.py
# renamed old files to _old.csv
mv ~/tmp/aparc*csv ~/data/baseline_prediction/
```

And again change the variable names to conform to the algorithms:

```r
library(gdata)
df = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
df = df[,1:2]
colnames(df) = c('MRN', 'mask.id')
for (i in c('', '.a2009s')) {
    for (j in c('', '_trimmed')) {
        print(sprintf('%s, %s', i, j))
        fname = sprintf('~/data/baseline_prediction/aparc%s%s.csv', i, j)
        data = read.csv(fname)
        data = merge(df, data, by.x='mask.id', by.y='maskid')
        data$mask.id = data$MRN
        var_names = grepl('TO', colnames(data))
        cnames = sapply(colnames(data)[var_names], function(x) sprintf('v_%s', x))
        colnames(data)[var_names] = cnames
        fname = sprintf('~/data/baseline_prediction/aparc%s%s_n215_10252018.RData.gz', i, j)
        save(data, file=fname, compress=T)
    }
}
```

# 2018-11-15 10:18:10

I was taking a look at the number of variables, and aparc has 3919 and aparc.a2009s has about 14K, so neither one is too much bigger than DTI. And that's even before I clean up the datasets by removing meaningless ROIs, so I don't think the number of variables is an issue. Let's clean them up a bit, and then maybe play a bit with within subject normalization. We can explore ICA and graph properties later.

Let's redo the fMRI datasets after cleaning up a bit of the variables that don't make much sense. We can also make the pcorr datasets. 

So, I started by running fmri/make_all_correlations.R and fmri/make_all_partial_correlations.R. Removing ROIs only matter for partial correlations, and we have already done that. For full correlations, we can just remove them when we construct the prediction datasets.

In other words, for the partial_correlation datasets all we need to do is add in the MRNs and rename variables:

```r
library(gdata)
df = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
df = df[,1:2]
colnames(df) = c('MRN', 'mask.id')
for (i in c('kendall', 'pearson', 'spearman')) {
    for (j in c('', '_trimmed')) {
        print(sprintf('%s, %s', i, j))
        fname = sprintf('~/data/baseline_prediction/rsfmri/partial_weighted_%s%s.csv',
                        i, j)
        data = read.csv(fname)
        colnames(data)[1] = 'mask.id'
        data = merge(df, data, by='mask.id')
        var_names = grepl('TO', colnames(data))
        cnames = sapply(colnames(data)[var_names], function(x) sprintf('v_%s', x))
        colnames(data)[var_names] = cnames
        fname = sprintf('~/data/baseline_prediction/aparc_pcorr_%s%s_n215_11152018.RData.gz',
            i, j)
        save(data, file=fname, compress=T)
    }
}
```

There are a couple subjects that have 100 or so of the 2280 partial connections as NA. I don't think that will screw up the algorithms too much... let's see.

For the full correlation it's a bit trickier, because we'll have both parcelations, and will need to remove any garbage ROIs. Let's see:

(Note that I left a few more ROIs in full correlation because I didn't want to remove some of the thalamic ones. I thought it would be necessary for partial correlations to remove the number of ROIs, but not as much for full correlations)

```r
rm_rois = c('CSF', 'Ventricle', 'Pallidum', 'Brain.Stem',
            'Accumbens.area', 'VentralDC', 'vessel', 'Cerebral',
            'choroid', 'Lat.Vent', 'White.Matter', 'hypointensities',
            '^CC', 'nknown', 'Chiasm', 'Cerebellum.Cortex', 'undetermined')
library(gdata)
df = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
df = df[,1:2]
colnames(df) = c('MRN', 'mask.id')
for (i in c('kendall', 'pearson', 'spearman')) {
    for (j in c('', '_trimmed')) {
        for (p in c('', '.a2009s')) {
            print(sprintf('%s, %s, %s', i, j, p))
            fname = sprintf('~/data/baseline_prediction/rsfmri/weighted_aparc%s_%s%s.csv',
                            p, i, j)
            data = read.csv(fname)
            colnames(data)[1] = 'mask.id'
            rm_me = c()
            for (r in rm_rois) {
                rm_me = c(rm_me, which(grepl(r, colnames(data)))) }
            data = data[, -unique(rm_me)]
            data = merge(df, data, by='mask.id')
            var_names = grepl('TO', colnames(data))
            cnames = sapply(colnames(data)[var_names], function(x) sprintf('v_%s', x))
            colnames(data)[var_names] = cnames
            fname = sprintf('~/data/baseline_prediction/aparc%s_corr_%s%s_n215_11152018.RData.gz',
                p, i, j)
            save(data, file=fname, compress=T)
        }
    }
}
```

So, in this clean version I end up with 12923 variables in aparc.a2009s and 3118
in aparc. For pcorr we have 2278. 

# 2018-11-16 12:45:02

Just to note that I tried GIFT (Calhoun's tool) with the ROI data, and also just
the NIFTI, and neither worked. The NIFTIs ran out of memory even with 120Gb for
one subject... not sure if it's worth going higher. The ROI data didn't work at
all, even after formatted into NIFTI.

Now I'm trying to run ICASSO on the time concatenated (after zscoring each time
series) of the ROIs. I get that by running:

```bash
Rscript ~/research_code/fmri/make_roi_niftis.R
sed -i -e "s/NA//g" tcat_aparc.csv
```

The resulting CSV is TRs by ROIs, so when loading it in Matlab we need to
transpose it:

```matlab
restoredefaultpath()
addpath('/data/NCR_SBRB/software/FastICA_25/')
addpath('/data/NCR_SBRB/software/icasso122/')
addpath('/data/NCR_SBRB/')
Ydd = dlmread(['/data/sudregp/baseline_prediction/rsfmri/roi_niftis/tcat_aparc.a2009s.csv'], ',', 1, 0);

sR=icassoEst('both', Ydd', 1000, 'lastEig', 15, 'g', 'pow3', 'approach', 'defl');
sR=icassoExp(sR);
[iq,A,W,S]=icassoResult(sR);
save(['/data/sudregp/baseline_prediction/rsfmri/roi_niftis/aparc.a2009s_1Kperms_15ics.mat'],'A','S','W','iq','sR','-v7.3')
```

NEED TO RUN OTHER FILES...

<!-- Of course we could choose just ICs above a certain iQ score. But for now, since
we're just operating with 15 ICs, let's use all of them and let ML decide which
ones are helpful or not.

Finally, we compute expression scores:

```python
import tables
import numpy as np
import statsmodels.formula.api as smf
from scipy import io
import pandas as pd

suffixes=['aparc', 'aparc.a2009s']
tsuffixes=[s + '_trimmed' for s in suffixes]
for suffix in suffixes + tsuffixes:
     print suffix
     file = tables.open_file('/data/sudregp/baseline_prediction/rsfmri/roi_niftis/%s_1Kperms_15ics.mat' % suffix)
     S = file.root.S[:]
     Ydd = t(np.genfromtxt('/data/sudregp/baseline_prediction/rsfmri/roi_niftis/tcat_%s.csv' % suffix, delimiter=',', skip_header=1))
     scores = []
     for s in range(Ydd.shape[0]):
          est = smf.OLS(Ydd[s,: ],S).fit()
          scores.append(est.params)
     exp_scores = np.array(scores)
     df = pd.DataFrame(exp_scores, columns=['ic%02d' % i for i in range(1, 16)])
     df.to_csv('/data/sudregp/baseline_prediction/rsfmri/roi_niftis/exp_scores_%s_1Kperms_15ics.csv' % suffix, index=False)
``` -->

# 2018-11-19 15:57:19

## Variability

I'm having a quick run with variability, just to see if the absolute difference
between connections (first and second half of TRs) can tells us anything. So, we
start by running:

```bash
ipython ~/research_code/fmri/collect_atlas_connectivityVariation_subjectSpace.py
```

which creates the trimmed/non-trimmed matrices for the 3 methods and both
parcellations. Then, the next step is to clean it up and create our data
matrices:

```r
rm_rois = c('CSF', 'Ventricle', 'Pallidum', 'Brain.Stem',
            'Accumbens.area', 'VentralDC', 'vessel', 'Cerebral',
            'choroid', 'Lat.Vent', 'White.Matter', 'hypointensities',
            '^CC', 'nknown', 'Chiasm', 'Cerebellum.Cortex', 'undetermined')
library(gdata)
df = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
df = df[,1:2]
colnames(df) = c('MRN', 'mask.id')
for (i in c('kendall', 'pearson', 'spearman')) {
    for (j in c('', '_trimmed')) {
        for (p in c('', '.a2009s')) {
            print(sprintf('%s, %s, %s', i, j, p))
            fname = sprintf('~/data/baseline_prediction/rsfmri/%sAbsDiff_aparc%s%s.csv',
                            i, p, j)
            data = read.csv(fname)
            colnames(data)[2] = 'mask.id'
            rm_me = c()
            for (r in rm_rois) {
                rm_me = c(rm_me, which(grepl(r, colnames(data)))) }
            data = data[, -unique(rm_me)]
            data = merge(df, data, by='mask.id')
            var_names = grepl('TO', colnames(data))
            cnames = sapply(colnames(data)[var_names], function(x) sprintf('v_%s', x))
            colnames(data)[var_names] = cnames
            fname = sprintf('~/data/baseline_prediction/aparc%s_absDiff_%s%s_n215_11202018.RData.gz',
                p, i, j)
            save(data, file=fname, compress=T)
        }
    }
}
```