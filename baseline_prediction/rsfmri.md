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

Now I just need to wait until it finishes running in the cluster to collect the
results...


