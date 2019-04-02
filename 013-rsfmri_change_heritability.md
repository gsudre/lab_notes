# 2019-03-27 16:24:43

I was brainstorming with Philip, and the first idea is to just run the
correlation matrices, calculate the slopes, and try to find heritability there.
The question then becomes what are the vertices in the matrix. We could go for
Freesurfer ROIs, but we might end up getting killed by multiple comparisons.
Another option is to derive the 900 spheres, classify them to belong to
different Yeo networks, and compute within and outside network mean
connectivity. I could also play with some ICA on the overall connectivity maps.
Finally, another idea would be to take prototypical templates of different
resting state networks in fMRI, do a spatial regressions, and just calculate
correlations for those time series (NANs for within-network correlation). Like
the Smith et al 10 resting state networks.

Something else I just thought about: we could have voxel to voxel (or sphere to
sphere) connectivity matrices, one per scan, then filter down to only
connections significant (nominally? FDR?) within scans, then compute slope only
between connections significant in both scans. That should give a good amount of
filtering, especially across subjects. Another quantity we could assess is
proportion of connections still stable? Or positive/negative changes in
connections? 

The first step is to check how many datasets we currently have. As usual, we
could use trimmedVSnontrimmed, as well as pearsonVSkendalVSspearman. **These are
all things to try later if the default (Pearson, nontrimmed) doesn't work out.**

For now, while we wait on all the mriqc parameters for our datasets, let's use
simply the number of clean TRs for QC, like we always do.

# 2019-03-29 17:05:31

```r
a = read.csv('~/data/heritability_change/resting_demo.csv')
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
m = merge(a, b, by='Mask.ID', all.y=F)
m = m[!is.na(m$usedTRs_fmri01), ]
# starting with 1370 scans

# restrict based on QC
minutes = 4
idx = m$usedTRs_fmri01 >= (minutes * 60 / 2.5)
m = m[idx,]
# down to 955 scans

keep_me = c()
for (s in unique(m$Medical.Record...MRN)) {
    subj_scans = m[m$Medical.Record...MRN==s, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                 format="%m/%d/%Y")
    if (length(dates) >= 2) {
        sdates = sort(dates)  # not sure why index.return is not working...
        # make sure there is at least 6 months between scans
        next_scan = 2
        while (((sdates[next_scan] - sdates[1]) < 180) && (next_scan < length(sdates))) {
            next_scan = next_scan + 1
        }
        first_scan_age = subj_scans[dates==sdates[1], 'age_at_scan']
        if (((sdates[next_scan] - sdates[1]) >= 180) && (first_scan_age < 26)) {
            idx1 = which(dates == sdates[1])
            keep_me = c(keep_me, which(m$Mask.ID == subj_scans[idx1, 'Mask.ID']))
            idx2 = which(dates == sdates[next_scan])
            keep_me = c(keep_me, which(m$Mask.ID == subj_scans[idx2, 'Mask.ID']))
        }
    }
}
m2 = m[keep_me, ]
# down to 368 scans
```

# 2019-04-01 11:08:43

OK, this is a good starting point, and that might be the number we use for
association. But let crop it a bit to see what are the heritability numbers so
far. We might have to end up going for 3 or 3.5 minutes.

```r
# make sure every family has at least two people
good_nuclear = names(table(m2$Nuclear.ID...FamilyIDs))[table(m2$Nuclear.ID...FamilyIDs) >= 4]
good_extended = names(table(m2$Extended.ID...FamilyIDs))[table(m2$Extended.ID...FamilyIDs) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, m2[which(m2$Nuclear.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
for (f in good_extended) {
    keep_me = c(keep_me, m2[which(m2$Extended.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
keep_me = unique(keep_me)

fam_subjs = c()
for (s in keep_me) {
    fam_subjs = c(fam_subjs, which(m2[, 'Medical.Record...MRN'] == s))
}
m3 = m2[fam_subjs, ]

write.csv(m2, file='~/data/heritability_change/rsfmri_3min_assoc_n462.csv',
          row.names=F)
```

Yeah, we're down to 164 scans (82 subjects) when requiring 4min of good data.
Let me run the code above for 3.5 and 3min, just to see what we get.

* 3.5 min: 414 scans left for association, 190 for heritability (95 subjects)
* 3 min: 462 scans left for association, 228 for heritability (114 subjects)

We might need to exclude some in the future because their relationship is too
far removed, and we will also lose a few numbers due to outliers within
connections (like we did for DTI). But 3min would be the threshold that best
matches the DTI numbers (133), so let's play with that and then we can check if
the results hold with more stringent thresholds.

OK, so here's what I'll try:

* Freesurfer ROIs (somewhat dependent on quality of Freesurfer segmentation)
* Put everything into common space, and use TTDaemon ROIs (dependent on
  transformation quality)
* Use the data from same space above into either voxelwise / spherewise/
  nertworkwise connectivity.

We start by computing the fMRI correlation tables:

```bash
net_dir=/Volumes/Shaw/MR_data_by_maskid/
cd ~/data/heritability_change/fmri_corr_tables
for maskid in `cut -d"," -f 1 ../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    3dresample \
        -master ${net_dir}/${m}/afni/${m}.rest.subjectSpace.results/errts.${m}.fanaticor+orig \
        -prefix ./rois.nii \
        -inset ${net_dir}/../freesurfer5.3_subjects/${m}/SUMA/aparc+aseg_REN_gm.nii.gz \
        -rmode NN -overwrite &&
    3drefit -labeltable ${net_dir}/../freesurfer5.3_subjects/${m}/SUMA/aparc+aseg_REN_all.niml.lt \
        ./rois.nii &&
    3dNetCorr                                       \
        -inset ${net_dir}/${m}/afni/${m}.rest.subjectSpace.results/errts.${m}.fanaticor+orig                    \
        -in_rois  ./rois.nii                       \
        -prefix  ${m}_aparc                           \
        -fish_z;
    rm rois.nii;
done
```

And we might as well do the one for more ROIs:

```bash
net_dir=/Volumes/Shaw/MR_data_by_maskid/
cd ~/data/heritability_change/fmri_corr_tables
for maskid in `cut -d"," -f 1 ../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    3dresample \
        -master ${net_dir}/${m}/afni/${m}.rest.subjectSpace.results/errts.${m}.fanaticor+orig \
        -prefix ./rois2.nii \
        -inset ${net_dir}/../freesurfer5.3_subjects/${m}/SUMA/aparc.a2009s+aseg_REN_gm.nii.gz \
        -rmode NN -overwrite &&
    3drefit -labeltable ${net_dir}/../freesurfer5.3_subjects/${m}/SUMA/aparc.a2009s+aseg_REN_all.niml.lt \
        ./rois2.nii &&
    3dNetCorr                                       \
        -inset ${net_dir}/${m}/afni/${m}.rest.subjectSpace.results/errts.${m}.fanaticor+orig                    \
        -in_rois  ./rois2.nii                       \
        -prefix  ${m}_aparc.a2009s                           \
        -fish_z;
    rm rois2.nii;
done
```

In the meanwhile, let's go ahead and convert all scans for which we don't have
the results in TT space yet:

```bash
net_dir=/Volumes/Shaw/MR_data_by_maskid/
cd ~/data/heritability_change/fmri_same_space/anat
export OMP_NUM_THREADS=4
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    echo $m;
    mri_convert /Volumes/Shaw/freesurfer5.3_subjects/${m}/mri/orig/001.mgz ./${m}.nii &&
    @SSwarper \
        -input ./${m}.nii \
        -base TT_N27_SSW.nii.gz -subid ${m} &&
    rm ${m}.nii;
done;
```

To read each one in R, we'll do something like this:

```r
a = read.table('~/data/heritability_change/fmri_corr_tables/0411_aparc_000.netcc', header=1)
# remove weird integer row
b = a[2:nrow(a),]
# split matrix into first set of rows as Rs, second set as Zs
rs = b[1:ncol(b),]
zs = b[(ncol(b)+1):nrow(b),]
# put correct names in the square matrix
rownames(rs) = colnames(rs)
rownames(zs) = colnames(zs)
```

# 2019-04-02 09:15:41

Continuing this work, let's check which scans didn't have SUMA, generate it,
then re-run the code above:

```bash
net_dir=/Volumes/Shaw/freesurfer5.3_subjects/
cd ~/data/heritability_change/fmri_corr_tables
for maskid in `cut -d"," -f 1 ../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    if [ ! -e ${net_dir}/${m}/SUMA/aparc.a2009s+aseg_REN_all.niml.lt ]; then
        echo $m >> redo.txt
    fi;
done
```

```bash
for m in `cat xaa`; do
    @SUMA_Make_Spec_FS -sid $m -NIFTI -fspath /Volumes/Shaw/freesurfer5.3_subjects/$m;
done
```

And then I re-ran the code above to create the correlation matrices for the
redo.txt subjects.


# check that zeros are censored by 3dNetCorr!
# check scans that don't have SUMA

