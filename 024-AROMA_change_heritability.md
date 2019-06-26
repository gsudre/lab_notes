# 2019-06-21 10:08:51

Let's give the first try with the AROMA pipelines in predicting heritability of
change. Luke and I ran regular AROMA and also AROMA+GSR for no threshold, then
.25 and .5mm, so we have a total of 6 pipelines to test. I'll start with the
actual connectivity matrix in the Power atlas, and then try other metrics, like
Luke's metrics or even MELODIC later. The data already comes out in the same
space, and it would be a nice parallel to the DTI voxelwise work.

Philip suggested doing the 2 best timepoints for each subject, and keep it the
same across pipelines. This way we don't deal with the issue of age change
across pipelines. To do that, we first need to compile a metric of how good the
scan is. I'll go with the percentage of spikes in the most stringent threshold
(.25) for now.

In other words, for all scans processed, grab the longitudinal ones, remove
anything with people >= 26, and pick the 2 best with at least 6 months between
them.

```r
a = read.csv('~/data/heritability_change/resting_demo_06212019.csv')
# remove adults and subjects with a single scan. This way we make sure everything for this study was processed
a = a[a$age_at_scan < 18, ]
idx = which(table(a$Medical.Record...MRN)>1)
long_subjs = names(table(a$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(a)) {
    if (a[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
a = a[keep_me,]
a = a[a$processed_AROMA == 'TRUE', ]

outliers = c()
# reading quality metric for all scans
for (m in a$Mask.ID) {
    fname = sprintf('/data/NCR_SBRB/tmp/p25/sub-%04d/sub-%04d_quality.csv', m, m)
    qual = read.csv(fname)
    outliers = c(outliers, qual$pctSpikesFD)
}
a$outliers = outliers

# we should also determine whether we're keeping only scans with a certain amount of time...
#
#

# keeping only subjects with two or more scans, at least 6 months in between scans
keep_me = c()
for (s in unique(a$Medical.Record...MRN)) {
    subj_scans = a[a$Medical.Record...MRN==s, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                 format="%m/%d/%Y")
    if (length(dates) >= 2) {
        best_scans = sort(subj_scans$outliers, index.return=T)
        # make sure there is at least 6 months between scans
        next_scan = 2
        while ((abs(dates[best_scans$ix[next_scan]] - dates[best_scans$ix[1]]) < 180) &&
                (next_scan < length(dates))) {
            next_scan = next_scan + 1
        }
        if (abs(dates[best_scans$ix[next_scan]] - dates[best_scans$ix[1]]) > 180) {
            idx1 = best_scans$ix[1]
            keep_me = c(keep_me, which(a$Mask.ID == subj_scans[idx1, 'Mask.ID']))
            idx2 = best_scans$ix[next_scan]
            keep_me = c(keep_me, which(a$Mask.ID == subj_scans[idx2, 'Mask.ID']))
        }
    }
}
a2 = a[keep_me, ]
print(sprintf('From %d to %d scans', nrow(a), nrow(a2)))
```

So, that's the people with 2 scans, but now let's see how many are in the same
families so we can run heritability:

```r
# make sure every family has at least two people
good_nuclear = names(table(a2$Nuclear.ID...FamilyIDs))[table(a2$Nuclear.ID...FamilyIDs) >= 4]
good_extended = names(table(a2$Extended.ID...FamilyIDs))[table(a2$Extended.ID...FamilyIDs) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, a2[which(a2$Nuclear.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
for (f in good_extended) {
    keep_me = c(keep_me, a2[which(a2$Extended.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
keep_me = unique(keep_me)

fam_subjs = c()
for (s in keep_me) {
    fam_subjs = c(fam_subjs, which(a2[, 'Medical.Record...MRN'] == s))
}
a3 = a2[fam_subjs, ]

# write.csv(a3, file='~/data/heritability_change/rsfmri_3min_assoc_n462.csv',
#           row.names=F)
```

OK, so we're down to 326 scans (163 subjects). But it's likely that not all
scans finished properly for a given scrubbing. So, we'll need to remove anyone
that didn't properly finish.

For association, we're at 612 scans (306 kids).

We start by collecting the fMRI correlation tables:

```r
nconn = 34716
data = matrix(nrow=nrow(a2), ncol=nconn)
for (m in 1:nrow(data)) {
    fname = sprintf('/data/NCR_SBRB/tmp/p25/sub-%04d/fcon/power264/sub-%04d_power264_network.txt',
                    a2[m,]$Mask.ID, a2[m,]$Mask.ID)
    if (file.exists(fname)) {
        data[m, ] = read.table(fname)[,1]
    }
}
data = cbind(a2$Mask.ID, data)
na_conns = rowSums(is.na(data))
data = data[na_conns < nconn, ]
colnames(data) = c('Mask.ID', sapply(1:nconn, function(x) sprintf('conn%d', x)))
# merge the data so that we can again only keep subjects that have 2 scans
m = merge(a2, data, by='Mask.ID', all.x=F)
idx = which(table(m$Medical.Record...MRN)>1)
long_subjs = names(table(m$Medical.Record...MRN))[idx]
keep_me = c()
mymrns = m$Medical.Record...MRN
for (i in 1:nrow(m)) {
    if (mymrns[i] %in% long_subjs) {
        keep_me = c(keep_me, i)
    }
}
m = m[keep_me,]
```

But we should also impose time thresholds for all scans, like 3min and 4min.
Let's see what our numbers look like then:

```r
pipelines = c('', '_p5', '_p25', '-GSR', '-GSR_p5', '-GSR_p25')
at_least_mins = c(0, 3, 4)  # needs to have at least these minutes of data

a = read.csv('~/data/heritability_change/resting_demo_06262019.csv')
cat(sprintf('Starting from %d scans\n', nrow(a)))
# remove adults and subjects with a single scan. This way we make sure everything for this study was processed
a = a[a$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(a)))
a = a[a$processed_AROMA == 'TRUE', ]
cat(sprintf('Down to %d to keep only scans that have been processed\n', nrow(a)))
idx = which(table(a$Medical.Record...MRN)>1)
long_subjs = names(table(a$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(a)) {
    if (a[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
a = a[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than 1 scan\n', nrow(a)))
for (p in pipelines) {
    pipe_dir = sprintf('/data/NCR_SBRB/xcpengine_output_AROMA%s/', p)
    cat(sprintf('Reading quality data from %s\n', pipe_dir))
    outliers = c()
    goodness = c()
    # reading quality metric for all scans
    for (m in a$Mask.ID) {
        fname = sprintf('%s/sub-%04d/sub-%04d_quality.csv', pipe_dir, m, m)
        qual = read.csv(fname)
        if (sum(names(qual)=='nVolCensored') == 0) {
            outliers = c(outliers, 0)
        }
        else {
            outliers = c(outliers, qual$nVolCensored)
        }
        # need to use a quality metric that works in all pipelines, regardless of censoring!
        if (sum(names(qual)=='relMeanRMSMotion') == 0) {
            cat(sprintf('WARNING!!! No relMeanRMSMotion for scan %04d!\n', m))
            goodness = c(goodness, 1000)
        }
        else {
            goodness = c(goodness, qual$relMeanRMSMotion)
        }
    }
    a$outliers = outliers
    a$goodness = goodness

    cat('Loading connectivity data...\n')
    nconn = 34716
    data = matrix(nrow=nrow(a), ncol=nconn)
    for (m in 1:nrow(data)) {
        fname = sprintf('%s/sub-%04d/fcon/power264/sub-%04d_power264_network.txt',
                        pipe_dir, a[m,]$Mask.ID, a[m,]$Mask.ID)
        if (file.exists(fname)) {
            data[m, ] = read.table(fname)[,1]
        }
    }
    data = cbind(a$Mask.ID, data)
    # remove scans that are NAs for all connections
    na_conns = rowSums(is.na(data))
    colnames(data) = c('Mask.ID', sapply(1:nconn, function(x) sprintf('conn%d', x)))

    data = data[na_conns < nconn, ]
    # only keep scans with at least some amount of time
    for (min_time in at_least_mins) {
        uncensored_time = (125 - a$outliers) * 2.5 / 60
        aGood = a[uncensored_time > min_time, ]
        cat(sprintf('\tDown to %d scans with good %d minutes\n', nrow(aGood),
                                                                 min_time))

        # merge the data so we can remove subjects with not enough time DOF
        m = merge(aGood, data, by='Mask.ID', all.x=T)
        cat(sprintf('\t\tDown to %d scans with connectivity data\n', nrow(m)))

        # keeping only the two best scans for each subject, at least 6 months apart
        keep_me = c()
        for (s in unique(m$Medical.Record...MRN)) {
            subj_scans = m[m$Medical.Record...MRN==s, ]
            dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                        format="%m/%d/%Y")
            if (length(dates) >= 2) {
                best_scans = sort(subj_scans$goodness, index.return=T)
                # make sure there is at least 6 months between scans
                next_scan = 2
                while ((abs(dates[best_scans$ix[next_scan]] - dates[best_scans$ix[1]]) < 180) &&
                        (next_scan < length(dates))) {
                    next_scan = next_scan + 1
                }
                if (abs(dates[best_scans$ix[next_scan]] - dates[best_scans$ix[1]]) > 180) {
                    idx1 = best_scans$ix[1]
                    keep_me = c(keep_me, which(m$Mask.ID == subj_scans[idx1, 'Mask.ID']))
                    idx2 = best_scans$ix[next_scan]
                    keep_me = c(keep_me, which(m$Mask.ID == subj_scans[idx2, 'Mask.ID']))
                }
            }
        }
        a2Good = m[keep_me, ]
        cat(sprintf('\t\tDown to %d scans only keeping two best ones 6-mo apart\n',
                    nrow(a2Good)))

        good_na_conns = rowSums(is.na(a2Good))
        for (sc in which(good_na_conns > 1500)) {
            cat(sprintf('WARNING!!! Scan %04d has %d uncovered connections (%.2f %%)\n',
                        a2Good[sc, 'Mask.ID'], good_na_conns[sc], good_na_conns[sc]/nconn*100))
        }

        fname = sprintf('~/data/heritability_change/rsfmri_AROMA%s_%dmin_best2scans.csv',
                        p, min_time)
        write.csv(a2Good, file=fname, row.names=F, na='', quote=F)
        # make sure every family has at least two people
        idx = table(a2Good$Nuclear.ID...FamilyIDs) >= 4
        good_nuclear = names(table(a2Good$Nuclear.ID...FamilyIDs))[idx]
        idx = table(a2Good$Extended.ID...FamilyIDs) >= 4
        good_extended = names(table(a2Good$Extended.ID...FamilyIDs))[idx]
        keep_me = c()
        for (f in good_nuclear) {
            keep_me = c(keep_me, a2Good[which(a2Good$Nuclear.ID...FamilyIDs == f),
                                    'Medical.Record...MRN'])
        }
        for (f in good_extended) {
            keep_me = c(keep_me, a2Good[which(a2Good$Extended.ID...FamilyIDs == f),
                                    'Medical.Record...MRN'])
        }
        keep_me = unique(keep_me)

        fam_subjs = c()
        for (s in keep_me) {
            fam_subjs = c(fam_subjs, which(a2Good[, 'Medical.Record...MRN'] == s))
        }
        a2GoodFam = a2Good[fam_subjs, ]
        cat(sprintf('\t\tDown to %d scans only keeping families\n',
                    nrow(a2GoodFam)))
        fname = sprintf('~/data/heritability_change/rsfmri_AROMA%s_%dmin_best2scansFams.csv',
                        p, min_time)
        write.csv(a2GoodFam, file=fname, row.names=F, na='', quote=F)
    }
}
```

# 2019-06-26 13:48:45

Let's compute the deltas and start running some heritability stuff, at least
with the non-scrubbed data.

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
m = read.csv('~/data/heritability_change/rsfmri_AROMA_0min_best2scans.csv')
df_var_names = colnames(m)[!grepl(colnames(m), pattern="conn")]
clin = read.csv('~/data/heritability_change/clinical_06262019.csv')
df = mergeOnClosestDate(m[, df_var_names], clin, unique(m$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
brain_var_names = colnames(m)[grepl(colnames(m), pattern="conn")]
df2 = merge(df, m[, c('Mask.ID', brain_var_names)], by='Mask.ID', all.x=F)

# make sure we still have two scans for everyone
rm_subjs = names(which(table(df2$Medical.Record...MRN)<2))
rm_me = df2$Medical.Record...MRN %in% rm_subjs
df2 = df2[!rm_me, ]

mres = df2
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    phen_cols = c(brain_var_names, 'SX_inatt', 'SX_HI')
    y = mres[idx[2], phen_cols] - mres[idx[1], phen_cols]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)

    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], 'SX_inatt'])
    row = c(row, mres[idx[base_DOA], 'SX_HI'])
    # DX1 is DSMV definition, DX2 will make SX >=4 as ADHD
    if (mres[idx[base_DOA], 'age_at_scan'] < 16) {
        if ((row[length(row)] >= 6) || (row[length(row)-1] >= 6)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    } else {
        if ((row[length(row)] >= 5) || (row[length(row)-1] >= 5)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    }
    if ((row[length(row)] >= 4) || (row[length(row)-1] >= 4)) {
        DX2 = 'ADHD'
    } else {
        DX2 = 'NV'
    }
    row = c(row, DX)
    row = c(row, DX2)
    res = rbind(res, row)
    print(nrow(res))
}
colnames(res) = c('ID', 'sex', brain_var_names, c('SX_inatt', 'SX_HI',
                                              'inatt_baseline',
                                              'HI_baseline', 'DX', 'DX2'))
# we only open this in R, so it's OK to be RData to load faster
save(res, file='~/data/heritability_change/rsfmri_AROMA_0min_best2scansSlopes_n306_06262019.RData')

# and remove outliers
res_clean = res
for (t in brain_var_names) {
    mydata = as.numeric(res_clean[, t])
    # identifying outliers
    ul = mean(mydata) + 3 * sd(mydata)
    ll = mean(mydata) - 3 * sd(mydata)
    bad_subjs = c(which(mydata < ll), which(mydata > ul))

    # remove within-variable outliers
    res_clean[bad_subjs, t] = NA
}
save(res_clean, file='~/data/heritability_change/rsfmri_AROMA_0min_best2scansSlopesClean_n306_06262019.RData')

# and make sure every family has at least two people
good_nuclear = names(table(m$Nuclear.ID...FamilyIDs))[table(m$Nuclear.ID...FamilyIDs) >= 4]
good_extended = names(table(m$Extended.ID...FamilyIDs))[table(m$Extended.ID...FamilyIDs) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, m[which(m$Nuclear.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
for (f in good_extended) {
    keep_me = c(keep_me, m[which(m$Extended.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
keep_me = unique(keep_me)

fam_subjs = c()
for (s in keep_me) {
    fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
}
res2 = res[fam_subjs, ]
res2_clean = res_clean[fam_subjs, ]

write.csv(res2, file='~/data/heritability_change/rsfmri_AROMA_0min_best2scansFamsSlopes_n163_06262019.csv', row.names=F, na='', quote=F)
write.csv(res2_clean, file='~/data/heritability_change/rsfmri_AROMA_0min_best2scansFamsSlopesClean_n163_06262019.csv', row.names=F, na='', quote=F)

# just need to run this once...
write.table(brain_var_names, file='~/data/heritability_change/power264_conns.txt',
            col.names=F, row.names=F, quote=F)
```

Then, I need to run the same thing for GSR and the other pipelines...

Finally, we do some SOLAR analysis just to see what's going on.

```bash
# bw interactive
module load solar
bash ~/research_code/run_solar_parallel.sh \
    rsfmri_AROMA_0min_best2scansFamsSlopesClean_n163_06262019 \
    ~/data/heritability_change/power264_conns.txt
```

# TODO
 * check that we're not using same DOA for two different scans!















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

The next step was to run a script to collect the 3dNetCorr matrices, based on
the snippet above: 

```r
source('~/research_code/fmri/collect_3dNetCorr_grids.R')
```

Note that I the script above generates matrices for aparc and aparc.a2009s. To
get Zs instead of Rs, since it's just a atanh() convertion, I'll do it within R
so I don't have to recode the function.

Now we follow a similar recipe as what we did in DTI: regress out movement per
scan, then calculate the slopes, and finally remove outliers using NAs, before
dumping it into SOLAR.

I did a quick check and using zeros or not in 3dNetCorr is not making much
difference. Basically, I output the time series and run their correlation with
and without the censored TRs. The differences (besides the significance of each
R, as we reduce the number of observations) was negligible.

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
m2 = read.csv('~/data/heritability_change/rsfmri_3min_assoc_n462.csv')
clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
b = read.csv('~/data/heritability_change/fmri_corr_tables/pearson_3min_n462_aparc.csv')
var_names = colnames(b)[2:ncol(b)]
b[, var_names] = atanh(b[, var_names])
df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id')

library(MASS)
mres = df2
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
for (t in var_names) {
    fm_str = sprintf('%s ~', t)
    fm_str = paste(fm_str,
                   'enormGoodTRs_fmri01 + I(enormGoodTRs_fmri01^2)')
    res.lm <- lm(as.formula(fm_str), data = mres, na.action=na.exclude)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    mres[, t] = residuals(step)
}
res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    for (t in var_names) {
        if (sum(is.na(mres[idx, t])) > 0) {
            # if any of the variables is NA, make the slope NA
            row = c(row, NA)
        } else {
            fm_str = sprintf('%s ~ age_at_scan', t)
            fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
            row = c(row, coefficients(fit)[2])
        }
    }
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], 'SX_inatt'])
    row = c(row, mres[idx[base_DOA], 'SX_HI'])
    # DX1 is DSMV definition, DX2 will make SX >=4 as ADHD
    if (mres[idx[base_DOA], 'age_at_scan'] < 16) {
        if ((row[length(row)] >= 6) || (row[length(row)-1] >= 6)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    } else {
        if ((row[length(row)] >= 5) || (row[length(row)-1] >= 5)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    }
    if ((row[length(row)] >= 4) || (row[length(row)-1] >= 4)) {
        DX2 = 'ADHD'
    } else {
        DX2 = 'NV'
    }
    row = c(row, DX)
    row = c(row, DX2)
    res = rbind(res, row)
    print(dim(res)[1])
}
colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI',
                                              'inatt_baseline',
                                              'HI_baseline', 'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_3min_assoc_n231_slopes.csv',
          row.names=F, na='', quote=F)

res_clean = res
# and remove outliers
for (t in var_names) {
    mydata = as.numeric(res[, t])
    # identifying outliers
    ul = mean(mydata) + 3 * sd(mydata)
    ll = mean(mydata) - 3 * sd(mydata)
    bad_subjs = c(which(mydata < ll), which(mydata > ul))

    # remove within-tract outliers
    res_clean[bad_subjs, t] = NA
}
write.csv(res_clean, file='~/data/heritability_change/rsfmri_3min_assoc_n231_slopesClean.csv',
          row.names=F, na='', quote=F)

# and make sure every family has at least two people
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
    fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
}
res2 = res[fam_subjs, ]
res2_clean = res_clean[fam_subjs, ]

write.csv(res2, file='~/data/heritability_change/rsfmri_3min_n114_slopes.csv',
          row.names=F, na='', quote=F)
write.csv(res2_clean, file='~/data/heritability_change/rsfmri_3min_n114_slopesClean.csv',
          row.names=F, na='', quote=F)
write.table(var_names, file='~/data/heritability_change/aparc.txt',
            row.names=F, col.names=F, quote=F)
```

And I ran the code above for the Z version as well, and aparc.2009s.

# 2019-04-03 09:41:52

Let's then run some SOLAR to see what we get. I'm currently kicked out of the
cluster, so I'll need to get creative about how to run all these variables. I'll
do something similar to run_solar_voxel.sh:

```bash
phen_file=rsfmri_3min_n114_slopesClean
tmp_dir=~/data/heritability_change
solar_dir=~/data/heritability_change
mkdir ${tmp_dir}/${phen_file}
mkdir /tmp/${phen_file}
for vox in `cat ~/data/heritability_change/aparc.txt`; do
    mkdir /tmp/${phen_file}/${vox};
    cp ${solar_dir}/pedigree.csv ${solar_dir}/procs.tcl ${solar_dir}/${phen_file}.csv /tmp/${phen_file}/${vox}/;
    cd /tmp/${phen_file}/${vox}/;
    solar run_phen_var $phen_file $vox;
    mv /tmp/${phen_file}/${vox}/i_${vox}/polygenic.out ${tmp_dir}/${phen_file}/${vox}_polygenic.out;
done
```

And run the same as above for the Clean version as well, just for comparison.

While we wait on results, let's run the association analysis.

```r
library(nlme)
a = read.csv('~/data/heritability_change/resting_demo.csv')
a$famID = a$Extended.ID...FamilyIDs
a[is.na(a$Extended.ID...FamilyIDs), 'famID'] = a[is.na(a$Extended.ID...FamilyIDs), 'Nuclear.ID...FamilyIDs']
famids = unique(a[, c('famID', 'Medical.Record...MRN')])
data = read.csv('~/data/heritability_change/rsfmri_3min_assoc_n231_slopesClean.csv')
data$sex = as.factor(data$sex)
data = merge(data, famids, by.x='ID', by.y='Medical.Record...MRN', all.x=T)
b = read.csv('~/data/heritability_change/fmri_corr_tables/pearson_3min_n462_aparc.csv')
var_names = colnames(b)[2:ncol(b)]
out_fname = '~/data/heritability_change/assoc_LME_3min_n231_pearsonSlopesClean.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline', 'DX', 'DX2')
targets = var_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-try(lme(as.formula(fm_str), data, ~1|famID, na.action=na.omit))
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
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX=='ADHD', ]
out_fname = '~/data/heritability_change/assoc_LME_3min_n231_pearsonSlopesClean_dx1.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = var_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|famID, na.action=na.omit))
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
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX2=='ADHD', ]
out_fname = '~/data/heritability_change/assoc_LME_3min_n231_pearsonSlopesClean_dx2.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = var_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|famID, na.action=na.omit))
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
}
write.csv(hold, out_fname, row.names=F)
```

# 2019-04-04 10:19:30

The copy and paste was getting a bit too much, so I created
fmri/run_change_association.R to make it more dynamic.

```r
source('~/research_code/fmri/run_change_association.R')
```

Now, it looks like there are two pairs in the 114 set who are being run in
SOLAR as unrelated. Let's figure out who they are and remove them from the
heritability set, just to see how the results hold:

```r
a = read.csv('~/data/heritability_change/rsfmri_3min_n114_slopes.csv')
a = a[a$ID != 7221745, ]
write.csv(a, file='~/data/heritability_change/rsfmri_3min_n113_slopes.csv',
          row.names=F, na='', quote=F)
```

And do the same for the Z and Clean sets, to re-run the SOLAR analysis.

Then, the goal is to compile the SOLAR results for the n113 sets (4 in aparc),
and see if any of them have intersections with the association results.

```r
a = read.csv('~/data/heritability_change/assoc_LME_3min_n231_pearsonSlopes_dx2.csv')
idx = a$term!='sex2' & a$term!='(Intercept)'
b = a[idx,]
good_tests = which(b$p.value < .05)
assoc_conns = unique(as.character(b[good_tests,]$target))
sol = read.csv('~/data/heritability_change/polygen_results_rsfmri_3min_n113_slopes.csv')
her_conns = as.character(sol[which(sol$h_pval < .05), 'phen'])
good_set = intersect(her_conns, assoc_conns)
length(good_set)
```

So, we have 44 connections out of the possible 3741 that are both heritable and
associated with some sort of clinical variable. These results at the moment
don't mean much... it's just the draft code to further the analysis. I should
try to find association within the same clinical variables that came up for DTI,
and even try some sort of multiple comparisons for either association,
heritability, or both.

```r
lh = grepl(var_names, pattern="^ctx.lh\\S+TO.ctx.lh", perl=T)
rh = grepl(var_names, pattern="^ctx.rh\\S+TO.ctx.rh", perl=T)
ctx_ps = sol[lh | rh, 'h_pval']
sum(ctx_ps < .05)
```

# 2019-04-25 14:07:13

I think that if we use Meff to correct the fMRI results we might get something
good. First, looking at the actual labels we're running, I want to remove the
left and right VentralDC, and the Brainstem. That leaves only cortical
structures, and a few of the subcortical ones. Let's play with that then.

The Clean datasets have more heritable results than the nonClean ones, so let's
stick with that for now. 

```r
a = read.csv('~/data/heritability_change/rsfmri_3min_n113_slopesClean.csv')
conns = colnames(a)[3:3743]
idx1 = grepl(conns, pattern='Ventral')
idx2 = grepl(conns, pattern='Brain.Stem')
cc = cor(a[, conns[!(idx1 | idx2)]], use='na.or.complete')
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

That gives a Meff of 43.15, which is p < .001 using initial alpha of .05. Only
one connection survives... Not good. What if we go the other way around, looking
at the regressions? 

Let's see then how many connections are still significant then. 
.05/meff
head(cc)
dim(cc)
dim(a)
head(colnames(cc))
tail(colnames(cc))

```
# TODO

