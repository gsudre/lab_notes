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
sol = read.csv('~/data/heritability_change/polygen_results_rsfmri_3min_n114_slopes.csv')
her_conns = as.character(sol[which(sol$h_pval < .05), 'phen'])
good_set = intersect(her_conns, assoc_conns)
length(good_set)
```

So, we have 44 connections out of the possible 3741 that are both heritable and
associated with some sort of clinical variable. These results at the moment
don't mean much... it's just the draft code to further the analysis. I'm still
waiting on SOLAR runs for the n113 sets. Then, I should try to find association
within the same clinical variables that came up for DTI, and even try some sort
of multiple comparisons for either association, heritability, or both.

# TODO

* Compile SOLAR results
* Clean up the file with 114 subjects to remove unrelated ones, but will need to
  use the unClean version for that... waiting for it to finish calculating
  slopes.
