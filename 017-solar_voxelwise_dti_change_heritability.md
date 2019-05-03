# 2019-05-02 15:29:04

I'm following the guidelines form a previous note in Evernote. The main script
to call is run_solar_voxel, and we need to make sure the solar directories ar
eproperly set up in the cluster.

But we also need to create a CSV that has voxelwise data. I don't want to fuss
with QC for now, maybe later. So, let's just take the 133 scans that gave the
tract result, and export those. Actually, we'll need the 566 scans fo get 283
slopes used in the association analysis. I can get a mask from that later, then.

Then, it's just a matter of dumping the data for all our IDs. But we first need
to save the IDs that will be used, especially because we'll need to calculate
the slopes for each voxel!

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# restrict based on QC
pct = m$missingVolumes / m$numVolumes
idx = m$norm.trans < 5 & m$norm.rot < .1 & pct < .15
m = m[idx,]
# down to 902 scans
keep_me = c()
for (s in unique(m$Medical.Record...MRN...Subjects)) {
    subj_scans = m[m$Medical.Record...MRN...Subjects==s, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                 format="%m/%d/%Y")
    if (length(dates) >= 2) {
        sdates = sort(dates)  # not sure why index.return is not working...
        # make sure there is at least 6 months between scans
        next_scan = 2
        while (((sdates[next_scan] - sdates[1]) < 180) && (next_scan < length(sdates))) {
            next_scan = next_scan + 1
        }
        first_scan_age = subj_scans[dates==sdates[1], 'age_at_scan...Scan...Subjects']
        if (((sdates[next_scan] - sdates[1]) >= 180) && (first_scan_age < 26)) {
            idx1 = which(dates == sdates[1])
            keep_me = c(keep_me, which(m$Mask.ID...Scan == subj_scans[idx1, 'Mask.ID...Scan']))
            idx2 = which(dates == sdates[next_scan])
            keep_me = c(keep_me, which(m$Mask.ID...Scan == subj_scans[idx2, 'Mask.ID...Scan']))
        }
    }
}
m2 = m[keep_me, ]
# down to 566 scans, as we're not using tract data at this point
write.table(sort(m2$Mask.ID...Scan), file='~/data/heritability_change/maskids_566.txt',
    col.names=F, row.names=F)
```

Those are the scans that will generate the slopes later. I then added a few
zeros manually.

As I copy them to the cluster, let's go ahead and get started on the mask:

```bash
sed "s/$/_tensor_diffeo.nii.gz/" maskids_566.txt > subjs_tensor.txt;
sed "s/^/..\/analysis_may2017\//" subjs_tensor.txt > subjs_tensor2.txt;
TVMean -in subjs_tensor2.txt -out mean_final_high_res.nii.gz
TVtool -in mean_final_high_res.nii.gz -fa
tbss_skeleton -i mean_final_high_res -o mean_fa_skeleton
3dcalc -a mean_fa_skeleton.nii.gz -prefix mean_566_fa_skeleton_mask.nii.gz \
    -expr 'step(a-.2)'
dti_dir=/mnt/shaw/dti_robust_tsa/analysis_may2017/
mkdir dti_voxels
for maskid in `cat maskids_566.txt`; do
     3dmaskdump -mask mean_566_fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_fa.txt ${dti_dir}/${maskid}_tensor_diffeo_fa.nii.gz;
     3dmaskdump -mask mean_566_fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_ad.txt ${dti_dir}/${maskid}_tensor_diffeo_ad.nii.gz;
     3dmaskdump -mask mean_566_fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_rd.txt ${dti_dir}/${maskid}_tensor_diffeo_rd.nii.gz;
done
```

# 2019-05-03 10:07:04

The mask that was created had way too many voxels. Over 30K! That's because it
was using lots of peripheral voxels. Can we just transport the FA
skeleton from TBSS to our space? Let's transform the FSL FA image to our group
image, and then use that transformation to transform the skeleton:

```bash
flirt -in /usr/local/neuro/fsl/data/standard/FMRIB58_FA_1mm.nii.gz \
    -ref mean_final_high_res_fa.nii.gz \
    -out FMRIB58_FA_IN_groupTemplate.nii.gz -omat FMRIB58_to_group.mat \
    -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 \
    -searchrz -90 90 -dof 12 -interp trilinear

flirt -in /usr/local/neuro/fsl/data/standard/FMRIB58_FA-skeleton_1mm.nii.gz \
    -ref mean_final_high_res_fa.nii.gz \
    -out FMRIB58_FA-skeleton_inGroup.nii.gz -applyxfm \
    -init FMRIB58_to_group.mat -interp nearestneighbour

3dcalc -a FMRIB58_FA-skeleton_inGroup.nii.gz -prefix fa_skeleton_mask.nii.gz \
    -expr 'step(a-.2)'

dti_dir=/mnt/shaw/dti_robust_tsa/analysis_may2017/
mkdir dti_voxels
for maskid in `cat maskids_566.txt`; do
     3dmaskdump -mask fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_fa.txt ${dti_dir}/${maskid}_tensor_diffeo_fa.nii.gz;
     3dmaskdump -mask fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_ad.txt ${dti_dir}/${maskid}_tensor_diffeo_ad.nii.gz;
     3dmaskdump -mask fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_rd.txt ${dti_dir}/${maskid}_tensor_diffeo_rd.nii.gz;
done
```

Then we construct the data files in R:

```r
maskids = read.table('/mnt/shaw/dti_robust_tsa/heritability/maskids_566.txt')[, 1]
nvox=14681
for (m in c('fa', 'ad', 'rd')) {
     print(m)
     dti_data = matrix(nrow=length(maskids), ncol=nvox)
     for (s in 1:nrow(dti_data)) {
          a = read.table(sprintf('/mnt/shaw/dti_robust_tsa/heritability/dti_voxels/%04d_%s.txt',
                                 maskids[s], m))
          dti_data[s, ] = a[,4]
     }
     dti_data = cbind(maskids, dti_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%05d', d)))
     colnames(dti_data) = cnames
     write.csv(dti_data, file=sprintf('/mnt/shaw/dti_robust_tsa/heritability/dti_%s_voxelwise_n566_03052019.csv', m), row.names=F)
}
```

Now we have to run the exact same code we ran for tracts, but for the voxels:

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# restrict based on QC
pct = m$missingVolumes / m$numVolumes
idx = m$norm.trans < 5 & m$norm.rot < .1 & pct < .15
m = m[idx,]
# down to 902 scans
keep_me = c()
for (s in unique(m$Medical.Record...MRN...Subjects)) {
    subj_scans = m[m$Medical.Record...MRN...Subjects==s, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                 format="%m/%d/%Y")
    if (length(dates) >= 2) {
        sdates = sort(dates)  # not sure why index.return is not working...
        # make sure there is at least 6 months between scans
        next_scan = 2
        while (((sdates[next_scan] - sdates[1]) < 180) && (next_scan < length(sdates))) {
            next_scan = next_scan + 1
        }
        first_scan_age = subj_scans[dates==sdates[1], 'age_at_scan...Scan...Subjects']
        if (((sdates[next_scan] - sdates[1]) >= 180) && (first_scan_age < 26)) {
            idx1 = which(dates == sdates[1])
            keep_me = c(keep_me, which(m$Mask.ID...Scan == subj_scans[idx1, 'Mask.ID...Scan']))
            idx2 = which(dates == sdates[next_scan])
            keep_me = c(keep_me, which(m$Mask.ID...Scan == subj_scans[idx2, 'Mask.ID...Scan']))
        }
    }
}
m2 = m[keep_me, ]

source('~/research_code/lab_mgmt/merge_on_closest_date.R')

clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN...Subjects),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN...Subjects')
b = read.csv('/Volumes/Shaw/dti_robust_tsa/heritability/dti_fa_voxelwise_n566_03052019.csv')
tract_names = colnames(b)[2:ncol(b)]
df2 = merge(df, b, by.x='Mask.ID...Scan', by.y='mask.id')

library(MASS)
mres = df2
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
for (t in tract_names) {
    print(t)
    fm_str = sprintf('%s ~', t)
    fm_str = paste(fm_str,
                   'norm.rot + I(norm.rot^2) + norm.trans + I(norm.trans^2) +',
                   'missingVolumes')
    res.lm <- lm(as.formula(fm_str), data = mres, na.action=na.exclude)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    mres[, t] = residuals(step)
}
res = c()
for (s in unique(mres$Medical.Record...MRN...Subjects)) {
    idx = which(mres$Medical.Record...MRN...Subjects == s)
    row = c(s, unique(mres[idx, 'Sex...Subjects']))
    for (t in tract_names) {
        if (sum(is.na(mres[idx, t])) > 0) {
            # if any of the tract values is NA, make the slope NA
            row = c(row, NA)
        } else {
            fm_str = sprintf('%s ~ age_at_scan...Scan...Subjects', t)
           fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
           row = c(row, coefficients(fit)[2])
        }
    }
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan...Scan...Subjects', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan...Scan...Subjects'])
    row = c(row, mres[idx[base_DOA], 'SX_inatt'])
    row = c(row, mres[idx[base_DOA], 'SX_HI'])
    # DX1 is DSMV definition, DX2 will make SX >=4 as ADHD
    if (mres[idx[base_DOA], 'age_at_scan...Scan...Subjects'] < 16) {
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
}
colnames(res) = c('ID', 'sex', tract_names, c('SX_inatt', 'SX_HI',
                                              'inatt_baseline',
                                              'HI_baseline', 'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/dti_fa_residNoSex_OLS_naSlopes283.csv',
          row.names=F, na='', quote=F)

# and remove outliers
res_clean = res
for (t in tract_names) {
    mydata = as.numeric(res_clean[, t])
    # identifying outliers
    ul = mean(mydata) + 3 * sd(mydata)
    ll = mean(mydata) - 3 * sd(mydata)
    bad_subjs = c(which(mydata < ll), which(mydata > ul))

    # remove within-tract outliers
    res_clean[bad_subjs, t] = NA
}
write.csv(res_clean, file='~/data/heritability_change/dti_fa_residNoSex_OLS_naSlopes283Clean.csv',
          row.names=F, na='', quote=F)
```

And we need to repeat the same procedure for AD and RD.

Next step is to keep only the 133 we'll use for heritability:

```r
# make sure every family has at least two people
good_nuclear = names(table(m2$Nuclear.ID...FamilyIDs))[table(m2$Nuclear.ID...FamilyIDs) >= 4]
good_extended = names(table(m2$Extended.ID...FamilyIDs))[table(m2$Extended.ID...FamilyIDs) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, m2[which(m2$Nuclear.ID...FamilyIDs == f),
                            'Medical.Record...MRN...Subjects'])
}
for (f in good_extended) {
    keep_me = c(keep_me, m2[which(m2$Extended.ID...FamilyIDs == f),
                            'Medical.Record...MRN...Subjects'])
}
keep_me = unique(keep_me)

fam_subjs = c()
for (s in keep_me) {
    fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
}
res2 = res[fam_subjs, ]
res2 = res2[res2[, 'ID'] != 7221745, ]
write.csv(res2, file='~/data/heritability_change/dti_fa_residNoSex_OLS_naSlopes133.csv',
          row.names=F, na='', quote=F)
res2 = res_clean[fam_subjs, ]
res2 = res2[res2[, 'ID'] != 7221745, ]
write.csv(res2, file='~/data/heritability_change/dti_fa_residNoSex_OLS_naSlopes133Clean.csv',
          row.names=F, na='', quote=F)
```

Then, it's time to start running heritability per voxel. I ahve the cluster tied
up in bedpost right now, so let's see how long it'll take to just run the main
analysis in an interactive node, or even in my own laptop:

```bash
phen_file=dti_fa_residNoSex_OLS_naSlopes133Clean
tmp_dir=~/data/heritability_change
solar_dir=~/data/heritability_change
mkdir ${tmp_dir}/${phen_file}
mkdir /tmp/${phen_file}
for v in 1:14681; do
    vox=`print %05d $v`;
    mkdir /tmp/${phen_file}/${vox};
    cp ${solar_dir}/pedigree.csv ${solar_dir}/procs.tcl ${solar_dir}/${phen_file}.csv /tmp/${phen_file}/${vox}/;
    cd /tmp/${phen_file}/${vox}/;
    ~/Download/solar84/solar run_phen_var $phen_file $vox;
    mv /tmp/${phen_file}/${vox}/i_${vox}/polygenic.out ${tmp_dir}/${phen_file}/${vox}_polygenic.out;
done
```

Then, we run the nonClean version for comparison, and both for AD and RD.
Finally, we need to put the results back into .nii format and clusterize them.
Hopefully we'll have some decent results, which we will be able to check for
associations (or, in worst case scenario, run voxelwise association), and do
permutations to assess cluster significance.

