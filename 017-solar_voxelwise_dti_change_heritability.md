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

Then we construct the data files in R:

<!-- ```r
pheno = read.csv('~/data/baseline_prediction/dti_gf_08152018_263timeDiff12mo.csv')
nvox=12022
for (m in c('fa', 'ad', 'rd')) {
     print(m)
     dti_data = matrix(nrow=nrow(pheno), ncol=nvox)
     for (s in 1:nrow(dti_data)) {
          a = read.table(sprintf('~/data/baseline_prediction/dti_voxels/%04d_%s.txt',pheno[s,]$mask.id, m))
          dti_data[s, ] = a[,4]
     }
     dti_data = cbind(pheno$mask.id, dti_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%05d', d)))
     colnames(dti_data) = cnames
     write.csv(dti_data, file=sprintf('~/data/baseline_prediction/dti_%s_voxelwise_n263_08152018.csv', m), row.names=F)
}
``` -->

<!-- ```r
clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN...Subjects),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN...Subjects')
b = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
tract_names = colnames(b)[2:ncol(b)]
df2 = merge(df, b, by.x='Mask.ID...Scan', by.y='id')

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

# and remove outliers
for (t in tract_names) {
    mydata = as.numeric(res[, t])
    # identifying outliers
    ul = mean(mydata) + 3 * sd(mydata)
    ll = mean(mydata) - 3 * sd(mydata)
    bad_subjs = c(which(mydata < ll), which(mydata > ul))

    # remove within-tract outliers
    res[bad_subjs, t] = NA
}
write.csv(res, file='~/data/heritability_change/dti_JHUtracts_residNoSex_OLS_naSlopes283.csv',
          row.names=F, na='', quote=F)
```
 -->


