# 2019-08-12 10:03:56

Let's try to re-run the Yeo masks and the MELODIC analysis on the non-Z results
of the dual regression, as recommended by the FSL folks
(https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/DualRegression/UserGuide). Maybe the
results will be better?

```bash
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_fancy.ica/
for m in `cat ../ids_1.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    for i in 3 22 18 1 39 27 10; do
        3dmaskdump -mask ../group_epi_mask_fancy.nii \
            -o dumps/${maskid}_IC${i}.txt dual/dr_stage2_${maskid}.nii.gz[${i}];
    done;
done
```

```bash
cd ~/data/heritability_change/xcp-36p_despike/yeo_masks/
for m in `cat ../ids_1.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    for i in {0..6}; do
        3dmaskdump -mask ../group_epi_mask_fancy.nii \
            -o dumps/${maskid}_net${i}.txt dual/dr_stage2_${maskid}.nii.gz[${i}];
    done;
done
```

Then, we collect our results in R:

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=231015
for (m in c(3, 22, 18, 1, 39, 27, 10)) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/groupmelodic_fancy.ica/dumps/%04d_IC%d.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/xcp-36p_despike/melodic_fancy_IC%d.rds', m)
     saveRDS(brain_data, file=fname)
}
```

Then, repeat the same for the yeo_masks dual regression.

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=231015
for (m in 0:6) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks/dumps/%04d_net%d.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks_fancy_net%d.rds', m)
     saveRDS(brain_data, file=fname)
}
```

Now it's just a matter of assigning MRNs, calculate slopes, and
prepare it for SOLAR voxelwise.

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv')
mydir='~/data/heritability_change/xcp-36p_despike/'
for (ic in c(3, 22, 18, 1, 39, 27, 10)) {
    fname = sprintf('%s/melodic_fancy_IC%d.rds', mydir, ic)
    b = readRDS(fname)
    var_names = colnames(b)[2:ncol(b)]
    df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id', all.x=F)

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
        y = mres[idx[2], var_names] - mres[idx[1], var_names]
        x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
        slopes = y / x
        row = c(row, slopes)
        for (t in c('SX_inatt', 'SX_HI', 'qc')) {
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
        print(nrow(res))
    }
    colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI', 'qc',
                                                'inatt_baseline',
                                                'HI_baseline', 'DX', 'DX2'))
    # we only open this in R, so it's OK to be RData to load faster
    fname = sprintf('%s/melodic_fancy_slopesNoZ_IC%d.rds', mydir, ic)
    saveRDS(res, file=fname)

    # and remove outliers
    res_clean = res
    for (t in var_names) {
        mydata = as.numeric(res_clean[, t])
        # identifying outliers
        ul = mean(mydata) + 3 * sd(mydata)
        ll = mean(mydata) - 3 * sd(mydata)
        bad_subjs = c(which(mydata < ll), which(mydata > ul))

        # remove within-variable outliers
        res_clean[bad_subjs, t] = NA
    }
    fname = sprintf('%s/melodic_fancy_slopesCleansNoZ_IC%d.rds', mydir, ic)
    saveRDS(res_clean, file=fname)

    # and make sure every family has at least two people
    good_nuclear = names(table(df2$Nuclear.ID...FamilyIDs))[table(df2$Nuclear.ID...FamilyIDs) >= 4]
    good_extended = names(table(df2$Extended.ID...FamilyIDs))[table(df2$Extended.ID...FamilyIDs) >= 4]
    keep_me = c()
    for (f in good_nuclear) {
        keep_me = c(keep_me, df2[which(df2$Nuclear.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    for (f in good_extended) {
        keep_me = c(keep_me, df2[which(df2$Extended.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    keep_me = unique(keep_me)

    fam_subjs = c()
    for (s in keep_me) {
        fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
    }
    res2 = res[fam_subjs, ]
    res2_clean = res_clean[fam_subjs, ]

    fname = sprintf('%s/melodic_fancy_slopesFamsNoZ_IC%d.csv', mydir, ic)
    write.csv(res2, file=fname, row.names=F, na='', quote=F)
    fname = sprintf('%s/melodic_fancy_slopesCleanFamsNoZ_IC%d.csv', mydir, ic)
    write.csv(res2_clean, file=fname, row.names=F, na='', quote=F)
}
```

And for the Yeo masks:

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv')
mydir='~/data/heritability_change/xcp-36p_despike/'
for (ic in 0:6) {
    fname = sprintf('%s/yeo_masks_fancy_net%d.rds', mydir, ic)
    b = readRDS(fname)
    var_names = colnames(b)[2:ncol(b)]
    df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id', all.x=F)

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
        y = mres[idx[2], var_names] - mres[idx[1], var_names]
        x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
        slopes = y / x
        row = c(row, slopes)
        for (t in c('SX_inatt', 'SX_HI', 'qc')) {
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
        print(nrow(res))
    }
    colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI', 'qc',
                                                'inatt_baseline',
                                                'HI_baseline', 'DX', 'DX2'))
    # we only open this in R, so it's OK to be RData to load faster
    fname = sprintf('%s/yeo_masks_fancy_slopesNoZ_net%d.rds', mydir, ic)
    saveRDS(res, file=fname)

    # and remove outliers
    res_clean = res
    for (t in var_names) {
        mydata = as.numeric(res_clean[, t])
        # identifying outliers
        ul = mean(mydata) + 3 * sd(mydata)
        ll = mean(mydata) - 3 * sd(mydata)
        bad_subjs = c(which(mydata < ll), which(mydata > ul))

        # remove within-variable outliers
        res_clean[bad_subjs, t] = NA
    }
    fname = sprintf('%s/yeo_masks_fancy_slopesCleansNoZ_net%d.rds', mydir, ic)
    saveRDS(res_clean, file=fname)

    # and make sure every family has at least two people
    good_nuclear = names(table(df2$Nuclear.ID...FamilyIDs))[table(df2$Nuclear.ID...FamilyIDs) >= 4]
    good_extended = names(table(df2$Extended.ID...FamilyIDs))[table(df2$Extended.ID...FamilyIDs) >= 4]
    keep_me = c()
    for (f in good_nuclear) {
        keep_me = c(keep_me, df2[which(df2$Nuclear.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    for (f in good_extended) {
        keep_me = c(keep_me, df2[which(df2$Extended.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    keep_me = unique(keep_me)

    fam_subjs = c()
    for (s in keep_me) {
        fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
    }
    res2 = res[fam_subjs, ]
    res2_clean = res_clean[fam_subjs, ]

    fname = sprintf('%s/yeo_masks_fancy_slopesFamsNoZ_net%d.csv', mydir, ic)
    write.csv(res2, file=fname, row.names=F, na='', quote=F)
    fname = sprintf('%s/yeo_masks_fancy_slopesCleanFamsNoZ_net%d.csv', mydir, ic)
    write.csv(res2_clean, file=fname, row.names=F, na='', quote=F)
}
```

I stopped this analysis because I've been playing a bit more with the gray mask.
But at least the csv files were created (desktop).