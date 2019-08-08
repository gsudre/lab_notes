# 2019-08-08 12:53:51

Let's give it a try with ReHo and ALFF. We could play with smoothed and
non-smoothed, as it's not taking that long to run things in the cluster based on
the MELODIC results.

```bash
#desktop
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
cd ~/data/heritability_change/xcp-36p_despike
mkdir reho
for maskid in `cat ids_1.txt`; do
    m=`printf %04d $maskid`;
    echo $m;
    for s in '' '_sm6' 'Z' 'Z_sm6'; do
        3dmaskdump -mask group_epi_mask_fancy.nii \
            -o reho/${m}${s}.txt $mydir/sub-${m}/norm/sub-${m}_reho${s}Std.nii.gz;
    done;
done
```

```bash
#desktop
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
cd ~/data/heritability_change/xcp-36p_despike
mkdir alff
for maskid in `cat ids_1.txt`; do
    m=`printf %04d $maskid`;
    echo $m;
    for s in '' '_sm6' 'Z' 'Z_sm6'; do
        3dmaskdump -mask group_epi_mask_fancy.nii \
            -o alff/${m}${s}.txt $mydir/sub-${m}/norm/sub-${m}_alff${s}Std.nii.gz;
    done;
done
```

Then, we collect our results in R:

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=231015
conn='reho'
for (m in c('', 'Z', '_sm6', 'Z_sm6')) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/%s/%04d%s.txt', conn, maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/xcp-36p_despike/%s%s.RData.gz', conn, m)
     save(brain_data, file=fname, compress=T)
}
```

And calculate slopes:

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv')
mydir='~/data/heritability_change/xcp-36p_despike/'
# don't use a loop so we can get them all done in parallel in BW
conn = 'reho'
suf = ''

fname = sprintf('%s/%s%s.RData.gz', mydir, conn, suf)
load(fname)
b = brain_data
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
fname = sprintf('%s/%s_fancy_slopes%s.rds', mydir, conn, suf)
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
fname = sprintf('%s/%s_fancy_slopesClean%s.rds', mydir, conn, suf)
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

fname = sprintf('%s/%s_fancy_slopesFam%s.csv', mydir, conn, suf)
write.csv(res2, file=fname, row.names=F, na='', quote=F)
fname = sprintf('%s/%s_fancy_slopesCleanFam%s.csv', mydir, conn, suf)
write.csv(res2_clean, file=fname, row.names=F, na='', quote=F)
```

We can later to p25 FD threshold if we get somewhat decent results.

<!-- Finally, set up the swarms:

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in 2 27 10 4 31 29 7; do
    phen_file=melodic_fancyp25_slopesCleanFam_IC${i};
    jname=fancyp25_${i}c;
    swarm_file=swarm.${jname};

    rm -f $swarm_file;
    for vlist in `ls $PWD/vlist*txt`; do  # getting full path to files
        echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
    done;
    swarm --gres=lscratch:10 -f $swarm_file --module solar -t 16 -g 5 \
            --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
            --partition quick,norm
done
``` -->