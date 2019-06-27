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
fname = sprintf('~/data/heritability_change/rsfmri_AROMA%s_%dmin_best2scansSlopes_n%d.RData',
                        p, min_time, nrow(res))
save(res, file=fname)

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
fname = sprintf('~/data/heritability_change/rsfmri_AROMA%s_%dmin_best2scansSlopesClean_n%d.RData',
                        p, min_time, nrow(res_clean))
save(res_clean, file=fname)

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

fname = sprintf('~/data/heritability_change/rsfmri_AROMA%s_%dmin_best2scansFamsSlopes_n%d.csv',
                        p, min_time, nrow(res2))
write.csv(res2, file=fname, row.names=F, na='', quote=F)
fname = sprintf('~/data/heritability_change/rsfmri_AROMA%s_%dmin_best2scansFamsSlopesClean_n%d.csv',
                        p, min_time, nrow(res2_clean))
write.csv(res2_clean, file=fname, row.names=F, na='', quote=F)

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

# 2019-06-27 10:21:51

And we can run it somewhat smoothly if we batch it:

```bash
cd ~/data/heritability_change/
rm swarm.aroma
for f in `/bin/ls *best2scansFamsSlopes*csv`; do
    phen=`echo $f | sed "s/\.csv//"`;
    echo "bash ~/research_code/run_solar_parallel.sh $phen " \
        "~/data/heritability_change/power264_conns.txt" >> swarm.aroma;
done
swarm --gres=lscratch:10 -f swarm.aroma --module solar -g 10 -t 32 \
    --logdir=trash_solaroma --job-name solaroma --time=8:00:00 --merge-output
```

The script to filter down the scan was getting too cumbersome, so I created
~/research_code/fmri/filter_aroma_scans.R and also
~/research_code/fmri/create_aroma_slopes.R.



# TODO
 * check that we're not using same DOA for two different scans!










