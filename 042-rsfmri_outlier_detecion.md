# 2019-09-17 12:15:23

We'll do something similar to what we did for note 040, but this time for fMRI.

```r
demo = read.csv('~/data/heritability_change/resting_demo_07032019.csv')
cat(sprintf('Starting from %d scans\n', nrow(demo)))

# keeping it to kids only to make sure everyone has been processed
demo = demo[a$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(demo)))

# let's grab QC metrics on everyone
# note that this only works for non-censoring pipelines!
mydir = '/Volumes/Shaw/rsfmri_36P/xcpengine_output_fc-36p_despike/'
qc_data = c()
for (s in demo$Mask.ID) {
    subj = sprintf('sub-%04d', s)
    # if it processed all the way
    std_fname = sprintf('%s/%s/norm/%s_std.nii.gz', mydir, subj, subj)
    if (file.exists(std_fname)) {
        subj_data = read.csv(sprintf('%s/%s/%s_quality.csv', mydir, subj, subj))
        qc_data = rbind(qc_data, subj_data)
    }
}
```

So, we started with 1783 scans, and went down to 1306 to keep only kids. Then,
only 793 were processed all the way to the end in the fc-36p_despike pipeline. 

```r
# # column names that are not constant
# qc_vars = names(which(apply(qc_data[,2:ncol(data)], 2, var, na.rm=TRUE)!=0))

# have some higly correlated qc variables, so let's remove the worse offenders (anything above abs(.8))
qc_vars = c('normCoverage', 'meanDV', 'pctSpikesDV',
            'motionDVCorrInit',
            'motionDVCorrFinal', "pctSpikesRMS", "relMeanRMSMotion")

library(solitude)
iso <- isolationForest$new()
iso$fit(qc_data[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]

library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(qc_data[, qc_vars], k = round(.5 * nrow(qc_data)))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
```

![](images/2019-09-19-13-12-15.png)

I think .5 seems like a reasonable threshold for IF, and 2 for LOF.

Now, we need to decide which dataset to use. I'll expand slowly, starting with
7x7, then 17x17, then roix7, and finally roix17. For each I can make posOnly as
well as everything, and I should also compute Fam datasets, since SOLAR seems to
be screwing that up.

```r
idx = scores_lof < 2 & scores_if < .5
nrois = 100

fname = sprintf('~/research_code/fmri/Schaefer2018_%dParcels_7Networks_order.txt',
                nrois)
nets = read.table(fname)
all_net_names = sapply(as.character(unique(nets[,2])),
                       function(y) strsplit(x=y, split='_')[[1]][3])
net_names = unique(all_net_names)
nnets = length(net_names)

# figure out which connection goes to which network
cat('Creating connection map...\n')
nverts = nrow(nets)
cnt = 1
conn_map = c()
for (i in 1:(nverts-1)) {
    for (j in (i+1):nverts) {
        conn = sprintf('conn%d', cnt)
        conn_map = rbind(conn_map, c(conn, all_net_names[i], all_net_names[j]))
        cnt = cnt + 1
    }
}

qc_data_clean = qc_data[idx, ]
fc = c()
for (s in qc_data_clean$id0) {
    fname = sprintf('%s/%s/fcon/schaefer%d/%s_schaefer%d_network.txt',
                                mydir, s, nrois, s, nrois)
    subj_data = read.table(fname)[, 1]
    fc = cbind(fc, subj_data)
}
fc = t(fc)
var_names = sapply(1:ncol(fc), function(x) sprintf('conn%d', x))
colnames(fc) = var_names
fcP = fc
fcP[fc<0] = NA

net_data = c()
net_dataP = c()
header = c()
for (i in 1:nnets) {
    for (j in i:nnets) {
        cat(sprintf('Evaluating connections from %s to %s\n',
                    net_names[i], net_names[j]))
        idx = (conn_map[,2]==net_names[i] | conn_map[,2]==net_names[j] |
            conn_map[,3]==net_names[i] | conn_map[,3]==net_names[j])
        res = apply(fc[, var_names[idx]], 1, mean, na.rm=T)
        net_data = cbind(net_data, res)
        res = apply(fcP[, var_names[idx]], 1, mean, na.rm=T)
        net_dataP = cbind(net_dataP, res)
        header = c(header, sprintf('conn_%sTO%s', net_names[i],
                                                net_names[j]))
    }
}
colnames(net_data) = header
rownames(net_data) = qc_data_clean$id0
colnames(net_dataP) = header
rownames(net_dataP) = qc_data_clean$id0
```

So we don't get confused, let's do the whole thing for the whole data, then we
redo it for positive only:

```r
iso <- isolationForest$new()
iso$fit(as.data.frame(net_data))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_data, k = round(.5 * nrow(net_data)))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='LOF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='IF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
```

![](images/2019-09-19-13-41-36.png)

I'll go with 1.25 for LOF and .55 for IF:

```r
idx = scores_lof < 1.25 & scores_if < .55
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_data)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

num_scans = 2  # number of scans to select
df$scores = scores_lof[idx]

# removing people with less than num_scans scans
idx = which(table(df$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(df$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(df)) {
    if (df[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
df = df[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than %d scans\n',
            nrow(df), num_scans))
keep_me = c()
for (s in unique(df$Medical.Record...MRN)) {
    found = F
    subj_idx = which(df$Medical.Record...MRN==s)
    subj_scans = df[subj_idx, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                    format="%m/%d/%Y")
    best_scans = sort(subj_scans$scores, index.return=T)
    # make sure they are at least 6 months apart. This is the idea:
    # grab the best X scans. Check the time difference between them.
    # Any time the time difference is not enough, remove the worse
    # scan and replace by the next in line. Keep doing this until
    # the time difference is enough between all scans, or we run out
    # of scans
    cur_scan = 1
    last_scan = num_scans
    cur_choice = best_scans$ix[cur_scan:last_scan]
    while (!found && last_scan <= nrow(subj_scans)) {
        time_diffs = abs(diff(dates[cur_choice]))
        if (all(time_diffs > 180)) {
            found = TRUE
        } else {
            # figure out which scan to remove. If there is more than one
            # to be removed, it will be taken care in the next iteration
            bad_diff = which.min(time_diffs)
            if (subj_scans$scores[cur_choice[bad_diff]] >
                subj_scans$scores[cur_choice[bad_diff + 1]]) {
                rm_scan = cur_choice[bad_diff]
            } else {
                rm_scan = cur_choice[bad_diff + 1]
            }
            last_scan = last_scan + 1
            if (last_scan <= nrow(subj_scans)) {
                cur_choice[cur_choice == rm_scan] = best_scans$ix[last_scan]
            }
        }
    }
    if (found) {
        keep_me = c(keep_me, subj_idx[cur_choice])
    }
}
filtered_data = df[keep_me, ]
```

So, in the end we still have 448 scans, which is much more than what we had
before. Let's try using that for heritability and also for association (just to
simplify the analysis), and let's see what we get.

Time to compute slopes. Let's re-use the code and keep some baseline metrics in
there as well:

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(filtered_data, clin,
                        unique(filtered_data$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
tract_names = header

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], c(tract_names, qc_vars)] - mres[idx[1], c(tract_names, qc_vars)]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], tract_names])
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
tract_base = sapply(tract_names, function(x) sprintf('%s_baseline', x))
colnames(res) = c('ID', 'sex', tract_names, qc_vars, c('SX_inatt', 'SX_HI',
                                              tract_base,
                                              'inatt_baseline',
                                              'HI_baseline',
                                              'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_7by7from100_OD_n224.csv',
          row.names=F, na='', quote=F)

junk = read.csv('~/data/heritability_change/rsfmri_7by7from100_OD_n224.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
tmp2 = merge(junk, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)
related = names(table(tmp2$FAMID))[table(tmp2$FAMID) >= 2]
keep_me = tmp2$FAMID %in% related
res2 = junk[keep_me, ]
write.csv(res2, file='~/data/heritability_change/rsfmri_7by7from100_OD_n224_Fam_n112.csv',
          row.names=F, na='', quote=F)
```

Let me run these through SOLAR, before I check the positive only dataset:

```bash
# interactive
cd ~/data/heritability_change
for t in "conn_VisTOVis" "conn_VisTOSomMot" "conn_VisTODorsAttn" \
    "conn_VisTOSalVentAttn" "conn_VisTOLimbic" "conn_VisTOCont" \
    "conn_VisTODefault" "conn_SomMotTOSomMot" "conn_SomMotTODorsAttn" \
    "conn_SomMotTOSalVentAttn" "conn_SomMotTOLimbic" "conn_SomMotTOCont" \
    "conn_SomMotTODefault" "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp rsfmri_7by7from100_OD_n224_Fam_n112 ${t};
done;
mv rsfmri_7by7from100_OD_n224_Fam_n112 ~/data/tmp/
cd ~/data/tmp/rsfmri_7by7from100_OD_n224_Fam_n112/
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py rsfmri_7by7from100_OD_n224_Fam_n112
```

Not much there at all. Let's try the positive-only dataset as well:

```r
iso <- isolationForest$new()
iso$fit(as.data.frame(net_dataP))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_dataP, k = round(.5 * nrow(net_data)))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='LOF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='IF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
```

![](images/2019-09-19-14-04-30.png)

```r
idx = scores_lof < 1.2 & scores_if < .55
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_dataP)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

num_scans = 2  # number of scans to select
df$scores = scores_lof[idx]

# removing people with less than num_scans scans
idx = which(table(df$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(df$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(df)) {
    if (df[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
df = df[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than %d scans\n',
            nrow(df), num_scans))
keep_me = c()
for (s in unique(df$Medical.Record...MRN)) {
    found = F
    subj_idx = which(df$Medical.Record...MRN==s)
    subj_scans = df[subj_idx, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                    format="%m/%d/%Y")
    best_scans = sort(subj_scans$scores, index.return=T)
    # make sure they are at least 6 months apart. This is the idea:
    # grab the best X scans. Check the time difference between them.
    # Any time the time difference is not enough, remove the worse
    # scan and replace by the next in line. Keep doing this until
    # the time difference is enough between all scans, or we run out
    # of scans
    cur_scan = 1
    last_scan = num_scans
    cur_choice = best_scans$ix[cur_scan:last_scan]
    while (!found && last_scan <= nrow(subj_scans)) {
        time_diffs = abs(diff(dates[cur_choice]))
        if (all(time_diffs > 180)) {
            found = TRUE
        } else {
            # figure out which scan to remove. If there is more than one
            # to be removed, it will be taken care in the next iteration
            bad_diff = which.min(time_diffs)
            if (subj_scans$scores[cur_choice[bad_diff]] >
                subj_scans$scores[cur_choice[bad_diff + 1]]) {
                rm_scan = cur_choice[bad_diff]
            } else {
                rm_scan = cur_choice[bad_diff + 1]
            }
            last_scan = last_scan + 1
            if (last_scan <= nrow(subj_scans)) {
                cur_choice[cur_choice == rm_scan] = best_scans$ix[last_scan]
            }
        }
    }
    if (found) {
        keep_me = c(keep_me, subj_idx[cur_choice])
    }
}
filtered_data = df[keep_me, ]

source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(filtered_data, clin,
                        unique(filtered_data$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
tract_names = header

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], c(tract_names, qc_vars)] - mres[idx[1], c(tract_names, qc_vars)]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], tract_names])
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
tract_base = sapply(tract_names, function(x) sprintf('%s_baseline', x))
colnames(res) = c('ID', 'sex', tract_names, qc_vars, c('SX_inatt', 'SX_HI',
                                              tract_base,
                                              'inatt_baseline',
                                              'HI_baseline',
                                              'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n200.csv',
          row.names=F, na='', quote=F)

data = read.csv('~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n200.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(data, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)
related = names(table(data$FAMID))[table(data$FAMID) >= 2]
keep_me = data$FAMID %in% related
data2 = data[keep_me, ]
write.csv(data2, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n200_Fam_n100.csv',
          row.names=F, na='', quote=F)
```

```bash
# interactive
cd ~/data/heritability_change
for t in "conn_VisTOVis" "conn_VisTOSomMot" "conn_VisTODorsAttn" \
    "conn_VisTOSalVentAttn" "conn_VisTOLimbic" "conn_VisTOCont" \
    "conn_VisTODefault" "conn_SomMotTOSomMot" "conn_SomMotTODorsAttn" \
    "conn_SomMotTOSalVentAttn" "conn_SomMotTOLimbic" "conn_SomMotTOCont" \
    "conn_SomMotTODefault" "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp rsfmri_7by7from100_OD_posOnly_n200_Fam_n100 ${t};
done;
mv rsfmri_7by7from100_OD_posOnly_n200_Fam_n100 ~/data/tmp/
cd ~/data/tmp/rsfmri_7by7from100_OD_posOnly_n200_Fam_n100/
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py rsfmri_7by7from100_OD_posOnly_n200_Fam_n100
```

The posOnly results are a bit better, but still nothing that good. We get
nominal significance from VAN to DMN and Limbic, and also from DAN to limbic,
but nothing that would survive.

Let's try the 17x17 approach:

```r
nrois = 100

fname = sprintf('~/research_code/fmri/Schaefer2018_%dParcels_17Networks_order.txt',
                nrois)
nets = read.table(fname)
all_net_names = sapply(as.character(unique(nets[,2])),
                       function(y) strsplit(x=y, split='_')[[1]][3])
net_names = unique(all_net_names)
nnets = length(net_names)

# figure out which connection goes to which network
cat('Creating connection map...\n')
nverts = nrow(nets)
cnt = 1
conn_map = c()
for (i in 1:(nverts-1)) {
    for (j in (i+1):nverts) {
        conn = sprintf('conn%d', cnt)
        conn_map = rbind(conn_map, c(conn, all_net_names[i], all_net_names[j]))
        cnt = cnt + 1
    }
}

net_data2 = c()
net_data2P = c()
header = c()
for (i in 1:nnets) {
    for (j in i:nnets) {
        cat(sprintf('Evaluating connections from %s to %s\n',
                    net_names[i], net_names[j]))
        idx = (conn_map[,2]==net_names[i] | conn_map[,2]==net_names[j] |
            conn_map[,3]==net_names[i] | conn_map[,3]==net_names[j])
        res = apply(fc[, var_names[idx]], 1, mean, na.rm=T)
        net_data2 = cbind(net_data2, res)
        res = apply(fcP[, var_names[idx]], 1, mean, na.rm=T)
        net_data2P = cbind(net_data2P, res)
        header = c(header, sprintf('conn_%sTO%s', net_names[i],
                                                net_names[j]))
    }
}
colnames(net_data2) = header
rownames(net_data2) = qc_data_clean$id0
colnames(net_data2P) = header
rownames(net_data2P) = qc_data_clean$id0
```

The same way we approached this before, let's start with everything, then we
jump into positive-only.

```r
iso <- isolationForest$new()
iso$fit(as.data.frame(net_data2))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_data2, k = round(.5 * nrow(net_data2)))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='LOF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='IF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
```

![](images/2019-09-19-14-34-54.png)

I'll go with 1.25 for LOF and .55 for IF:

```r
idx = scores_lof < 1.25 & scores_if < .55
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_data2)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

num_scans = 2  # number of scans to select
df$scores = scores_lof[idx]

# removing people with less than num_scans scans
idx = which(table(df$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(df$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(df)) {
    if (df[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
df = df[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than %d scans\n',
            nrow(df), num_scans))
keep_me = c()
for (s in unique(df$Medical.Record...MRN)) {
    found = F
    subj_idx = which(df$Medical.Record...MRN==s)
    subj_scans = df[subj_idx, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                    format="%m/%d/%Y")
    best_scans = sort(subj_scans$scores, index.return=T)
    # make sure they are at least 6 months apart. This is the idea:
    # grab the best X scans. Check the time difference between them.
    # Any time the time difference is not enough, remove the worse
    # scan and replace by the next in line. Keep doing this until
    # the time difference is enough between all scans, or we run out
    # of scans
    cur_scan = 1
    last_scan = num_scans
    cur_choice = best_scans$ix[cur_scan:last_scan]
    while (!found && last_scan <= nrow(subj_scans)) {
        time_diffs = abs(diff(dates[cur_choice]))
        if (all(time_diffs > 180)) {
            found = TRUE
        } else {
            # figure out which scan to remove. If there is more than one
            # to be removed, it will be taken care in the next iteration
            bad_diff = which.min(time_diffs)
            if (subj_scans$scores[cur_choice[bad_diff]] >
                subj_scans$scores[cur_choice[bad_diff + 1]]) {
                rm_scan = cur_choice[bad_diff]
            } else {
                rm_scan = cur_choice[bad_diff + 1]
            }
            last_scan = last_scan + 1
            if (last_scan <= nrow(subj_scans)) {
                cur_choice[cur_choice == rm_scan] = best_scans$ix[last_scan]
            }
        }
    }
    if (found) {
        keep_me = c(keep_me, subj_idx[cur_choice])
    }
}
filtered_data = df[keep_me, ]

source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(filtered_data, clin,
                        unique(filtered_data$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
tract_names = header

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], c(tract_names, qc_vars)] - mres[idx[1], c(tract_names, qc_vars)]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], tract_names])
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
tract_base = sapply(tract_names, function(x) sprintf('%s_baseline', x))
colnames(res) = c('ID', 'sex', tract_names, qc_vars, c('SX_inatt', 'SX_HI',
                                              tract_base,
                                              'inatt_baseline',
                                              'HI_baseline',
                                              'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_17by17from100_OD_n230.csv',
          row.names=F, na='', quote=F)

junk = read.csv('~/data/heritability_change/rsfmri_17by17from100_OD_n230.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
tmp2 = merge(junk, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)
related = names(table(tmp2$FAMID))[table(tmp2$FAMID) >= 2]
keep_me = tmp2$FAMID %in% related
res2 = junk[keep_me, ]
write.csv(res2, file='~/data/heritability_change/rsfmri_17by17from100_OD_n230_Fam_n120.csv',
          row.names=F, na='', quote=F)
```

```bash
# interactive
cd ~/data/heritability_change
for t in `cat ~/data/heritability_change/tract_names_17x17.txt`; do
        solar run_phen_var_OD_xcp rsfmri_17by17from100_OD_n230_Fam_n120 ${t};
done;
mv rsfmri_17by17from100_OD_n230_Fam_n120 ~/data/tmp/
cd ~/data/tmp/rsfmri_17by17from100_OD_n230_Fam_n120/
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py rsfmri_17by17from100_OD_n230_Fam_n120
```

![](images/2019-09-19-14-55-58.png)

We started getting some weird, too high heritability values again...

And as usual we check the positive-only as well:

```r
iso <- isolationForest$new()
iso$fit(as.data.frame(net_data2P))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_data2P, k = round(.5 * nrow(net_data2P)))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='LOF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='IF Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
```

![](images/2019-09-19-14-42-14.png)

```r
idx = scores_lof < 1.25 & scores_if < .55
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_data2P)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

num_scans = 2  # number of scans to select
df$scores = scores_lof[idx]

# removing people with less than num_scans scans
idx = which(table(df$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(df$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(df)) {
    if (df[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
df = df[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than %d scans\n',
            nrow(df), num_scans))
keep_me = c()
for (s in unique(df$Medical.Record...MRN)) {
    found = F
    subj_idx = which(df$Medical.Record...MRN==s)
    subj_scans = df[subj_idx, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                    format="%m/%d/%Y")
    best_scans = sort(subj_scans$scores, index.return=T)
    # make sure they are at least 6 months apart. This is the idea:
    # grab the best X scans. Check the time difference between them.
    # Any time the time difference is not enough, remove the worse
    # scan and replace by the next in line. Keep doing this until
    # the time difference is enough between all scans, or we run out
    # of scans
    cur_scan = 1
    last_scan = num_scans
    cur_choice = best_scans$ix[cur_scan:last_scan]
    while (!found && last_scan <= nrow(subj_scans)) {
        time_diffs = abs(diff(dates[cur_choice]))
        if (all(time_diffs > 180)) {
            found = TRUE
        } else {
            # figure out which scan to remove. If there is more than one
            # to be removed, it will be taken care in the next iteration
            bad_diff = which.min(time_diffs)
            if (subj_scans$scores[cur_choice[bad_diff]] >
                subj_scans$scores[cur_choice[bad_diff + 1]]) {
                rm_scan = cur_choice[bad_diff]
            } else {
                rm_scan = cur_choice[bad_diff + 1]
            }
            last_scan = last_scan + 1
            if (last_scan <= nrow(subj_scans)) {
                cur_choice[cur_choice == rm_scan] = best_scans$ix[last_scan]
            }
        }
    }
    if (found) {
        keep_me = c(keep_me, subj_idx[cur_choice])
    }
}
filtered_data = df[keep_me, ]

source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(filtered_data, clin,
                        unique(filtered_data$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
tract_names = header

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], c(tract_names, qc_vars)] - mres[idx[1], c(tract_names, qc_vars)]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], tract_names])
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
tract_base = sapply(tract_names, function(x) sprintf('%s_baseline', x))
colnames(res) = c('ID', 'sex', tract_names, qc_vars, c('SX_inatt', 'SX_HI',
                                              tract_base,
                                              'inatt_baseline',
                                              'HI_baseline',
                                              'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_17by17from100_OD_posOnly_n231.csv',
          row.names=F, na='', quote=F)

junk = read.csv('~/data/heritability_change/rsfmri_17by17from100_OD_posOnly_n231.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
tmp2 = merge(junk, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)
related = names(table(tmp2$FAMID))[table(tmp2$FAMID) >= 2]
keep_me = tmp2$FAMID %in% related
res2 = junk[keep_me, ]
write.csv(res2, file='~/data/heritability_change/rsfmri_17by17from100_OD_posOnly_n231_Fam_n124.csv',
          row.names=F, na='', quote=F)
write.table(tract_names, file='~/data/heritability_change/tract_names_17x17.txt',
            row.names=F, col.names=F, quote=F)
```

```bash
# interactive
cd ~/data/heritability_change
for t in `cat ~/data/heritability_change/tract_names_17x17.txt`; do
        solar run_phen_var_OD_xcp rsfmri_17by17from100_OD_posOnly_n231_Fam_n124 ${t};
done;
mv rsfmri_17by17from100_OD_posOnly_n231_Fam_n124 ~/data/tmp/
cd ~/data/tmp/rsfmri_17by17from100_OD_posOnly_n231_Fam_n124/
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py rsfmri_17by17from100_OD_posOnly_n231_Fam_n124
```

I'm getting those weird values again... not sure what's causing it. Maybe I need
to clean the data within phenotype? Or run the entire cohort instead?

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_7by7from100_OD_n224;
for t in "conn_VisTOVis" "conn_VisTOSomMot" "conn_VisTODorsAttn" \
    "conn_VisTOSalVentAttn" "conn_VisTOLimbic" "conn_VisTOCont" \
    "conn_VisTODefault" "conn_SomMotTOSomMot" "conn_SomMotTODorsAttn" \
    "conn_SomMotTOSalVentAttn" "conn_SomMotTOLimbic" "conn_SomMotTOCont" \
    "conn_SomMotTODefault" "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_7by7from100_OD_posOnly_n200;
for t in "conn_VisTOVis" "conn_VisTOSomMot" "conn_VisTODorsAttn" \
    "conn_VisTOSalVentAttn" "conn_VisTOLimbic" "conn_VisTOCont" \
    "conn_VisTODefault" "conn_SomMotTOSomMot" "conn_SomMotTODorsAttn" \
    "conn_SomMotTOSalVentAttn" "conn_SomMotTOLimbic" "conn_SomMotTOCont" \
    "conn_SomMotTODefault" "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_17by17from100_OD_n230;
for t in `cat ~/data/heritability_change/tract_names_17x17.txt`; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_17by17from100_OD_posOnly_n231;
for t in `cat ~/data/heritability_change/tract_names_17x17.txt`; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

If looking at positive only, we do have some results. At least the numbers are
not going crazy high anymore. And the DTI results using the entire cohort are
also decent. First, for 7by7:

![](images/2019-09-19-16-20-20.png)

Then 17 by 17:

![](images/2019-09-19-16-20-44.png)

I'm just not sure what to do to keep these results after corrections... maybe
FDR? Or even Meff?

# 2019-09-23 13:30:26

Phikip suggested I could leave just VAN, DAN, DMN, and cognitive. Maybe limbic
too, but the first four make up the current model of ADHD, with DMN linking to
the task-positive networks. Even with that, FDR doesn't survive, and I doubt that
Meff is going to get below 5...

If I increase my sample size, I wonder if I can get better results. Let me see
if I automatize the elbow function I can get more people. Then, it's just a
matter of shoing that the chosen result still remain under more strict
thresholds.

```r
demo = read.csv('~/data/heritability_change/resting_demo_07032019.csv')
cat(sprintf('Starting from %d scans\n', nrow(demo)))

# keeping it to kids only to make sure everyone has been processed
demo = demo[a$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(demo)))

# let's grab QC metrics on everyone
# note that this only works for non-censoring pipelines!
mydir = '/Volumes/Shaw/rsfmri_36P/xcpengine_output_fc-36p_despike/'
qc_data = c()
for (s in demo$Mask.ID) {
    subj = sprintf('sub-%04d', s)
    # if it processed all the way
    std_fname = sprintf('%s/%s/norm/%s_std.nii.gz', mydir, subj, subj)
    if (file.exists(std_fname)) {
        subj_data = read.csv(sprintf('%s/%s/%s_quality.csv', mydir, subj, subj))
        qc_data = rbind(qc_data, subj_data)
    }
}

# have some higly correlated qc variables, so let's remove the worse offenders (anything above abs(.8))
qc_vars = c('normCoverage', 'meanDV', 'pctSpikesDV',
            'motionDVCorrInit',
            'motionDVCorrFinal', "pctSpikesRMS", "relMeanRMSMotion")

library(solitude)
iso <- isolationForest$new()
iso$fit(qc_data[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]

library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(qc_data[, qc_vars], k = round(.5 * nrow(qc_data)))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
```

In the past I chose scores_lof < 1.25 and scores_if < .55. LEt's see where the
elbow is:

```r
# to get distance function
source('~/research_code/get_elbow.R')
x = sort(scores_lof)
y = 1:length(x)
dists = vector()
for (i in 1:length(x)) {
    dists = c(dists, distancePointSegment(x[i], y[i], x[1], y[1], max(x), max(y)))
}
pos = which.max(dists)
print(sprintf('Elbow LOF at %.2f', x[pos]))
x = sort(scores_if)
y = 1:length(x)
dists = vector()
for (i in 1:length(x)) {
    dists = c(dists, distancePointSegment(x[i], y[i], x[1], y[1], max(x), max(y)))
}
pos = which.max(dists)
print(sprintf('Elbow LOF at %.2f', x[pos]))
```

OK, now that I have the elbows, let's define the datasets again. I'll just do
posOnly, as that's what I got best results with:

```r
idx = scores_lof < 1.99 & scores_if < .48
nrois = 100

fname = sprintf('~/research_code/fmri/Schaefer2018_%dParcels_7Networks_order.txt',
                nrois)
nets = read.table(fname)
all_net_names = sapply(as.character(unique(nets[,2])),
                       function(y) strsplit(x=y, split='_')[[1]][3])
net_names = unique(all_net_names)
nnets = length(net_names)

# figure out which connection goes to which network
cat('Creating connection map...\n')
nverts = nrow(nets)
cnt = 1
conn_map = c()
for (i in 1:(nverts-1)) {
    for (j in (i+1):nverts) {
        conn = sprintf('conn%d', cnt)
        conn_map = rbind(conn_map, c(conn, all_net_names[i], all_net_names[j]))
        cnt = cnt + 1
    }
}

qc_data_clean = qc_data[idx, ]
fc = c()
for (s in qc_data_clean$id0) {
    fname = sprintf('%s/%s/fcon/schaefer%d/%s_schaefer%d_network.txt',
                                mydir, s, nrois, s, nrois)
    subj_data = read.table(fname)[, 1]
    fc = cbind(fc, subj_data)
}
fc = t(fc)
var_names = sapply(1:ncol(fc), function(x) sprintf('conn%d', x))
colnames(fc) = var_names
fcP = fc
fcP[fc<0] = NA

net_dataP = c()
header = c()
for (i in 1:nnets) {
    for (j in i:nnets) {
        cat(sprintf('Evaluating connections from %s to %s\n',
                    net_names[i], net_names[j]))
        idx = (conn_map[,2]==net_names[i] | conn_map[,2]==net_names[j] |
            conn_map[,3]==net_names[i] | conn_map[,3]==net_names[j])
        res = apply(fcP[, var_names[idx]], 1, mean, na.rm=T)
        net_dataP = cbind(net_dataP, res)
        header = c(header, sprintf('conn_%sTO%s', net_names[i],
                                                net_names[j]))
    }
}
colnames(net_dataP) = header
rownames(net_dataP) = qc_data_clean$id0
```

And we need to re-measure elbows.

```r
iso <- isolationForest$new()
iso$fit(as.data.frame(net_dataP))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_dataP, k = round(.5 * nrow(net_dataP)))

# to get distance function
source('~/research_code/get_elbow.R')
x = sort(scores_lof)
y = 1:length(x)
dists = vector()
for (i in 1:length(x)) {
    dists = c(dists, distancePointSegment(x[i], y[i], x[1], y[1], max(x), max(y)))
}
pos = which.max(dists)
elbow_lof=x[pos]
print(sprintf('Elbow LOF at %.2f', elbow_lof))
x = sort(scores_if)
y = 1:length(x)
dists = vector()
for (i in 1:length(x)) {
    dists = c(dists, distancePointSegment(x[i], y[i], x[1], y[1], max(x), max(y)))
}
pos = which.max(dists)
elbow_if=x[pos]
print(sprintf('Elbow IF at %.2f', elbow_if))

par(mfrow=c(2,2))
sscores = sort(scores_lof)
hist(sscores, breaks=20, main='Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Local Outlier Factor')
abline(v=elbow_lof, col='red')
sscores = sort(scores_if)
hist(sscores, breaks=20, main='Anomaly score distribution')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main='Isolation Forest')
abline(v=elbow_if, col='red')
```

![](images/2019-09-23-14-07-15.png)

Hum, I don't think I'll end up with more people... but maybe better data?

```r
idx = scores_lof < 1.34 & scores_if < .46
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_dataP)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

num_scans = 2  # number of scans to select
df$scores = scores_lof[idx]

# removing people with less than num_scans scans
idx = which(table(df$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(df$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(df)) {
    if (df[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
df = df[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than %d scans\n',
            nrow(df), num_scans))
keep_me = c()
for (s in unique(df$Medical.Record...MRN)) {
    found = F
    subj_idx = which(df$Medical.Record...MRN==s)
    subj_scans = df[subj_idx, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                    format="%m/%d/%Y")
    best_scans = sort(subj_scans$scores, index.return=T)
    # make sure they are at least 6 months apart. This is the idea:
    # grab the best X scans. Check the time difference between them.
    # Any time the time difference is not enough, remove the worse
    # scan and replace by the next in line. Keep doing this until
    # the time difference is enough between all scans, or we run out
    # of scans
    cur_scan = 1
    last_scan = num_scans
    cur_choice = best_scans$ix[cur_scan:last_scan]
    while (!found && last_scan <= nrow(subj_scans)) {
        time_diffs = abs(diff(dates[cur_choice]))
        if (all(time_diffs > 180)) {
            found = TRUE
        } else {
            # figure out which scan to remove. If there is more than one
            # to be removed, it will be taken care in the next iteration
            bad_diff = which.min(time_diffs)
            if (subj_scans$scores[cur_choice[bad_diff]] >
                subj_scans$scores[cur_choice[bad_diff + 1]]) {
                rm_scan = cur_choice[bad_diff]
            } else {
                rm_scan = cur_choice[bad_diff + 1]
            }
            last_scan = last_scan + 1
            if (last_scan <= nrow(subj_scans)) {
                cur_choice[cur_choice == rm_scan] = best_scans$ix[last_scan]
            }
        }
    }
    if (found) {
        keep_me = c(keep_me, subj_idx[cur_choice])
    }
}
filtered_data = df[keep_me, ]

source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(filtered_data, clin,
                        unique(filtered_data$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
tract_names = header

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], c(tract_names, qc_vars)] - mres[idx[1], c(tract_names, qc_vars)]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], tract_names])
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
tract_base = sapply(tract_names, function(x) sprintf('%s_baseline', x))
colnames(res) = c('ID', 'sex', tract_names, qc_vars, c('SX_inatt', 'SX_HI',
                                              tract_base,
                                              'inatt_baseline',
                                              'HI_baseline',
                                              'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n178.csv',
          row.names=F, na='', quote=F)

data = read.csv('~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n178.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(data, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)
related = names(table(data$FAMID))[table(data$FAMID) >= 2]
keep_me = data$FAMID %in% related
data2 = data[keep_me, ]
write.csv(data2, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n178_Fam_n75.csv',
          row.names=F, na='', quote=F)
```

And finally we run SOLAR, but only for the interesting connections:

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_7by7from100_OD_posOnly_n178;
for t in "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

Wow... that wasn't good at all. Removed all my results... quite weird.

Let's go the other way around. Instead of using the elbow, let's determine the
percentile of the distribution. 

# 2019-09-24 09:15:37

```r
demo = read.csv('~/data/heritability_change/resting_demo_07032019.csv')
cat(sprintf('Starting from %d scans\n', nrow(demo)))

# keeping it to kids only to make sure everyone has been processed
demo = demo[a$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(demo)))

# let's grab QC metrics on everyone
# note that this only works for non-censoring pipelines!
mydir = '/Volumes/Shaw/rsfmri_36P/xcpengine_output_fc-36p_despike/'
qc_data = c()
for (s in demo$Mask.ID) {
    subj = sprintf('sub-%04d', s)
    # if it processed all the way
    std_fname = sprintf('%s/%s/norm/%s_std.nii.gz', mydir, subj, subj)
    if (file.exists(std_fname)) {
        subj_data = read.csv(sprintf('%s/%s/%s_quality.csv', mydir, subj, subj))
        qc_data = rbind(qc_data, subj_data)
    }
}

# have some higly correlated qc variables, so let's remove the worse offenders (anything above abs(.8))
qc_vars = c('normCoverage', 'meanDV', 'pctSpikesDV',
            'motionDVCorrInit',
            'motionDVCorrFinal', "pctSpikesRMS", "relMeanRMSMotion")

library(solitude)
iso <- isolationForest$new()
iso$fit(qc_data[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]

library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(qc_data[, qc_vars], k = round(.5 * nrow(qc_data)))

qtile = .95
par(mfrow=c(2,2))
sscores = sort(scores_lof)
thresh_lof = quantile(scores_lof, qtile)
hist(sscores, breaks=20, main='Local Outlier Factor: anomaly score distribution')
abline(v=thresh_lof, col='red')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main=sprintf('LOF %.2f pctile threshold = %.2f', qtile, thresh_lof))
abline(v=thresh_lof, col='red')
sscores = sort(scores_if)
thresh_if = quantile(scores_if, qtile)
hist(sscores, breaks=20, main='Isolation Forest: anomaly score distribution')
abline(v=thresh_if, col='red')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main=sprintf('IF %.2f pctile threshold = %.2f', qtile, thresh_if))
abline(v=thresh_if, col='red')
```

![](images/2019-09-24-09-51-17.png)

I think we'll end up with more people now, just not sure if better data:

```r
idx = scores_lof < 2.19 & scores_if < .54
nrois = 100

fname = sprintf('~/research_code/fmri/Schaefer2018_%dParcels_7Networks_order.txt',
                nrois)
nets = read.table(fname)
all_net_names = sapply(as.character(unique(nets[,2])),
                       function(y) strsplit(x=y, split='_')[[1]][3])
net_names = unique(all_net_names)
nnets = length(net_names)

# figure out which connection goes to which network
cat('Creating connection map...\n')
nverts = nrow(nets)
cnt = 1
conn_map = c()
for (i in 1:(nverts-1)) {
    for (j in (i+1):nverts) {
        conn = sprintf('conn%d', cnt)
        conn_map = rbind(conn_map, c(conn, all_net_names[i], all_net_names[j]))
        cnt = cnt + 1
    }
}

qc_data_clean = qc_data[idx, ]
fc = c()
for (s in qc_data_clean$id0) {
    fname = sprintf('%s/%s/fcon/schaefer%d/%s_schaefer%d_network.txt',
                                mydir, s, nrois, s, nrois)
    subj_data = read.table(fname)[, 1]
    fc = cbind(fc, subj_data)
}
fc = t(fc)
var_names = sapply(1:ncol(fc), function(x) sprintf('conn%d', x))
colnames(fc) = var_names
fcP = fc
fcP[fc<0] = NA

net_dataP = c()
header = c()
for (i in 1:nnets) {
    for (j in i:nnets) {
        cat(sprintf('Evaluating connections from %s to %s\n',
                    net_names[i], net_names[j]))
        idx = (conn_map[,2]==net_names[i] | conn_map[,2]==net_names[j] |
            conn_map[,3]==net_names[i] | conn_map[,3]==net_names[j])
        res = apply(fcP[, var_names[idx]], 1, mean, na.rm=T)
        net_dataP = cbind(net_dataP, res)
        header = c(header, sprintf('conn_%sTO%s', net_names[i],
                                                net_names[j]))
    }
}
colnames(net_dataP) = header
rownames(net_dataP) = qc_data_clean$id0

iso <- isolationForest$new()
iso$fit(as.data.frame(net_dataP))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_dataP, k = round(.5 * nrow(net_dataP)))

qtile = .95
par(mfrow=c(2,2))
sscores = sort(scores_lof)
thresh_lof = quantile(scores_lof, qtile)
hist(sscores, breaks=20, main='Local Outlier Factor: anomaly score distribution')
abline(v=thresh_lof, col='red')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main=sprintf('LOF %.2f pctile threshold = %.2f', qtile, thresh_lof))
abline(v=thresh_lof, col='red')
sscores = sort(scores_if)
thresh_if = quantile(scores_if, qtile)
hist(sscores, breaks=20, main='Isolation Forest: anomaly score distribution')
abline(v=thresh_if, col='red')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main=sprintf('IF %.2f pctile threshold = %.2f', qtile, thresh_if))
abline(v=thresh_if, col='red')
```

![](images/2019-09-24-09-55-57.png)

```r
idx = scores_lof < 1.57 & scores_if < .54
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_dataP)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

num_scans = 2  # number of scans to select
df$scores = scores_lof[idx]

# removing people with less than num_scans scans
idx = which(table(df$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(df$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(df)) {
    if (df[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
df = df[keep_me,]
cat(sprintf('Down to %d to keep only subjects with more than %d scans\n',
            nrow(df), num_scans))
keep_me = c()
for (s in unique(df$Medical.Record...MRN)) {
    found = F
    subj_idx = which(df$Medical.Record...MRN==s)
    subj_scans = df[subj_idx, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                    format="%m/%d/%Y")
    best_scans = sort(subj_scans$scores, index.return=T)
    # make sure they are at least 6 months apart. This is the idea:
    # grab the best X scans. Check the time difference between them.
    # Any time the time difference is not enough, remove the worse
    # scan and replace by the next in line. Keep doing this until
    # the time difference is enough between all scans, or we run out
    # of scans
    cur_scan = 1
    last_scan = num_scans
    cur_choice = best_scans$ix[cur_scan:last_scan]
    while (!found && last_scan <= nrow(subj_scans)) {
        time_diffs = abs(diff(dates[cur_choice]))
        if (all(time_diffs > 180)) {
            found = TRUE
        } else {
            # figure out which scan to remove. If there is more than one
            # to be removed, it will be taken care in the next iteration
            bad_diff = which.min(time_diffs)
            if (subj_scans$scores[cur_choice[bad_diff]] >
                subj_scans$scores[cur_choice[bad_diff + 1]]) {
                rm_scan = cur_choice[bad_diff]
            } else {
                rm_scan = cur_choice[bad_diff + 1]
            }
            last_scan = last_scan + 1
            if (last_scan <= nrow(subj_scans)) {
                cur_choice[cur_choice == rm_scan] = best_scans$ix[last_scan]
            }
        }
    }
    if (found) {
        keep_me = c(keep_me, subj_idx[cur_choice])
    }
}
filtered_data = df[keep_me, ]

source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(filtered_data, clin,
                        unique(filtered_data$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
tract_names = header

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], c(tract_names, qc_vars)] - mres[idx[1], c(tract_names, qc_vars)]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], tract_names])
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
tract_base = sapply(tract_names, function(x) sprintf('%s_baseline', x))
colnames(res) = c('ID', 'sex', tract_names, qc_vars, c('SX_inatt', 'SX_HI',
                                              tract_base,
                                              'inatt_baseline',
                                              'HI_baseline',
                                              'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n241.csv',
          row.names=F, na='', quote=F)

data = read.csv('~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n241.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(data, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)
related = names(table(data$FAMID))[table(data$FAMID) >= 2]
keep_me = data$FAMID %in% related
data2 = data[keep_me, ]
write.csv(data2, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n241_Fam_n125.csv',
          row.names=F, na='', quote=F)
```

And finally we run SOLAR, but only for the interesting connections:

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_7by7from100_OD_posOnly_n241;
for t in "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

Hum, that didn't do it either... concerning that there is a sweet spot of data
that this works. I can also go for .9 qtile jsut to see what we get. But maybe
there's something in the Fam only dataset?

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_7by7from100_OD_posOnly_n241_Fam_n125;
for t in "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

Still nothing... maybe median? Or do another pass of QC on slopes?

```r
data = read.csv('~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n241.csv')
tract_names = colnames(data)[16:30]

iso <- isolationForest$new()
iso$fit(as.data.frame(data[, tract_names]))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(data[, tract_names], k = round(.5 * nrow(data)))

qtile = .95
par(mfrow=c(2,2))
sscores = sort(scores_lof)
thresh_lof = quantile(scores_lof, qtile)
hist(sscores, breaks=20, main='Local Outlier Factor: anomaly score distribution')
abline(v=thresh_lof, col='red')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main=sprintf('LOF %.2f pctile threshold = %.2f', qtile, thresh_lof))
abline(v=thresh_lof, col='red')
sscores = sort(scores_if)
thresh_if = quantile(scores_if, qtile)
hist(sscores, breaks=20, main='Isolation Forest: anomaly score distribution')
abline(v=thresh_if, col='red')
plot(sscores, 1:length(sscores), ylab='Remaining observations',
     xlab='Score threshold', main=sprintf('IF %.2f pctile threshold = %.2f', qtile, thresh_if))
abline(v=thresh_if, col='red')
```

![](images/2019-09-24-10-14-59.png)

I did .95 there just for consistency, but I think we need to eventually be a bit
more stringent. Let's see what happens with these data though.

```r
idx = scores_lof < 2.3 & scores_if < .59
data2 = data[idx, ]
write.csv(data2, file='~/data/heritability_change/rsfmri_7by7from100_OD_posOnly_n241_QC_n226.csv',
          row.names=F, na='', quote=F)
```

```bash
# interactive
cd ~/data/heritability_change
phen=rsfmri_7by7from100_OD_posOnly_n241_QC_n226;
for t in "conn_DorsAttnTODorsAttn" \
    "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
    "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
    "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
    "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
    "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
    "conn_DefaultTODefault"; do
        solar run_phen_var_OD_xcp ${phen} ${t};
done;
mv ${phen} ~/data/tmp/
cd ~/data/tmp/${phen}
for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
python ~/research_code/compile_solar_multivar_results.py ${phen}
```

I put the R lines above in ~/research_code/fmri/make_outlier_detection_slopes.R.
I then created a few iterations of phenotype files... let's see if results come
up:

```bash
# interactive
for OD in 80 85 90 95; do
    for suf in '' '_Fam' '_QC'; do
        cd ~/data/heritability_change
        phen=rsfmri_7by7from100_OD0.${OD}${suf};
        for t in "conn_DorsAttnTODorsAttn" \
            "conn_DorsAttnTOSalVentAttn" "conn_DorsAttnTOLimbic" "conn_DorsAttnTOCont" \
            "conn_DorsAttnTODefault" "conn_SalVentAttnTOSalVentAttn" \
            "conn_SalVentAttnTOLimbic" "conn_SalVentAttnTOCont" \
            "conn_SalVentAttnTODefault" "conn_LimbicTOLimbic" "conn_LimbicTOCont" \
            "conn_LimbicTODefault" "conn_ContTOCont" "conn_ContTODefault" \
            "conn_DefaultTODefault"; do
                solar run_phen_var_OD_xcp ${phen} ${t};
        done;
        mv ${phen} ~/data/tmp/
        cd ~/data/tmp/${phen}
        for p in `/bin/ls`; do cp $p/polygenic.out ${p}_polygenic.out; done
        python ~/research_code/compile_solar_multivar_results.py ${phen}
    done;
done
```

And I also did it for median, just for kicks. And removed the positive
constraints as well. Does any of that work?

Let's merge all files to see if we get anything:

```bash
cd ~/data/tmp
echo "file,phen,n,h2r,h_pval,h2r_se,c2,c2_pval,high_kurtosis" > output.csv;
for f in `ls polygen_results_rsfmri_7by7from100_OD*.csv`; do
    # skip header
    for line in `tail -n +2 $f`; do
        echo $f,$line >> output.csv;
    done
done
```

Hum... there might be some results. But limbic doesn't seem to be helping. I'll
need to redo the data if I go with any of the QC sets!

After looking with lots of filters in Excel, I narrowed it down to .85 OD (by
far!), and then it can be one of these iterations:

![](images/2019-09-24-14-31-26.png)

![](images/2019-09-24-14-31-41.png)

![](images/2019-09-24-14-31-58.png)

Now I need to choose wha to go with. Maybe based on regressions? I might need to
remove Limbic as well to make some of them work out. Let's see what's best.

It looks like polygen_results_rsfmri_7by7from100_OD0.85_median.csv. with 143
subjects, has a strong showing. For example, doing FDR including limbic we get
the first 9 connections with q<.05:

```
conn_DorsAttnTOCont
conn_DorsAttnTODefault
conn_DorsAttnTOLimbic
conn_DorsAttnTOSalVentAttn
conn_SalVentAttnTOLimbic
conn_SalVentAttnTOCont
conn_DorsAttnTODorsAttn
conn_SalVentAttnTOSalVentAttn
conn_LimbicTODefault
```

Without limbic it's the first 7:

```
conn_DorsAttnTOCont
conn_DorsAttnTODefault
conn_DorsAttnTOSalVentAttn
conn_SalVentAttnTOCont
conn_DorsAttnTODorsAttn
conn_SalVentAttnTOSalVentAttn
conn_SalVentAttnTODefault
```

Just to keep score, we start with 1306 kiddie scans, but actually only 764 were
processed in this pipeline (due to several reasons). Some were just because the
scan was terrible, so we could count those into our exclusion numbers as well.
But even if not, after the qc_data outlier pass we have 564 left. After QCing on
the data, we have 435 scans left. Finally, that's 330 if we keep only subjects
with 2 or more scans, and 286 scans to calculate the slopes. So, based on data
quality, we removed 43% of the data. Even more, if we count the ones that we
could not even process. And, in the end, we only used 286 out of the 435
considered good, but that wasn't due to data quality anymore.

**NOTE that we do the data QC based on the 7x7 networks! That wasn;t the
original idea, but I forgot to change the script. We can still justify it by
saying that we want to QC on the entire brain, but we are only interested in
analyzing the 5 networks related to ADHD. Worst case scenario I re-run the whole
thing including only the networks we care about.**

OK, how about association to ADHD?

```r
library(nlme)
data = read.csv('~/data/heritability_change/rsfmri_7by7from100_OD0.85_median.csv')
tmp = read.csv('~/data/heritability_change/pedigree.csv')
data = merge(data, tmp[, c('ID', 'FAMID')], by='ID', all.x=T, all.y=F)

tract_names = colnames(data)[16:30]

# let's run just the residuals after stepwise regression
mres = data
library(MASS)
for (t in tract_names) {
    print(t)
    fm_str = sprintf('%s ~', t)
    fm_str = paste(fm_str, 'sex + normCoverage + meanDV + pctSpikesDV + motionDVCorrInit + motionDVCorrFinal + pctSpikesRMS + relMeanRMSMotion')
    res.lm <- lm(as.formula(fm_str), data = data)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    mres[, t] = residuals(step)
}

out_fname = '~/data/heritability_change/assoc_LME_rsfmri_7by7from100_OD0.85_median.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline', 'DX', 'DX2')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s', i, j)
        model1<-try(lme(as.formula(fm_str), mres, ~1|FAMID, na.action=na.omit))
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

data2 = mres[mres$DX=='ADHD', ]
out_fname = gsub(x=out_fname, pattern='.csv', '_dx1.csv')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|FAMID, na.action=na.omit))
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

data2 = mres[mres$DX2=='ADHD', ]
out_fname = gsub(x=out_fname, pattern='dx1', 'dx2')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|FAMID, na.action=na.omit))
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

We have something in DX2:

![](images/2019-09-24-15-59-05.png)

And VAN to Cog is heritable... nothing in DX1. The association with HI slopes
gets a bit stronger including NVs:

![](images/2019-09-24-16-01-12.png)

Now we should plot them to make sure it's not driven by outliers.

# TODO
 * maybe do it again only using for OD the networks we're actually using in the analysis?
 * run associations using SOLAR-selected covariates
