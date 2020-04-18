# 2020-04-16 19:37:36

Based on my chat with Philip:

```
I think the substantive things will be (1) looking at the overlap between ENIGMA QC protocols and ours- how many subjects would be retained by both; (2) combining R+L (at least for DTI)- maybe for rsfMRI; (3) doing the flow chart of how we rejected so many - so how many were rejected due to qc, how many sibs we lost because of this - having one bad scan sometimes meant we lost that subject and also their subs---and how many due to being not the best per family etc
```

And there's a summary file in
data/heritability_change_rev/reviewer_SOBP_comments.docx. Also, Philip thinks
this is a good idea to sumarize the models:

![](images/2020-04-16-19-40-48.png)

Let's first tacle the ENIGMA QC comparison:

Numbers DTI:

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# restrict based on QC
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
m = m[m$"age_at_scan...Scan...Subjects" < 18, ]
m = m[m$"goodVolumes" <= 61, ]
m = m[m$"numVolumes" < 80, ]

m$FAMID = m$Extended.ID...FamilyIDs
idx = is.na(m$FAMID)
m[idx,]$FAMID = m[idx,]$Nuclear.ID...FamilyIDs
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/heritability_change/clinical_09182019.csv')
df = mergeOnClosestDate(m, clin,
                        unique(m$Medical.Record...MRN...Subjects),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN...Subjects')
mres = df
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
mres$DX = NA
for (r in 1:nrow(mres)) {
    if (mres[r, 'age_at_scan...Scan...Subjects'] < 16) {
        if ((mres[r, 'SX_HI'] >= 6) || (mres[r, 'SX_inatt'] >= 6)) {
            mres[r, 'DX'] = 'ADHD'
        } else {
            mres[r, 'DX'] = 'NV'
        }
    } else {
        if ((mres[r, 'SX_HI'] >= 5) || (mres[r, 'SX_inatt'] >= 5)) {
            mres[r, 'DX'] = 'ADHD'
        } else {
            mres[r, 'DX'] = 'NV'
        }
    }
}
m=mres

print(sprintf('Start with < 18 scans: %d scans, %d subjects, %d families, %d ADHD',
              nrow(m), length(unique(m$Medical.Record...MRN...Subjects)),
              length(unique(m$FAMID)), sum(m$DX=='ADHD')))
```
[1] "Start with < 18 scans: 954 scans, 299 subjects, 215 families, 424 ADHD"

```r
library(solitude)
iso <- isolationForest$new()
iso$fit(m[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]
library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(m[, qc_vars], k = round(.5 * nrow(m)))

qtile=.95
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

print(sprintf('After qc_vars OD: %d scans, %d subjects, %d families, %d ADHD',
              nrow(m[idx,]), length(unique(m[idx,]$Medical.Record...MRN...Subjects)),
              length(unique(m[idx,]$FAMID)), sum(m[idx,]$DX=='ADHD')))
```

[1] "After qc_vars OD: 879 scans, 297 subjects, 214 families, 380 ADHD"

```r
tracts = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
# somehow I have two entries for 1418?
x = duplicated(tracts$id)
data = merge(m[idx,], tracts[!x, ], by.x='Mask.ID...Scan', by.y='id')
tract_names = colnames(tracts)[grepl(colnames(tracts), pattern="^ad") | 
                                grepl(colnames(tracts), pattern="^rd")]

iso <- isolationForest$new()
iso$fit(data[, tract_names])
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(data[, tract_names], k = round(.5 * nrow(data)))

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

print(sprintf('After tract OD: %d scans, %d subjects, %d families, %d ADHD',
              nrow(data[idx,]),
              length(unique(data[idx,]$Medical.Record...MRN...Subjects)),
              length(unique(data[idx,]$FAMID)),
              sum(data[idx,]$DX=='ADHD')))
```

[1] "After tract OD: 815 scans, 295 subjects, 213 families, 354 ADHD"

```r
num_scans = 2  # number of scans to select
data$scores = scores_lof
a = data[idx, ]
# removing people with less than num_scans scans
idx = which(table(a$Medical.Record...MRN)>=num_scans)
long_subjs = names(table(a$Medical.Record...MRN))[idx]
keep_me = c()
for (m in 1:nrow(a)) {
    if (a[m, ]$Medical.Record...MRN %in% long_subjs) {
        keep_me = c(keep_me, m)
    }
}
a = a[keep_me,]

print(sprintf('After removing kids with only 1 scan: %d scans, %d subjects, %d families, %d ADHD',
              nrow(a),
              length(unique(a$Medical.Record...MRN...Subjects)),
              length(unique(a$FAMID)),
              sum(a$DX=='ADHD')))
```

[1] "After removing kids with only 1 scan: 774 scans, 254 subjects, 186
families, 337 ADHD"

```r
keep_me = c()
for (s in unique(a$Medical.Record...MRN)) {
    found = F
    subj_idx = which(a$Medical.Record...MRN==s)
    subj_scans = a[subj_idx, ]
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
filtered_data = a[keep_me, ]

a = filtered_data
print(sprintf('After keeping 2 best scans: %d scans, %d subjects, %d families, %d ADHD',
              nrow(a),
              length(unique(a$Medical.Record...MRN...Subjects)),
              length(unique(a$FAMID)),
              sum(a$DX=='ADHD')))
```

[1] "After keeping 2 best scans: 504 scans, 252 subjects, 185 families, 213
ADHD"

Let's give the mean and sd of each qc var and tract_var, and then make
histograms of the tract_vars ala ENIGMA. Or, even better, show that the scans
that stayed are better than the ones that were removed?

```r
for (v in qc_vars) {
    print(sprintf('%s: %.2f (+- %.2f)', v, mean(a[, v]), sd(a[, v])))
}
```
[1] "meanX.trans: 0.18 (+- 0.13)"
[1] "meanY.trans: 0.38 (+- 0.23)"
[1] "meanZ.trans: 0.82 (+- 0.77)"
[1] "meanX.rot: 0.01 (+- 0.01)"
[1] "meanY.rot: 0.01 (+- 0.00)"
[1] "meanZ.rot: 0.00 (+- 0.00)"
[1] "goodVolumes: 59.67 (+- 1.55)"

```r
mres$kept = F
mres[mres$Mask.ID...Scan %in% a$Mask.ID...Scan, 'kept'] = T
for (v in qc_vars) {
    print(v)
    print(t.test(mres[mres$kept==T,v], mres[mres$kept==F,v]))
}
```

```
[1] "meanX.trans"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = -3.901, df = 628.93, p-value = 0.0001061
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.08156248 -0.02694185
sample estimates:
mean of x mean of y
0.1752001 0.2294523

[1] "meanY.trans"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = -4.597, df = 532.29, p-value = 5.36e-06
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.23135628 -0.09282376
sample estimates:
mean of x mean of y
0.3809131 0.5430031

[1] "meanZ.trans"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = -4.4124, df = 601.88, p-value = 1.211e-05
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.5710976 -0.2192998
sample estimates:
mean of x mean of y
0.8158318 1.2110305

[1] "meanX.rot"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = -5.1882, df = 587.54, p-value = 2.927e-07
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.010936408 -0.004930128
sample estimates:
 mean of x  mean of y
0.01289131 0.02082458

[1] "meanY.rot"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = -3.4931, df = 717.9, p-value = 0.0005068
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.0022381101 -0.0006274993
sample estimates:
  mean of x   mean of y
0.005377990 0.006810794

[1] "meanZ.rot"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = -4.6084, df = 582.23, p-value = 4.989e-06
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.003341304 -0.001344337
sample estimates:
  mean of x   mean of y
0.004142585 0.006485405

[1] "goodVolumes"

        Welch Two Sample t-test

data:  mres[mres$kept == T, v] and mres[mres$kept == F, v]
t = 6.9631, df = 561.16, p-value = 9.338e-12
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 1.029921 1.839285
sample estimates:
mean of x mean of y
  59.6746   58.2400
```

```r
a = read.csv('~/data/heritability_change_rev/dti_JHUtracts_ADRDonly_OD0.95_twoTimePoints_noOtherDX.csv')
tract_names = colnames(a)[grepl(colnames(a), pattern="^ad")]
par(mfrow=c(5, 4))
for (v in tract_names) {
    t_str = sprintf('%s: %.2f (+- %.2f)', v, mean(a[, v]), sd(a[, v]))
    hist(a[, v], breaks=25, main=t_str)
}
tract_names = colnames(a)[grepl(colnames(a), pattern="^rd")]
par(mfrow=c(5, 4))
for (v in tract_names) {
    t_str = sprintf('%s: %.2f (+- %.2f)', v, mean(a[, v]), sd(a[, v]))
    hist(a[, v], breaks=25, main=t_str)
}
```

![](images/2020-04-16-21-10-59.png)

![](images/2020-04-16-21-11-49.png)

# 2020-04-17 14:09:50

Some fMRI reconstruction:

```r
qtile = .95

demo = read.csv('~/data/heritability_change/resting_demo_07032019.csv')
cat(sprintf('Starting from %d scans\n', nrow(demo)))

# keeping it to kids only to make sure everyone has been processed
demo = demo[demo$age_at_scan < 18, ]
cat(sprintf('Down to %d to keep < 18 only\n', nrow(demo)))

print(sprintf('Start with < 18 scans: %d scans, %d subjects',
              nrow(demo), length(unique(demo$Medical.Record...MRN))))
```

[1] "Start with < 18 scans: 1306 scans, 497 subjects"

```r
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

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)

idx = scores_lof < thresh_lof & scores_if < thresh_if

qc_data_clean = qc_data[idx, ]
qc_data_clean$mask.id = as.numeric(gsub(qc_data_clean$id0,
                                        pattern='sub-', replacement=''))

df = merge(qc_data_clean, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

print(sprintf('After QC OD: %d scans, %d subjects',
              nrow(df), length(unique(df$Medical.Record...MRN))))
```

[1] "After QC OD: 662 scans, 296 subjects"

```r
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

fc[abs(fc) > .21] = 1
fc[fc < 1] = 0

net_dataP = c()
header = c()
for (i in 1:nnets) {
    for (j in i:nnets) {
        cat(sprintf('Evaluating connections from %s to %s\n',
                    net_names[i], net_names[j]))
        idx = (conn_map[,2]==net_names[i] & conn_map[,3]==net_names[j]) |
            (conn_map[,3]==net_names[i] & conn_map[,2]==net_names[j])
        res = apply(fc[, var_names[idx]], 1, sum, na.rm=T)
        net_dataP = cbind(net_dataP, res)
        header = c(header, sprintf('conn_%sTO%s', net_names[i],
                                                net_names[j]))
    }
}
colnames(net_dataP) = header
rownames(net_dataP) = qc_data_clean$id0

var_names = c("conn_DorsAttnTODorsAttn", "conn_DorsAttnTOSalVentAttn",
              "conn_DorsAttnTOCont", "conn_DorsAttnTODefault", "conn_SalVentAttnTOSalVentAttn", "conn_SalVentAttnTOCont",
              "conn_SalVentAttnTODefault", "conn_ContTOCont",
              "conn_ContTODefault", "conn_DefaultTODefault")

iso <- isolationForest$new()
iso$fit(as.data.frame(net_dataP[, var_names]))
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(net_dataP[, var_names], k = round(.5 * nrow(net_dataP)))

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)

idx = scores_lof < thresh_lof & scores_if < thresh_if
data = cbind(qc_data_clean[, c('id0', qc_vars)], net_dataP)
data = data[idx, ]

data$mask.id = as.numeric(gsub(data$id0, pattern='sub-', replacement=''))

df = merge(data, demo, by.x='mask.id', by.y='Mask.ID', all.x=T, all.y=F)

print(sprintf('After brain var OD: %d scans, %d subjects',
              nrow(df), length(unique(df$Medical.Record...MRN))))
```

[1] "After brain var OD: 610 scans, 291 subjects"

```r
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
print(sprintf('After removing subjects with one scan only: %s scans, %d subjects',
              nrow(df), length(unique(df$Medical.Record...MRN))))
```

[1] "After removing subjects with one scan only: 546 scans, 227 subjects"

```r
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
print(sprintf('After keeping 2 best: %s scans, %d subjects',
              nrow(filtered_data),
              length(unique(filtered_data$Medical.Record...MRN))))
```

[1] "After keeping 2 best: 452 scans, 226 subjects"

And we can run the same t-tests as we did for DTI:

```r
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
qc_data$mask.id = as.numeric(gsub(qc_data$id0, pattern='sub-', replacement=''))
qc_data$kept = F
qc_data[qc_data$mask.id %in% filtered_data$mask.id, 'kept'] = T
for (v in qc_vars) {
    print(v)
    print(t.test(qc_data[qc_data$kept==T,v],
                 qc_data[qc_data$kept==F,v]))
}
```

```
[1] "normCoverage"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = 3.2079, df = 401.28, p-value = 0.001444
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 0.001423151 0.005928369
sample estimates:
mean of x mean of y
0.9788721 0.9751963

[1] "meanDV"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = -2.174, df = 551.74, p-value = 0.03013
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.041358961 -0.002096491
sample estimates:
mean of x mean of y
 1.081117  1.102845

[1] "pctSpikesDV"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = -3.3329, df = 500.57, p-value = 0.000923
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.029178832 -0.007535926
sample estimates:
 mean of x  mean of y
0.06299038 0.08134776

[1] "motionDVCorrInit"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = -3.3506, df = 532.89, p-value = 0.0008634
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.07298044 -0.01903399
sample estimates:
mean of x mean of y
0.5444068 0.5904140

[1] "motionDVCorrFinal"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = -0.90066, df = 538.32, p-value = 0.3682
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.03318850  0.01232214
sample estimates:
mean of x mean of y
0.3334610 0.3438941

[1] "pctSpikesRMS"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = -3.8762, df = 482.28, p-value = 0.0001208
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1286172 -0.0420862
sample estimates:
mean of x mean of y
0.2947551 0.3801068

[1] "relMeanRMSMotion"

        Welch Two Sample t-test

data:  qc_data[qc_data$kept == T, v] and qc_data[qc_data$kept == F, v]
t = -5.0948, df = 299.38, p-value = 6.191e-07
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4759653 -0.2107258
sample estimates:
mean of x mean of y
0.3155604 0.6589060
```

## Testing AD/RD variable normality

Let's run some test of normality in our AD/RD variables like ENIGMA:

```r
a = read.csv('~/data/heritability_change_rev/dti_JHUtracts_ADRDonly_OD0.95_twoTimePoints_noOtherDX.csv')
tract_names = colnames(a)[grepl(colnames(a), pattern="^ad")]
for (v in tract_names) {
    print(shapiro.test(a[, v]))
}
```

Actually, that didn't turn out favorably. Several variables are not normal,
which is fine. There are no huge outliers either, and the ENIGMA protocol has no
formal test of normality. Only visual inspection of the histograms, which we
did.


# TODO
 * make table comparing good scans to bad scans in terms of QC variables