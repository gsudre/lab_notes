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