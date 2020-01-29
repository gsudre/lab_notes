# 2020-01-29 15:24:12

While we wait for mode to compute (069), let's start from where we are in terms
of models to use, and do the same for the anatomical data. I'll use our common
QC parameters for now, but still use the OD pipeline for scan selection. 

The idea then, so far, is to use just one per family so we don't have to worry
about mixed models, and figure out the best predictors for the logistic
regression between improvers and non-improvers. Then, we check how NVs compare
to them in those good variables.

We could potentially use the data we're trhwoing away for a pseudo-validation
later. But we'll see. 

```r
setwd('~/data/baseline_prediction/prs_start/')
clin_long = read.csv('long_clin_01062020_lt16.csv')
clin_long$SX_total = clin_long$SX_inatt + clin_long$SX_hi

winsorize = function(x, cut = 0.01){
  cut_point_top <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x, cut, na.rm = T)
  i = which(x >= cut_point_top) 
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom) 
  x[j] = cut_point_bottom
  return(x)
}

df = data.frame(MRN=unique(clin_long$MRN))
for (r in 1:nrow(df)) {
    subj_data = clin_long[clin_long$MRN==df$MRN[r], ]
    for (sx in c('inatt', 'hi', 'total')) {
        fit = lm(as.formula(sprintf('SX_%s ~ age', sx)), data=subj_data)
        df[r, sprintf('slope_%s', sx)] = fit$coefficients['age']
        base_row = which.min(subj_data$age)
        df[r, sprintf('base_%s', sx)] = subj_data[base_row, sprintf('SX_%s', sx)]
        last_row = which.max(subj_data$age)
        df[r, sprintf('last_%s', sx)] = subj_data[last_row, sprintf('SX_%s', sx)]
        df[r, 'base_age'] = subj_data[base_row, 'age']
        df[r, 'last_age'] = subj_data[last_row, 'age']
        df[r, 'sex'] = subj_data[last_row, 'sex']
    }
}
for (min_sx in c(0, 3, 4, 6)) {
    idx = df$base_inatt>=min_sx | df$base_hi>=min_sx
    for (sx in c('inatt', 'hi', 'total')) {
        df[, sprintf('slope_%s_GE%d_wp05', sx, min_sx)] = NA
        junk = winsorize(df[idx, sprintf('slope_%s', sx)], cut=.05)
        df[idx, sprintf('slope_%s_GE%d_wp05', sx, min_sx)] = junk
    }
}

demo = read.csv('prs_demo.csv')
# just to get FAMID, sex already there
df = merge(df, subset(demo, select=-sex), by='MRN')

# selecting best kid in family
df$bestInFamily = F
nvisits = table(clin_long$MRN)
df = merge(df, as.matrix(nvisits),
                 by.x='MRN', by.y=0)
colnames(df)[ncol(df)] = 'nvisits'
for (f in unique(df$FAMID)) {
    fam_rows = which(df$FAMID == f)
    fam_data = df[fam_rows,]
    if (nrow(fam_data) == 1) {
        df[fam_rows,]$bestInFamily = T
    } else {
        stotal = sort(fam_data$slope_total, index.return=T, decreasing=T)
        # if there's a tie
        if (stotal$x[1] == stotal$x[2]) {
            # print(sprintf('Tie in slope for %d', f))
            svisits = sort(fam_data$nvisits, index.return=T, decreasing=T)
            if (svisits$x[1] == svisits$x[2]) {
                print(sprintf('Tie in number of visits for %d', f))
                print(fam_data[fam_data$nvisits==svisits$x[1], ]$MRN)
            } else {
                df[fam_rows[svisits$ix[1]], ]$bestInFamily = T
            }
        } else {
            df[fam_rows[stotal$ix[1]], ]$bestInFamily = T
        }
    }
}

df[df$MRN==4585574, ]$bestInFamily = T
df[df$MRN==4925051, ]$bestInFamily = T
df[df$MRN==7079035, ]$bestInFamily = T
df[df$MRN==7378993, ]$bestInFamily = T
# chosen because of overall best MPRAGE QC
df[df$MRN==4640378, ]$bestInFamily = T
# chosen because of overall best MPRAGE QC
df[df$MRN==7218965, ]$bestInFamily = T
```

At this point we have just the df of everyone with prs (393 subjects), and the
selections for best in family. Now, let's merge in the anatomical data.

```r
qc = read.csv('~/data/baseline_prediction/prs_start/prs_and_mprage_qc.csv')
brain_meta = merge(df, qc, by='MRN', all.x=F, all.y=F)

# out of the 1849 scans scored, 1400 are in our prs set

# restrict based on QC
qc_vars = c("mprage_score", "ext_avg", "int_avg")
brain_meta = brain_meta[brain_meta$"age_at_scan" < 18, ]
na_scores_idx = is.na(brain_meta$mprage_score) | is.na(brain_meta$ext_avg) |
                is.na(brain_meta$int_avg)
brain_meta = brain_meta[!na_scores_idx, ]
# down to 1290 after removing adults and anyone without a score

qtile=.95
library(solitude)
iso <- isolationForest$new()
iso$fit(brain_meta[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]
library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(brain_meta[, qc_vars], k = round(.5 * nrow(brain_meta)))
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if
```

<!-- tracts = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
# somehow I have two entries for 1418?
x = duplicated(tracts$id)
jhu_data = merge(dti_meta[idx,], tracts[!x, ], by.x='Mask.ID...Scan', by.y='id')
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])

iso <- isolationForest$new()
iso$fit(jhu_data[, tract_names])
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(jhu_data[, tract_names], k = round(.5 * nrow(jhu_data)))

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

clean_jhu_data = jhu_data[idx, ]

# down to 800 scans when only scans at .95 in both criteria are used

# selecting earliest scan for each subject, regardless of score, as we're assu ing every scan now is good
keep_me = c()
for (s in unique(clean_jhu_data$Medical.Record...MRN...Subjects)) {
    subj_rows = which(clean_jhu_data$Medical.Record...MRN...Subjects == s)
    subj_data = clean_jhu_data[subj_rows, ]
    min_subj_row = which.min(subj_data$age_at_scan...Scan...Subjects)
    keep_me = c(keep_me, subj_rows[min_subj_row])
}
data_dti = clean_jhu_data[keep_me, ]
# finished with 277 scans when using baseline for each subject

b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
# m has all 1020 processed scans
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# dti_meta has all scans for everyone who has PRS data (954 scans)
dti_meta = merge(m, df, by.x="Medical.Record...MRN...Subjects", by.y='MRN',
          all.x=F, all.y=F)

# restrict based on QC
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
dti_meta = dti_meta[dti_meta$"age_at_scan...Scan...Subjects" < 18, ]
dti_meta = dti_meta[dti_meta$"goodVolumes" <= 61, ]
dti_meta = dti_meta[dti_meta$"numVolumes" < 80, ]

# down to 928 scans that obey criteria and have PRS -->

# TODO
* Can we do QC using mriqc parameters?