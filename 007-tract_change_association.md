# 2019-03-13 09:19:01

This is a complement to the 005 note, but here we're using the entire cohort to
get tract associations... let's see if it helps.

First, let's run through DTI-TK anyone that we can and would potentially help
the analysis.

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190308.csv')
a = read.csv('~/data/heritability_change/ready_1045.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# keep only subjects with more than one scan processed in DTITK
keep_me = c()
for (s in unique(m$Medical.Record...MRN)) {
    if (sum(m$Medical.Record...MRN == s) > 1) {
        keep_me = c(keep_me, which(m$Medical.Record...MRN == s))
    }
}
m = m[keep_me, ]
# down to 1021 because not all 1045 had two scans ready form TORTOISE 

# restrict based on QC
m = m[!is.na(m$norm.trans), ]
idx = (m$mean_fa < (mean(m$mean_fa)+3*sd(m$mean_fa)) &
       m$mean_fa > (mean(m$mean_fa)-3*sd(m$mean_fa)) &
       m$mean_ad < (mean(m$mean_ad)+3*sd(m$mean_ad)) &
       m$mean_ad > (mean(m$mean_ad)-3*sd(m$mean_ad)) &
       m$mean_rd < (mean(m$mean_rd)+3*sd(m$mean_rd)) &
       m$mean_rd > (mean(m$mean_rd)-3*sd(m$mean_rd)))
m = m[idx,]
# down to 1009

# select only the first two for each person (younger than 26 y.o.)
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
# down to 624 scans
```

Now we need to merge those good scans with their clinical data:

```r
clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN...Subjects),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN...Subjects')
b = read.csv('~/data/heritability_change/dti_mean_phenotype_624.csv')
tract_names = colnames(b)[2:ncol(b)]
df2 = merge(df, b, by.x='Mask.ID...Scan', by.y='file')
```

Time to residualize the data and compute the slopes, including SX slopes.

```r
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
    res.lm <- lm(as.formula(fm_str), data = mres)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    mres[, t] = residuals(step)
}
res = c()
for (s in unique(mres$Medical.Record...MRN...Subjects)) {
    idx = which(mres$Medical.Record...MRN...Subjects == s)
    row = c(s, unique(mres[idx, 'Sex...Subjects']))
    for (t in tract_names) {
        fm_str = sprintf('%s ~ age_at_scan...Scan...Subjects', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ])
        row = c(row, coefficients(fit)[2])
    }
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan...Scan...Subjects', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ])
        row = c(row, coefficients(fit)[2])
    }
    base_DOA = which.min(mres[idx, 'age_at_scan...Scan...Subjects'])
    row = c(row, mres[idx[base_DOA], 'SX_inatt'])
    row = c(row, mres[idx[base_DOA], 'SX_HI'])
    res = rbind(res, row)
}
colnames(res) = c('ID', 'sex', tract_names, c('SX_inatt', 'SX_HI',
                                              'inatt_baseline',
                                              'HI_baseline'))
# write.csv(res, file='~/data/heritability_change/dti_tracts_residNoSex_OLS_152.csv',
#           row.names=F, quote=F)

```