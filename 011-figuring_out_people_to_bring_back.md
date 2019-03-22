# 2019-03-22 18:23:45

Philip wanted me to check who we should bring back to increase our heritability
numbers. Let's see:

```r
b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# restrict based on QC
pct = m$missingVolumes / m$numVolumes
idx = m$norm.trans < 5 & m$norm.rot < .1 & pct < .15
m = m[idx,]

mBadQC = m[!idx, ]

# down to 902 scans
keep_me = c()
mOneGoodScan = c()
mNoGoodDates = c()
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
        } else {
            mNoGoodDates = rbind(mNoGoodDates, subj_scans)
        }
    } else {
        mOneGoodScan = rbind(mOneGoodScan, subj_scans)
    }
}
m2 = m[keep_me, ]
# down to 566 scans

# let's save everything to visually check later
write.csv(mBadQC, file='~/data/heritability_change/bad_qc.csv',
          row.names=F, na='', quote=F)
write.csv(mOneGoodScan, file='~/data/heritability_change/one_good_scan.csv',
          row.names=F, na='', quote=F)
write.csv(mNoGoodDates, file='~/data/heritability_change/no_good_dates.csv',
          row.names=F, na='', quote=F)
```