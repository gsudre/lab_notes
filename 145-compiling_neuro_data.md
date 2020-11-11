# 2020-11-09 20:27:15

Since we had this long break because of COVID, all Freesurfer and DTI data
finished being processed. It makes sense to create a big csv with all the data
in it.

## MPRAGE

We start from Labmatrix to get the best mprage QC, SID, and age_scan. We also
make sure to only get the scans that were processed.

Then:

```bash
# ncrshell01
bash ~/research_code/lab_mgmt/get_freesurfer_roi_data.sh ~/tmp/ids10.txt
```

Then we organize everything in Excel, and to create the summary metrics we do:

```r
rois = read.csv('~/research_Code/REGIONAL_ANALYSES_FREESURFER.csv')
brain_data = read.csv('~/data/all_freesurfer.csv')
brain_vars = colnames(brain_data)[14:287]
for (part in c('sublobar', 'lobar', 'theoryDriven')) {
    for (roi in unique(rois[, part])) {
        labels = rois[which(rois[, part]==roi), 'region']
        to_avg = c()
        for (l in labels) {
            to_avg = c(to_avg,
                    brain_vars[grepl(brain_vars, pattern=sprintf("^%s", l))])
        }
        # only use variable if it's selected initially and defined
        if (length(to_avg) > 0 && sum(is.na(brain_data[, to_avg])) == 0 &&
            nchar(roi) > 0) {
            if (length(to_avg) == 1) {
                brain_data[, sprintf('%s_%s', part, roi)] = brain_data[, to_avg]
            } else {
                brain_data[, sprintf('%s_%s', part, roi)] = rowMeans(brain_data[, to_avg])
            }
        }
    }
}
write.csv(brain_data, file='~/data/all_freesurfer_11122020.csv', row.names=F, quote=F)
```

## DTI


