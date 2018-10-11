# 2018-10-11 11:35:41

We agreed that I needed to check if the subjects were already remitted or not by
the time of the phenotype acquisition. 

Note that this only applies to the persistent VS remitted phenotype, because
once we start grouping into improvers and non-improvers based on slope, then
there's not really a cut-off anymore.

Starting with fMRI:

```r
load('~/data/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz')
data = data[, !grepl(pattern = '^v', colnames(data))]  # we don't need the actual data
clin = read.csv('~/data/baseline_prediction/long_clin_0918.csv')
df = merge(data, clin, by='MRN')
library(gdata)
scans = read.xls('~/data/baseline_prediction/rsfmri_09182018.xlsx')
scans_base = scans[scans$baseline_any, ]
scans_base$scan_date = as.character(scans_base$record.date.collected...Scan)
df2 = merge(df, scans_base, by.x='MRN', by.y='Medical.Record...MRN...Subjects')

source('~/research_code/lab_mgmt/combine_sx_spreadsheets.R')
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df3 = mergeOnClosestDate(df2, sx, unique(df2$MRN), x.date='scan_date')
df3$cur_group = NA
# 1=NV, 2=remission, 3=persistent
adhd_now = df3$SX_inatt >=6 | df3$SX_hi >= 6
adhd_now[84]=F
adhd_past = df3$SX_inatt_baseline >=6 | df3$SX_HI_baseline >= 6
df3[adhd_past & adhd_now,]$cur_group = 3
df3[adhd_past & !adhd_now,]$cur_group = 2
df3[!adhd_past & !adhd_now,]$cur_group = 1
df3[!adhd_past & adhd_now,]$cur_group = 3 # new onset become persistent
df3[df3$diag_group=='new_onset',]$cur_group = 3 # new onset become persistent
> sum(df3$diag_group2==2 & df3$cur_group==2)
[1] 5
> df3[df3$diag_group2==2 & df3$cur_group==2, 'MRN']
[1] 4439600 4441990 4755674 4836649 7036723
```

So, we have 5 out of the 215 that had already remitted by the time we scanned
them. Should we remove them from the analysis? Let's look at other datasets:

```r
load('~/data/baseline_prediction/dti_ad_voxelwise_n272_09212018.RData.gz')
data = data[, !grepl(pattern = '^v', colnames(data))]  # we don't need the actual data
clin = read.csv('~/data/baseline_prediction/long_clin_0918.csv')
df = merge(data, clin, by='MRN')
library(gdata)
scans = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'dti')
scans$scan_date = as.character(scans$record.date.collected...Scan)
df2 = merge(df, scans, by.x='mask.id', by.y='Mask.ID')

source('~/research_code/lab_mgmt/combine_sx_spreadsheets.R')
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df3 = mergeOnClosestDate(df2, sx, unique(df2$MRN), x.date='scan_date')
df3$cur_group = NA
# 1=NV, 2=remission, 3=persistent
adhd_now = df3$SX_inatt >=6 | df3$SX_hi >= 6
adhd_now[55]=F
adhd_past = df3$SX_inatt_baseline >=6 | df3$SX_HI_baseline >= 6
df3[adhd_past & adhd_now,]$cur_group = 3
df3[adhd_past & !adhd_now,]$cur_group = 2
df3[!adhd_past & !adhd_now,]$cur_group = 1
df3[!adhd_past & adhd_now,]$cur_group = 3 # new onset become persistent
df3[df3$diag_group=='new_onset',]$cur_group = 3 # new onset become persistent
> sum(df3$diag_group2==2 & df3$cur_group==2)
[1] 6
> df3[df3$diag_group2==2 & df3$cur_group==2, 'MRN']
[1] 4439600 4441990 4755674 4915392 7036723 4451521
```

```r
load('~/data/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz')
data = data[, !grepl(pattern = '^v', colnames(data))]  # we don't need the actual data
clin = read.csv('~/data/baseline_prediction/long_clin_0918.csv')
df = merge(data, clin, by='MRN')
library(gdata)
scans = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'dti')
scans$scan_date = as.character(scans$record.date.collected...Scan)
df2 = merge(df, scans, by.x='mask.id', by.y='Mask.ID')

source('~/research_code/lab_mgmt/combine_sx_spreadsheets.R')
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df3 = mergeOnClosestDate(df2, sx, unique(df2$MRN), x.date='scan_date')
df3$cur_group = NA
# 1=NV, 2=remission, 3=persistent
adhd_now = df3$SX_inatt >=6 | df3$SX_hi >= 6
adhd_now[27]=F
adhd_past = df3$SX_inatt_baseline >=6 | df3$SX_HI_baseline >= 6
df3[adhd_past & adhd_now,]$cur_group = 3
df3[adhd_past & !adhd_now,]$cur_group = 2
df3[!adhd_past & !adhd_now,]$cur_group = 1
df3[!adhd_past & adhd_now,]$cur_group = 3 # new onset become persistent
df3[df3$diag_group=='new_onset',]$cur_group = 3 # new onset become persistent
> sum(df3$diag_group2==2 & df3$cur_group==2)
[1] 4
> df3[df3$diag_group2==2 & df3$cur_group==2, 'MRN']
[1] 4441990 4836649 4915392 4451521
```

So, the numbers are still quite close, and it's the same subjects across
modalities. Just for confirmation, let's look at structural:

```r
load('~/data/baseline_prediction/struct_rois_09192018_260timeDiff12mo.RData.gz')
data = data[, !grepl(pattern = '^v', colnames(data))]  # we don't need the actual data
clin = read.csv('~/data/baseline_prediction/long_clin_0918.csv')
df = merge(data, clin, by='MRN')
library(gdata)
scans = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'mprage')
scans$scan_date = as.character(scans$record.date.collected...Scan)
df2 = merge(df, scans, by.x='mask.id', by.y='Mask.ID...Scan')

source('~/research_code/lab_mgmt/combine_sx_spreadsheets.R')
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df3 = mergeOnClosestDate(df2, sx, unique(df2$MRN), x.date='scan_date')
df3$cur_group = NA
# 1=NV, 2=remission, 3=persistent
adhd_now = df3$SX_inatt >=6 | df3$SX_hi >= 6
adhd_now[3]=F
adhd_now[66]=F
adhd_past = df3$SX_inatt_baseline >=6 | df3$SX_HI_baseline >= 6
df3[adhd_past & adhd_now,]$cur_group = 3
df3[adhd_past & !adhd_now,]$cur_group = 2
df3[!adhd_past & !adhd_now,]$cur_group = 1
df3[!adhd_past & adhd_now,]$cur_group = 3 # new onset become persistent
df3[df3$diag_group=='new_onset',]$cur_group = 3 # new onset become persistent
> sum(df3$diag_group2==2 & df3$cur_group==2)
[1] 6
> df3[df3$diag_group2==2 & df3$cur_group==2, 'MRN']
[1] 4432605 4441990 4439600 4975406 4836649 4915392
```

So, it's good to know these numbers, but let's not do anything about it now.
Mostly because it could just be a robustness analysis to remove them later, or
in a different scenario they'd be used anyways if we split the groups into
improvers and non-improvers.