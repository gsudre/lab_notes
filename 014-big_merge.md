# 2019-04-08 10:54:25

A few notes on how I merged the DLPFC big files for Philip. Basically, the idea is to
create a "skinny" version of merge. It looks for the IDs that are in both files,
then grabs the indexes for IDs in each matrix, slice them, and just do a cbind.
We finally rename columns to know what came from where.

```r
library(data.table)
grex = fread('~/philip/GREX_COMMON-2.csv', header = T, sep = ',')
pc = read.csv('~/philip/PC_common.csv')
brain = read.csv('~/philip/imaging_combined_short_for_gustavo.csv')
ages = read.csv('~/philip/Age_gender_common.csv')
brain2 = read.csv('~/philip/NHGROI_good_qc_choosing_anatomic_scan_one_perfamily.csv')
pnc_code = read.csv('~/philip/family_ID.csv')

# it will be easy to split this by cohorts
cohort = 'genr'
idx = grex$cohort==cohort  # had to do it this way... otherwise it wouldn't work!
grex.g = grex[idx,]
pc.g = pc[pc$cohort==cohort,]
brain.g = brain[brain$cohort==cohort,]
ages.g = ages[ages$cohort==cohort, c('ID', 'AGE_ASSESSMENT_USE', 'Sex.y')]
colnames(ages.g) = c('ID2',  'AGE_ASSESSMENT_USE', 'Sex')
m1.g = merge(grex.g, pc.g, by='ID2')
m1.g = merge(m1.g, ages.g, by='ID2', all.x=T)
m2.g = merge(m1.g, brain.g, by.x='ID2', by.y='ID')

cohort = 'pnc'
idx = grex$cohort==cohort  # had to do it this way... otherwise it wouldn't work!
grex.p = grex[idx,]
pc.p = pc[pc$cohort==cohort,]
brain.p = brain[brain$cohort==cohort,]
ages.p = ages[ages$cohort==cohort, c('ID', 'AGE_ASSESSMENT_USE', 'Sex.y')]
colnames(ages.p) = c('ID',  'AGE_ASSESSMENT_USE', 'Sex')
# merge it here while matrices are still manageable
pc.p = merge(pc.p, ages.p, by='ID')
# we'll need a skinny merge here, because there's too much data for a regular merge
ids1 = as.character(grex.p$ID)
ids2 = as.character(pc.p$ID)
in_both = intersect(unique(ids1), unique(ids2))
ninter = length(in_both)
idx1 = vector(length=ninter)
idx2 = vector(length=ninter)
for (i in 1:ninter) {
    idx1[i] = which(grepl(ids1, pattern=in_both[i]))
    idx2[i] = which(grepl(ids2, pattern=in_both[i]))
}
m1.p = cbind(grex.p[idx1, ], pc.p[idx2, ])
ids3 = as.character(brain.p$dbGaP_Subject_ID.x)
in_both2 = intersect(in_both, unique(ids3))
ninter = length(in_both2)
idx1 = vector(length=ninter)
idx2 = vector(length=ninter)
for (i in 1:ninter) {
    idx1[i] = which(grepl(in_both, pattern=in_both2[i]))
    idx2[i] = which(grepl(ids3, pattern=in_both2[i]))
}
m2.p = cbind(m1.p[idx1, ], brain.p[idx2, ])

cohort = 'nhgri'
idx = grex$cohort==cohort  # had to do it this way... otherwise it wouldn't work!
grex.n = grex[idx,]
pc.n = pc[pc$cohort==cohort,]
ages.n = ages[ages$cohort==cohort, c('ID', 'AGE_ASSESSMENT_USE', 'Sex.y')]
colnames(ages.n) = c('ID',  'AGE_ASSESSMENT_USE', 'Sex')
pc.n = merge(pc.n, ages.n, by='ID')
m1.n = merge(grex.n, pc.n, by='ID')
brain.n = brain[brain$cohort==cohort,]
brain.n = merge(brain.n, brain2, by.x='ID2', by.y='MASKID')
m2.n = merge(m1.n, brain.n, by.x='ID', by.y='MRN')

# now we merge across cohorts both m1 and m2
m1.g$mergeID = as.character(m1.g$ID2)
m2.g$mergeID = as.character(m2.g$ID2)
m1.p$mergeID = sapply(in_both, function(x) sprintf('P%s', x))
m2.p$mergeID = sapply(in_both2, function(x) sprintf('P%s', x))
m1.n$mergeID = sapply(m1.n$ID, function(x) sprintf('N%d', x))
m2.n$mergeID = sapply(m2.n$ID, function(x) sprintf('N%d', x))

# special plyr function to autofill with NAs columns that only exist in one matrix
library(plyr)
m1 = rbind.fill(m1.g, m1.n)
m1 = rbind.fill(m1, m1.p)
write.csv(m1, file='~/philip/grex_pc.csv', row.names=F, quote=F)

m2 = rbind.fill(m2.g, m2.n)
m2 = rbind.fill(m2, m2.p)
write.csv(m2, file='~/philip/grex_pc_anat.csv', row.names=F, quote=F)
```