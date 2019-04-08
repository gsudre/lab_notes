# 2019-04-08 10:54:25

A few notes on how I merged the DLPFC big files for Philip. Basically, the idea is to
create a "skinny" version of merge. It looks for the IDs that are in both files,
then grabs the indexes for IDs in each matrix, slice them, and just do a cbind.
We finally rename columns to know what came from where.

```r
library(data.table)
grex = fread('~/philip/GREX_COMMON.csv', header = T, sep = ',')
pc = read.csv('~/philip/PC_common.csv')
brain = read.csv('~/philip/imaging_combined_short_for_gustavo.csv')
ages = read.csv('~/philip/Age_gender_common.csv')
brain2 = read.csv('~/philip/NHGROI_good_qc_choosing_anatomic_scan_one_perfamily.csv')

ids1 = as.character(grex$ID)
ids1[is.na(grex$ID)] = as.character(grex[is.na(grex$ID),]$ID2)
ids2 = as.character(pc$ID)
ids2[is.na(pc$ID)] = as.character(pc[is.na(pc$ID),]$ID2)
in_both = intersect(unique(ids1), unique(ids2))
ninter = length(in_both)
idx1 = vector(length=ninter)
idx2 = vector(length=ninter)
for (i in 1:ninter) {
    idx1[i] = which(grepl(ids1, pattern=in_both[i]))
    idx2[i] = which(grepl(ids2, pattern=in_both[i]))
}
m1 = cbind(grex[idx1, ], pc[idx2, ])
cnames = sapply(colnames(grex), function(x) sprintf('grex.%s', x))
colnames(m1)[1:ncol(grex)] = cnames
cnames = sapply(colnames(pc), function(x) sprintf('pc.%s', x))
colnames(m1)[(ncol(grex)+1):(ncol(grex)+ncol(pc))] = cnames
m1 = cbind(m1, ids1[idx1])
colnames(m1)[ncol(m1)] = 'ID'
write.csv(m1, file='~/philip/grex_pc.csv', row.names=F, quote=F)

ncrbrain = brain[brain$cohort=='nhgri', ]
b = merge(ncrbrain, brain2, by.x='ID2', by.y='MASKID')
library(plyr)
b2 = rbind.fill(brain[brain$cohort!='nhgri', ], b)
ids3 = as.character(b2$MRN)
ids3[b2$cohort=='genr'] = as.character(b2[b2$cohort=='genr',]$ID)
ids3[b2$cohort=='pnc'] = as.character(b2[b2$cohort=='pnc',]$ID2)


ids3 = as.character(ages$ID)
in_both = intersect(unique(m1$ID), unique(ids3))
ninter = length(in_both)
idx1 = vector(length=ninter)
idx2 = vector(length=ninter)
for (i in 1:ninter) {
    idx1[i] = which(grepl(ids1, pattern=in_both[i]))
    idx2[i] = which(grepl(ids2, pattern=in_both[i]))
}


```