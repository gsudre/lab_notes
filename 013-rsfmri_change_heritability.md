# 2019-03-27 16:24:43

I was brainstorming with Philip, and the first idea is to just run the
correlation matrices, calculate the slopes, and try to find heritability there.
The question then becomes what are the vertices in the matrix. We could go for
Freesurfer ROIs, but we might end up getting killed by multiple comparisons.
Another option is to derive the 900 spheres, classify them to belong to
different Yeo networks, and compute within and outside network mean
connectivity. I could also play with some ICA on the overall connectivity maps. Finally, another idea would be to take prototypical templates of different resting state networks in fMRI, do a spatial regressions, and just calculate correlations for those time series (NANs for within-network correlation). Like the Smith et al 10 resting state networks.

Something else I just thought about: we could have voxel to voxel (or sphere to sphere) connectivity matrices, one per scan, then filter down to only connections significant (nominally? FDR?) within scans, then compute slope only between connections significant in both scans. That should give a good amount of filtering, especially across subjects. Another quantity we could assess is proportion of connections still stable? Or positive/negative changes in connections? 

The first step is to check how many datasets we currently have. As usual, we
could use trimmedVSnontrimmed, as well as pearsonVSkendalVSspearman. **These are
all things to try later if the default (Pearson, nontrimmed) doesn't work out.**

For now, while we wait on all the mriqc parameters for our datasets, let's use
simply the number of clean TRs for QC, like we always do.

```r
m = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
m = m[!is.na(m$usedTRs_fmri01), ]
# starting with 1370 scans

# restrict based on QC
minutes = 4
idx = m$usedTRs_fmri01 >= (minutes * 60 / 2.5)
m = m[idx,]
# down to 955 scans

# keep_me = c()
# for (s in unique(m$Medical.Record...MRN...Subjects)) {
#     subj_scans = m[m$Medical.Record...MRN...Subjects==s, ]
#     dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
#                                  format="%m/%d/%Y")
#     if (length(dates) >= 2) {
#         sdates = sort(dates)  # not sure why index.return is not working...
#         # make sure there is at least 6 months between scans
#         next_scan = 2
#         while (((sdates[next_scan] - sdates[1]) < 180) && (next_scan < length(sdates))) {
#             next_scan = next_scan + 1
#         }
#         first_scan_age = subj_scans[dates==sdates[1], 'age_at_scan...Scan...Subjects']
#         if (((sdates[next_scan] - sdates[1]) >= 180) && (first_scan_age < 26)) {
#             idx1 = which(dates == sdates[1])
#             keep_me = c(keep_me, which(m$Mask.ID...Scan == subj_scans[idx1, 'Mask.ID...Scan']))
#             idx2 = which(dates == sdates[next_scan])
#             keep_me = c(keep_me, which(m$Mask.ID...Scan == subj_scans[idx2, 'Mask.ID...Scan']))
#         }
#     }
# }
# m2 = m[keep_me, ]
# # down to 566 scans, as we're not using tract data at this point

# clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
# df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN...Subjects),
#                          x.date='record.date.collected...Scan',
#                          x.id='Medical.Record...MRN...Subjects')
# b = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
# tract_names = colnames(b)[2:ncol(b)]
# df2 = merge(df, b, by.x='Mask.ID...Scan', by.y='id')

# library(MASS)
# mres = df2
# mres$SX_HI = as.numeric(as.character(mres$SX_hi))
# mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
# for (t in tract_names) {
#     print(t)
#     fm_str = sprintf('%s ~', t)
#     fm_str = paste(fm_str,
#                    'norm.rot + I(norm.rot^2) + norm.trans + I(norm.trans^2) +',
#                    'missingVolumes')
#     res.lm <- lm(as.formula(fm_str), data = mres, na.action=na.exclude)
#     step <- stepAIC(res.lm, direction = "both", trace = F)
#     mres[, t] = residuals(step)
# }
# res = c()
# for (s in unique(mres$Medical.Record...MRN...Subjects)) {
#     idx = which(mres$Medical.Record...MRN...Subjects == s)
#     row = c(s, unique(mres[idx, 'Sex...Subjects']))
#     for (t in tract_names) {
#         if (sum(is.na(mres[idx, t])) > 0) {
#             # if any of the tract values is NA, make the slope NA
#             row = c(row, NA)
#         } else {
#             fm_str = sprintf('%s ~ age_at_scan...Scan...Subjects', t)
#            fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
#            row = c(row, coefficients(fit)[2])
#         }
#     }
#     for (t in c('SX_inatt', 'SX_HI')) {
#         fm_str = sprintf('%s ~ age_at_scan...Scan...Subjects', t)
#         fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
#         row = c(row, coefficients(fit)[2])
#     }
#     # grabbing inatt and HI at baseline
#     base_DOA = which.min(mres[idx, 'age_at_scan...Scan...Subjects'])
#     row = c(row, mres[idx[base_DOA], 'SX_inatt'])
#     row = c(row, mres[idx[base_DOA], 'SX_HI'])
#     # DX1 is DSMV definition, DX2 will make SX >=4 as ADHD
#     if (mres[idx[base_DOA], 'age_at_scan...Scan...Subjects'] < 16) {
#         if ((row[length(row)] >= 6) || (row[length(row)-1] >= 6)) {
#             DX = 'ADHD'
#         } else {
#             DX = 'NV'
#         }
#     } else {
#         if ((row[length(row)] >= 5) || (row[length(row)-1] >= 5)) {
#             DX = 'ADHD'
#         } else {
#             DX = 'NV'
#         }
#     }
#     if ((row[length(row)] >= 4) || (row[length(row)-1] >= 4)) {
#         DX2 = 'ADHD'
#     } else {
#         DX2 = 'NV'
#     }
#     row = c(row, DX)
#     row = c(row, DX2)
#     res = rbind(res, row)
# }
# colnames(res) = c('ID', 'sex', tract_names, c('SX_inatt', 'SX_HI',
#                                               'inatt_baseline',
#                                               'HI_baseline', 'DX', 'DX2'))

# # and remove outliers
# for (t in tract_names) {
#     mydata = as.numeric(res[, t])
#     # identifying outliers
#     ul = mean(mydata) + 3 * sd(mydata)
#     ll = mean(mydata) - 3 * sd(mydata)
#     bad_subjs = c(which(mydata < ll), which(mydata > ul))

#     # remove within-tract outliers
#     res[bad_subjs, t] = NA
# }

# # and make sure every family has at least two people
# good_nuclear = names(table(m2$Nuclear.ID...FamilyIDs))[table(m2$Nuclear.ID...FamilyIDs) >= 4]
# good_extended = names(table(m2$Extended.ID...FamilyIDs))[table(m2$Extended.ID...FamilyIDs) >= 4]
# keep_me = c()
# for (f in good_nuclear) {
#     keep_me = c(keep_me, m2[which(m2$Nuclear.ID...FamilyIDs == f),
#                             'Medical.Record...MRN...Subjects'])
# }
# for (f in good_extended) {
#     keep_me = c(keep_me, m2[which(m2$Extended.ID...FamilyIDs == f),
#                             'Medical.Record...MRN...Subjects'])
#     # print(f)
#     # print(m2[which(m2$Extended.ID...FamilyIDs == f),
#     #                         'Medical.Record...MRN...Subjects'])
# }
# keep_me = unique(keep_me)

# fam_subjs = c()
# for (s in keep_me) {
#     fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
# }
# res2 = res[fam_subjs, ]

# res2 = res2[res2[, 'ID'] != 7221745, ]
# write.csv(res2, file='~/data/heritability_change/dti_JHUtracts_residNoSex_OLS_naSlopes133.csv',
#           row.names=F, na='', quote=F)
# ```
# ```