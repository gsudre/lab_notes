library(RcppCNPy)
gf = read.csv('~/data/prs/fmri_processed_prs_unique_08092017.csv')
idx = gf$all=='TRUE'
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_all_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_all', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(gf, mydata, by.x=1, by.y=1, all.x=T)
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_all_trimmed_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_all_trimmed', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)

idx = gf$X2min=='TRUE'
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_2min_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_2min', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_2min_trimmed_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_2min_trimmed', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)

idx = gf$X3min=='TRUE'
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_3min_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_3min', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_3min_trimmed_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_3min_trimmed', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)

idx = gf$X4min=='TRUE'
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_4min_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_4min', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_4min_trimmed_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_4min_trimmed', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)

idx = gf$elbow..110.TRs.=='TRUE'
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_elbow_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_elbow', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)
fmat <- npyLoad("/Users/sudregp/data/prs/exp_scores_elbow_trimmed_1Kperms_15ics.npy")
colnames(fmat) = lapply(1:15, function(d) sprintf('IC%02d_elbow_trimmed', d))
mydata = cbind(gf[idx,]$Medical.Record...MRN...Subjects, fmat)
colnames(mydata)[1] = 'mrn'
all_data = merge(all_data, mydata, by.x=1, by.y=1, all.x=T)

write.csv(all_data, file='~/data/prs/resting_08112017.csv', row.names=F)