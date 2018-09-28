# 2018-09-27 10:51:48

Today I also ran some very crude tests and it doesn't look like movement is
driving the results. The best test will be to pick a model, and then just run it
under different seeds, trying to predict whatever the target was using just the
movement variables, age, and sex. Of course, before we do that we'll need to
decide on the model. But still, this current test gives me hope. Here's some
details:

```r
winsorize = function(x, cut = 0.01){
  cut_point_top <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x, cut, na.rm = T)
  i = which(x >= cut_point_top) 
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom) 
  x[j] = cut_point_bottom
  return(x)
}
# starting h2o
library(h2o)
if (Sys.info()['sysname'] == 'Darwin') {
  max_mem = '16G'
} else {
  max_mem = paste(Sys.getenv('SLURM_MEM_PER_NODE'),'m',sep='')
}
h2o.init(ip='localhost', nthreads=future::availableCores(), max_mem_size=max_mem)
clin_fname='/data/NCR_SBRB/baseline_prediction/long_clin_0918.csv'
data_fname='/data/NCR_SBRB/baseline_prediction/dti_tracts_n272_09212018.RData.gz'
print('Loading files')
# merging phenotype and clinical data
clin = read.csv(clin_fname)
load(data_fname)  #variable is data
# remove constant variables that screw up PCA and univariate tests
print('Removing constant variables')
feat_var = apply(data, 2, var, na.rm=TRUE)
idx = feat_var != 0  # TRUE for features with 0 variance (constant)
# categorical variables give NA variance, but we want to keep them
idx[is.na(idx)] = TRUE
data = data[, idx]
nNAs = colSums(is.na(data))  # number of NAs in each variable
# remove variables that are all NAs
data = data[, nNAs < nrow(data)]
print(sprintf('Features remaining: %d (%d with NAs)', ncol(data)-1, sum(nNAs>0)))
print('Merging files')
df = merge(clin, data, by='MRN')
print('Looking for data columns')
dti = read.csv('/data/NCR_SBRB/baseline_prediction/dti_long_09272018.csv')
m = merge(df, dti, by='mask.id')
df = m
x = c('norm.trans', 'norm.rot', 'goodVolumes', 'age_at_scan', 'Sex')

target='nonew_ADHDonly_diag_group2'
if (grepl('nonew', target)) {
  df = df[df$diag_group != 'new_onset', ]
  df$diag_group = factor(df$diag_group)
  target = sub('nonew_', '', target)
}
if (grepl('ADHDonly', target)) {
  df = df[df$diag_group != 'unaffect', ]
  df$diag_group = factor(df$diag_group)
  target = sub('ADHDonly_', '', target)
}
if (grepl('ADHDNOS', target)) {
  df = df[df$DX != 'NV', ]
  target = sub('ADHDNOS_', '', target)
  if (grepl('groupOLS', target) || grepl('grouprandom', target)) {
    df[, target] = 'nonimprovers'
    slope = sub('group', '', target)
    df[df[, slope] < 0, target] = 'improvers'
    df[, target] = as.factor(df[, target])
  }
}
myseed=42
set.seed(myseed)
idx = sample(1:nrow(df), nrow(df), replace=F)
mylim = floor(.10 * nrow(df))
data.test = df[idx[1:mylim], ]
data.train = df[idx[(mylim + 1):nrow(df)], ]
print(sprintf('Using %d samples for training, %d for testing.',
              nrow(data.train),
              nrow(data.test)))
print('Converting to H2O')
dtrain = as.h2o(data.train[, c(x, target)])
dtest = as.h2o(data.test[, c(x, target)])
if (grepl(pattern = 'group', target)) {
  dtrain[, target] = as.factor(dtrain[, target])
  dtest[, target] = as.factor(dtest[, target])
}
dtrain[, 'Sex'] = as.factor(dtrain[, 'Sex'])
dtest[, 'Sex'] = as.factor(dtest[, 'Sex'])
print(sprintf('Running model on %d features', length(x)))
aml <- h2o.automl(x = x, y = target, training_frame = dtrain,
                  seed=myseed,
                  validation_frame=dtest,
                  max_runtime_secs = NULL,
                  max_models = NULL)
aml@leaderboard
if (grepl(pattern = 'group', target)) {
  print('Class distribution:')
  print(as.vector(h2o.table(dtrain[,target])['Count'])/nrow(dtrain))
} else {
  preds = rep(mean(dtrain[,target]), nrow(dtrain))
  m = h2o.make_metrics(as.h2o(preds), dtrain[, target])
  print('MSE prediction mean:')
  print(m@metrics$MSE)
}
``` 

My best model without age and Sex had .68 AUC, and the class ratio is .69. Now, let's add sex and
age to this. Yep, still only up to .71. Nice.


# 2018-09-28 11:37:18

Let's do the same as above, but now for motion in the rsfmri data:

```r
winsorize = function(x, cut = 0.01){
  cut_point_top <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x, cut, na.rm = T)
  i = which(x >= cut_point_top) 
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom) 
  x[j] = cut_point_bottom
  return(x)
}
# starting h2o
library(h2o)
if (Sys.info()['sysname'] == 'Darwin') {
  max_mem = '16G'
} else {
  max_mem = paste(Sys.getenv('SLURM_MEM_PER_NODE'),'m',sep='')
}
h2o.init(ip='localhost', nthreads=future::availableCores(), max_mem_size=max_mem)
clin_fname='/data/NCR_SBRB/baseline_prediction/long_clin_0918.csv'
data_fname='/data/NCR_SBRB/baseline_prediction/aparc.a2009s_n215_09182018.RData.gz'
print('Loading files')
# merging phenotype and clinical data
clin = read.csv(clin_fname)
load(data_fname)  #variable is data
# remove constant variables that screw up PCA and univariate tests
print('Removing constant variables')
feat_var = apply(data, 2, var, na.rm=TRUE)
idx = feat_var != 0  # TRUE for features with 0 variance (constant)
# categorical variables give NA variance, but we want to keep them
idx[is.na(idx)] = TRUE
data = data[, idx]
nNAs = colSums(is.na(data))  # number of NAs in each variable
# remove variables that are all NAs
data = data[, nNAs < nrow(data)]
print(sprintf('Features remaining: %d (%d with NAs)', ncol(data)-1, sum(nNAs>0)))
print('Merging files')
df = merge(clin, data, by='MRN')
print('Looking for data columns')
library(gdata)
rsfmri = read.xls('/data/NCR_SBRB/baseline_prediction/rsfmri_09282018.xlsx')
rsfmri = rsfmri[rsfmri$baseline_3min,]
m = merge(df, rsfmri, by.x='MRN', by.y='Medical.Record...MRN...Subjects')
df = m
x = c('enorm', 'used_TRs', 'age_at_scan', 'Sex')

target='nonew_ADHDonly_diag_group2'
if (grepl('nonew', target)) {
  df = df[df$diag_group != 'new_onset', ]
  df$diag_group = factor(df$diag_group)
  target = sub('nonew_', '', target)
}
if (grepl('ADHDonly', target)) {
  df = df[df$diag_group != 'unaffect', ]
  df$diag_group = factor(df$diag_group)
  target = sub('ADHDonly_', '', target)
}
if (grepl('ADHDNOS', target)) {
  df = df[df$DX != 'NV', ]
  target = sub('ADHDNOS_', '', target)
  if (grepl('groupOLS', target) || grepl('grouprandom', target)) {
    df[, target] = 'nonimprovers'
    slope = sub('group', '', target)
    df[df[, slope] < 0, target] = 'improvers'
    df[, target] = as.factor(df[, target])
  }
}
myseed=42
set.seed(myseed)
idx = sample(1:nrow(df), nrow(df), replace=F)
mylim = floor(.10 * nrow(df))
data.test = df[idx[1:mylim], ]
data.train = df[idx[(mylim + 1):nrow(df)], ]
print(sprintf('Using %d samples for training, %d for testing.',
              nrow(data.train),
              nrow(data.test)))
print('Converting to H2O')
dtrain = as.h2o(data.train[, c(x, target)])
dtest = as.h2o(data.test[, c(x, target)])
if (grepl(pattern = 'group', target)) {
  dtrain[, target] = as.factor(dtrain[, target])
  dtest[, target] = as.factor(dtest[, target])
}
dtrain[, 'Sex'] = as.factor(dtrain[, 'Sex'])
dtest[, 'Sex'] = as.factor(dtest[, 'Sex'])
print(sprintf('Running model on %d features', length(x)))
aml <- h2o.automl(x = x, y = target, training_frame = dtrain,
                  seed=myseed,
                  validation_frame=dtest,
                  max_runtime_secs = NULL,
                  max_models = NULL)
aml@leaderboard
if (grepl(pattern = 'group', target)) {
  print('Class distribution:')
  print(as.vector(h2o.table(dtrain[,target])['Count'])/nrow(dtrain))
} else {
  preds = rep(mean(dtrain[,target]), nrow(dtrain))
  m = h2o.make_metrics(as.h2o(preds), dtrain[, target])
  print('MSE prediction mean:')
  print(m@metrics$MSE)
}
``` 

For fMRI I got AUC of .58, with class ratio of .63.