# 2019-02-08 17:33:57

I just had an idea for a new plan, which could work for any new group of
datasets we decide to use, and also include NVs. 

Start by checking which regions of the brain show baseline differences between
NV and ADHD. At first, try ADHDNOS in both groups, and removing new_onset
yes/no. We can then take the biggest X clusters (just to avoid the time
consuming permutation analysis at first), and ask whether they predict which
kids get better (within-ADHD analysis only). For the clusters that do, is there
a relationship to NVs between kids that recover and the ones that don't? Note
that the analysis coud be ML or not.

Then we still have to consider most of the brain, which likely showed no baseline
difference between NVs and ADHDs. Could some of those regions be good predictors
of who gets better? Run the analysis just like before for the entire brain, but
only using ADHDs (and possible ADHD_NOS). 

Whenever we talk about recovery, we could use perVSrem, impVsnonimp,
deterVSmildVSmarked, or just the slope regression.

The first step is to residualize our datasets, just so we don't have to redo it
all the time:

## structural

```r
data_fname = 'struct_volume_11142018_260timeDiff12mo.RData.gz'
base_name = '/data/NCR_SBRB/'
print('Loading files')
# merging phenotype and clinical data
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/%s', base_name, data_fname))  #variable is data
print('Merging files')
df = merge(clin, data, by='MRN')
print('Looking for data columns')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
qc = read.csv(sprintf('%s/baseline_prediction/master_qc.csv', base_name))
df = merge(df, qc, by.x='mask.id', by.y='Mask.ID')
library(gdata)
mprage = read.xls(sprintf('%s/baseline_prediction/long_scans_08072018.xlsx', base_name),
                  sheet='mprage')
df = merge(df, mprage, by.x='mask.id', by.y='Mask.ID...Scan')
library(nlme)
library(MASS)
print(length(x))
mydata = df[, c('Sex...Subjects', 'ext_avg_freesurfer5.3',
                    'int_avg_freesurfer5.3', 'mprage_QC',
                    'age_at_scan', 'nuclearFamID')]
for (v in x) {
    cat(data_fname, v, length(x), '\n')
    mydata$y = df[,v]
    fm = as.formula("y ~ Sex...Subjects + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC + age_at_scan + I(age_at_scan^2)")
    fit = try(lme(fm, random=~1|nuclearFamID, data=mydata, na.action=na.omit, method='ML'))
    if (length(fit) > 1) {
        step = try(stepAIC(fit, direction = "both", trace = F))
        if (length(step) > 1) {
            mydata$y = residuals(step)
        } else {
            mydata$y = residuals(fit)
        }
    }
    df[, v] = mydata$y
}
data = df[, colnames(data)]
out_fname = sprintf('%s/baseline_prediction/resids_%s', base_name, data_fname)
save(data, file=out_fname, compress=T)
```

## dti

```r
data_fname = 'dti_ad_voxelwise_n272_09212018.RData.gz'
base_name = '/data/NCR_SBRB/'
print('Loading files')
# merging phenotype and clinical data
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/%s', base_name, data_fname))  #variable is data
print('Merging files')
df = merge(clin, data, by='MRN')
print('Looking for data columns')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
dti = read.csv(sprintf('%s/baseline_prediction/dti_long_09272018.csv', base_name))
df = merge(df, dti, by='mask.id')

library(nlme)
library(MASS)
print(length(x))
mydata = df[, c('Sex', 'norm.rot', 'norm.trans', 'age_at_scan', 'nuclearFamID')]
for (v in x) {
    cat(data_fname, v, length(x), '\n')
    mydata$y = df[, v]
    fm = as.formula("y ~ Sex + norm.rot + I(norm.rot^2) + norm.trans + I(norm.trans^2) + age_at_scan + I(age_at_scan^2)")
    fit = try(lme(fm, random=~1|nuclearFamID, data=mydata, na.action=na.omit, method='ML'))
    if (length(fit) > 1) {
        step = try(stepAIC(fit, direction = "both", trace = F))
        if (length(step) > 1) {
            mydata$y = residuals(step)
        } else {
            mydata$y = residuals(fit)
        }
    }
    df[, v] = mydata$y
}
data = df[, colnames(data)]
out_fname = sprintf('%s/baseline_prediction/resids_%s', base_name, data_fname)
save(data, file=out_fname, compress=T)
```

## melodic

```r
data_fname = 'melodic_inter_IC11_12142018.RData.gz'
base_name = '/data/NCR_SBRB/'
print('Loading files')
# merging phenotype and clinical data
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/%s', base_name, data_fname))  #variable is data
print('Merging files')
df = merge(clin, data, by='MRN')
print('Looking for data columns')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
qc = read.csv(sprintf('%s/baseline_prediction/master_qc.csv', base_name))
df = merge(df, qc, by.x='mask.id', by.y='Mask.ID')
library(gdata)
mprage = read.xls(sprintf('%s/baseline_prediction/long_scans_08072018.xlsx', base_name),
                  sheet='mprage')
df = merge(df, mprage, by.x='mask.id', by.y='Mask.ID...Scan')
library(nlme)
library(MASS)
print(length(x))
mydata = df[, c('Sex...Subjects', 'enormGoodTRs_fmri01', 'age_at_scan', 'nuclearFamID')]
for (v in x) {
    cat(data_fname, v, length(x), '\n')
    mydata$y = df[,v]
    fm = as.formula("y ~ Sex...Subjects + enormGoodTRs_fmri01 + I(enormGoodTRs_fmri01^2) + age_at_scan + I(age_at_scan^2)")
    fit = try(lme(fm, random=~1|nuclearFamID, data=mydata, na.action=na.omit, method='ML'))
    if (length(fit) > 1) {
        step = try(stepAIC(fit, direction = "both", trace = F))
        if (length(step) > 1) {
            mydata$y = residuals(step)
        } else {
            mydata$y = residuals(fit)
        }
    }
    df[, v] = mydata$y
}
data = df[, colnames(data)]
out_fname = sprintf('%s/baseline_prediction/resids_%s', base_name, data_fname)
save(data, file=out_fname, compress=T)
```

I did all the above, for all 3 neuroimaging domains and their files. It took
overnight, so careful when running it again! I used an interactive node with
several CPUs, but no implicit parallelization.





