# 2018-12-13 15:32:10

Now that we have most of the descriptives results, let's combine them into a
dataset and fire up some ML model. I'll try figuring out what new datasets we
can use for testing later. Or, the second option would be to reduce the current
dataset a bit while the results still hold. But let's try the first approach.

```r
# struct
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
load('/data/NCR_SBRB/baseline_prediction/struct_volume_11142018_260timeDiff12mo.RData.gz')
df = merge(clin, data, by='MRN')
x = colnames(df)[grepl(pattern = '^v_rh', colnames(df))]
a = read.table('~/tmp/clusters.txt')[,1]
idx = which(a==1)
HI_vol_rh = rowMeans(df[, x[idx]])
x = colnames(df)[grepl(pattern = '^v_lh', colnames(df))]
# created new cluster file based on notes
a = read.table('~/tmp/clusters.txt')[,1]
idx = which(a==1)
inatt_vol_lh = rowMeans(df[, x[idx]])
struct = cbind(df$MRN, HI_vol_rh, inatt_vol_lh)
colnames(struct)[1] = 'MRN'
# DTI
load('/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n272_09212018.RData.gz')
a = read.table('~/tmp/out.txt')[,4]
idx = which(a==1)
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
df = merge(clin, data, by='MRN')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
inatt_AD_clu1 = rowMeans(df[, x[idx]])
idx = which(a==2)
inatt_AD_clu2 = rowMeans(df[, x[idx]])
# more cluster work in bash based on notes
load('/data/NCR_SBRB/baseline_prediction/dti_rd_voxelwise_n272_09212018.RData.gz')
a = read.table('~/tmp/out.txt')[,4]
idx = which(a==1)
df = merge(clin, data, by='MRN')
HI_RD_clu1 = rowMeans(df[, x[idx]])
dti = cbind(df$MRN, HI_RD_clu1, inatt_AD_clu1, inatt_AD_clu2)
colnames(dti)[1] = 'MRN'

m = merge(struct, dti, by='MRN', all.x=T, all.y=T)

# now, all the rest
load('/data/NCR_SBRB/baseline_prediction/cog_all_09242018.RData.gz')
m = merge(m, data[, c('MRN', "v_Raw_SS_total", "v_Raw_SSB")], by='MRN', all.x=T, all.y=T)
colnames(m)[7] = 'HI_cog_RawSStotal'
colnames(m)[8] = 'HI_cog_RawSSB'
load('/data/NCR_SBRB/baseline_prediction/adhd200_10042018.RData.gz')
m = merge(m, data[, c('MRN', "v_Age")], by='MRN', all.x=T, all.y=T)
colnames(m)[9] = 'inatt_adhd200_Age'
load('/data/NCR_SBRB/baseline_prediction/clinics_binary_sx_baseline_10022018.RData.gz')
m = merge(m, data[, c('MRN', "v_SX_HI", "vCateg_fidgety", "vCateg_waiting.turn")], by='MRN', all.x=T, all.y=T)
for (i in 10:12) { colnames(m)[i] = sprintf('HI_clinBin_%s', colnames(m)[i]) }
m = merge(m, data[, c('MRN', "v_SX_inatt",          "v_SX_HI",             "vCateg_diff.organ", "vCateg_avoids",       "vCateg_loses",        "vCateg_easily.distr", "vCateg_forgetful","vCateg_waiting.turn")], by='MRN', all.x=T, all.y=T)
for (i in 13:20) { colnames(m)[i] = sprintf('inatt_clinBin_%s', colnames(m)[i]) }
load('/data/NCR_SBRB/baseline_prediction/geno3_prs_09192018.RData.gz')
m = merge(m, data[, c('MRN', "v_PROFILES.0.0001.profile", "v_PROFILES.0.0005.profile")], by='MRN', all.x=T, all.y=F)
for (i in 21:22) { colnames(m)[i] = sprintf('inatt_geno3prs_%s', colnames(m)[i]) }
m = merge(m, data[, c('MRN', "v_PROFILES.0.00001.profile")], by='MRN', all.x=T, all.y=F)
for (i in 23:23) { colnames(m)[i] = sprintf('HI_geno3prs_%s', colnames(m)[i]) }

data=m
save(data, file='/data/NCR_SBRB/baseline_prediction/combined_descriptives_12122018.RData.gz', compressed=T)
```

Now, let's make sure these all still correlate with the targets (without covariates).

```r
library(nlme)
library(lme4)
df = merge(clin, data, by='MRN')
idx2 = df$diag_group != 'new_onset' & df$DX != 'NV'
for (sx in c('HI', 'inatt')) {
    for (i in which(grepl(sprintf('^%s_', sx), colnames(df)))) {
        if (is.factor(df[idx2, i])) {
            fm = as.formula(sprintf("%s ~ OLS_%s_slope + (1|nuclearFamID)", colnames(df)[i], sx))
            print(fm)
            fit = try(glmer(fm, data=df[idx2, ], na.action=na.omit, family = binomial(link = "logit")))
            print(summary(fit)$coefficients[2,])
        } else {
            fm = as.formula(sprintf("%s ~ OLS_%s_slope", colnames(df)[i], sx))
            print(fm)
            fit = try(lme(fm, random=~1|nuclearFamID, data=df[idx2, ], na.action=na.omit))
            print(summary(fit)$tTable[2,])
        }
    }
}
```

Yep, looking good! Just need to add MELODIC clusters later. But for now, we can
start playing with data combination. Or, at least figure out the subjects we can
use for testing.

# 2018-12-17 16:50:21

Let's add the MELODIC clusters.

```r
# struct
clin = read.csv('/data/NCR_SBRB/baseline_prediction/long_clin_11302018.csv')
load('/data/NCR_SBRB/baseline_prediction/melodic_inter_IC31_12142018.RData.gz')
a = read.table('~/tmp/out.txt')[,4]
idx = which(a==1)
df = merge(clin, data, by='MRN')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
inatt_melodic_limbic = rowMeans(df[, x[idx]])
# some cluster work in bash
load('/data/NCR_SBRB/baseline_prediction/melodic_inter_IC2_12142018.RData.gz')
a = read.table('~/tmp/out.txt')[,4]
idx = which(a==1)
df = merge(clin, data, by='MRN')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
inatt_melodic_DMN = rowMeans(df[, x[idx]])
# some cluster work in bash
load('/data/NCR_SBRB/baseline_prediction/melodic_inter_IC11_12142018.RData.gz')
a = read.table('~/tmp/out.txt')[,4]
idx = which(a==1)
df = merge(clin, data, by='MRN')
x = colnames(df)[grepl(pattern = '^v', colnames(df))]
inatt_melodic_VAN = rowMeans(df[, x[idx]])
melodic = cbind(df$MRN, inatt_melodic_limbic, inatt_melodic_DMN,
                inatt_melodic_VAN)
colnames(melodic)[1] = 'MRN'
load('/data/NCR_SBRB/baseline_prediction/combined_descriptives_12122018.RData.gz')
m = merge(data, melodic, by='MRN', all.x=T, all.y=T)
data=m
save(data, file='/data/NCR_SBRB/baseline_prediction/combined_descriptives_12172018.RData.gz', compressed=T)
```

I also re-ran the code ago and the regressions still look good!