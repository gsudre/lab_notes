# 2020-02-03 08:53:18

At this point, we have a set of variables we want to include in the big linear
model. Even though the actual set of variables is still up for debate (see note
74), I have an idea of where it's going. 

Now, the idea is to create two pie charts, one for all group prediction and one
for improvers VS nonimproves. Each pie chart represents each predictor's
contribution to the overall model. This is a nice tutorial using the dominance
package:

https://cran.r-project.org/web/packages/dominanceanalysis/vignettes/da-logistic-regression.html

Let's first make sure it works with a mixel effect logistic regression, and then
also with a multinomial logistic regression. I might end up having to cut that
to either keep it to binomial and/or dropping the mixed effect, if the model
doesn't converge. But the first question is wheter or not we should impute our
data.

```r
setwd('~/data/baseline_prediction/prs_start/')
clin_long = read.csv('long_clin_01062020_lt16.csv')
clin_long$SX_total = clin_long$SX_inatt + clin_long$SX_hi

winsorize = function(x, cut = 0.01){
  cut_point_top <- quantile(x, 1 - cut, na.rm = T)
  cut_point_bottom <- quantile(x, cut, na.rm = T)
  i = which(x >= cut_point_top) 
  x[i] = cut_point_top
  j = which(x <= cut_point_bottom) 
  x[j] = cut_point_bottom
  return(x)
}

df = data.frame(MRN=unique(clin_long$MRN))
for (r in 1:nrow(df)) {
    subj_data = clin_long[clin_long$MRN==df$MRN[r], ]
    for (sx in c('inatt', 'hi', 'total')) {
        fit = lm(as.formula(sprintf('SX_%s ~ age', sx)), data=subj_data)
        df[r, sprintf('slope_%s', sx)] = fit$coefficients['age']
        base_row = which.min(subj_data$age)
        df[r, sprintf('base_%s', sx)] = subj_data[base_row, sprintf('SX_%s', sx)]
        last_row = which.max(subj_data$age)
        df[r, sprintf('last_%s', sx)] = subj_data[last_row, sprintf('SX_%s', sx)]
        df[r, 'base_age'] = subj_data[base_row, 'age']
        df[r, 'base_DOA'] = subj_data[base_row, 'DOA']
        df[r, 'last_age'] = subj_data[last_row, 'age']
        df[r, 'sex'] = subj_data[last_row, 'sex']
    }
}
for (min_sx in c(0, 3, 4, 6)) {
    idx = df$base_inatt>=min_sx | df$base_hi>=min_sx
    for (sx in c('inatt', 'hi', 'total')) {
        df[, sprintf('slope_%s_GE%d_wp05', sx, min_sx)] = NA
        junk = winsorize(df[idx, sprintf('slope_%s', sx)], cut=.05)
        df[idx, sprintf('slope_%s_GE%d_wp05', sx, min_sx)] = junk
    }
}

demo = read.csv('prs_demo.csv')
# just to get FAMID, sex already there
df = merge(df, subset(demo, select=-sex), by='MRN')

# selecting best kid in family
df$bestInFamily = F
nvisits = table(clin_long$MRN)
df = merge(df, as.matrix(nvisits),
                 by.x='MRN', by.y=0)
colnames(df)[ncol(df)] = 'nvisits'
for (f in unique(df$FAMID)) {
    fam_rows = which(df$FAMID == f)
    fam_data = df[fam_rows,]
    if (nrow(fam_data) == 1) {
        df[fam_rows,]$bestInFamily = T
    } else {
        stotal = sort(fam_data$slope_total, index.return=T, decreasing=T)
        # if there's a tie
        if (stotal$x[1] == stotal$x[2]) {
            # print(sprintf('Tie in slope for %d', f))
            svisits = sort(fam_data$nvisits, index.return=T, decreasing=T)
            if (svisits$x[1] == svisits$x[2]) {
                print(sprintf('Tie in number of visits for %d', f))
                print(fam_data[fam_data$nvisits==svisits$x[1], ]$MRN)
            } else {
                df[fam_rows[svisits$ix[1]], ]$bestInFamily = T
            }
        } else {
            df[fam_rows[stotal$ix[1]], ]$bestInFamily = T
        }
    }
}

df[df$MRN==4585574, ]$bestInFamily = T
df[df$MRN==4925051, ]$bestInFamily = T
df[df$MRN==7079035, ]$bestInFamily = T
df[df$MRN==7378993, ]$bestInFamily = T
# chosen because of overall best MPRAGE QC
df[df$MRN==4640378, ]$bestInFamily = T
# chosen because of overall best MPRAGE QC
df[df$MRN==7218965, ]$bestInFamily = T

prs = read.csv('/Volumes/NCR/reference/merged_NCR_1KG_PRS_12192019.csv')
data = merge(df, prs, by='MRN', all.x=F, all.y=F)

library(nlme)
library(MASS)

covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
sx = 'hi'
min_sx = 6
bv = 'ADHD_PRS0.001000'
if (sx == 'inatt') {
    thresh = 0
} else if (sx == 'hi') {
    thresh = -.5
}
phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data[, phen] = NA
data[which(data[, phen_slope] < thresh), phen] = 'imp'
data[which(data[, phen_slope] >= thresh), phen] = 'nonimp'
data[, phen] = factor(data[, phen], ordered=F)
data[, phen] = relevel(data[, phen], ref='imp')
use_me = !is.na(data[, phen])

this_data = data[use_me, c(phen, 'FAMID', bv, covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$sex = data[use_me, 'sex']
tmp_covars = c(covars, 'sex')
fm_str = paste(bv, " ~ ordered +",
                    paste(tmp_covars, collapse='+'),
                    sep="")
fit = lme(as.formula(fm_str), ~1|FAMID, data=this_data, method='ML')
step = stepAIC(fit, direction='both', trace=F,
                scope = list(lower = ~ ordered))


fm_str = paste(phen, " ~ ", bv, ' + ',
               paste(tmp_covars, collapse='+'), ' + (1|FAMID)',
               sep="")
fit = glmer(as.formula(fm_str), data=this_data,
            family=binomial(link='logit'))
```

I'm getting singular fit in impVSnonimp even for PRS, which is the condition
where we'd have most data. So, I cannot run lmer there... I'll have to either
run one kid per family, or assume the family term doesn't contribute that much.
But since I cannot even run the model, that's difficult. Let me see if a
multinomial logistic mixed model runs. I can compare the contribution of the
family term that way. 

Or I can run just the model to predict the binomial case for best in family,
then for everyone, and show that the results don't change that much.

```r
fm_str = paste(phen, " ~ ", bv, ' + ',
               paste(tmp_covars, collapse='+'),
               sep="")
tmp_covars = c(covars, 'sex')

use_me = !is.na(data[, phen])
this_data = data[use_me, c(phen, bv, covars)]
this_data[, 2:ncol(this_data)] = scale(this_data[, 2:ncol(this_data)])
this_data$sex = data[use_me, 'sex']
fit_all = glm(as.formula(fm_str), data=this_data, family=binomial(link='logit'))

use_me = !is.na(data[, phen]) & data$bestInFamily
this_data = data[use_me, c(phen, bv, covars)]
this_data[, 2:ncol(this_data)] = scale(this_data[, 2:ncol(this_data)])
this_data$sex = data[use_me, 'sex']
tmp_covars = c(covars, 'sex')
fit_bif = glm(as.formula(fm_str), data=this_data, family=binomial(link='logit'))

# make sure to step to remove only covariates!
step = stepAIC(fit, direction='both', trace=F,
                scope = list(lower = as.formula(sprintf('~ %s', bv))))
```

So, at this point we almost have a 4x4 matrix of conditions: imputedVSnonimputed
and bestInFamilyVSall. In all cases, we still need to decide on the predictors,
but we'll be running a binomial logistic regression only. 

Let's merge everything keeping NAs, and we can decide it later. It's also not a
bad idea to save this matrix for future use...

```r
qc = read.csv('~/data/baseline_prediction/prs_start/prs_and_mprage_qc.csv')
brain_meta = merge(df, qc, by='MRN', all.x=F, all.y=F)
qc_vars = c("mprage_score", "ext_avg", "int_avg")
brain_meta = brain_meta[brain_meta$"age_at_scan" < 18, ]
na_scores_idx = is.na(brain_meta$mprage_score) | is.na(brain_meta$ext_avg) |
                is.na(brain_meta$int_avg)
brain_meta = brain_meta[!na_scores_idx, ]
qtile=.95
library(solitude)
iso <- isolationForest$new()
iso$fit(brain_meta[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]
library(dbscan)
scores_lof = lof(brain_meta[, qc_vars], k = round(.5 * nrow(brain_meta)))
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if
all_brain_data = read.table('~/data/baseline_prediction/merged_rois.txt', header=T)
x = duplicated(all_brain_data$lh.aparc.area)
brain_data = merge(brain_meta[idx,], all_brain_data[!x, ], by.x='maskid',
                   by.y='lh.aparc.area', all.x=F, all.y=F)
rois = read.csv('~/data/baseline_prediction/REGIONAL_ANALYSES_FREESURFER.csv')
brain_vars = colnames(brain_data)[grepl(colnames(brain_data), pattern="_thickness$")]
part = 'lobar'
new_brain_vars = c()
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
            brain_data[, roi] = brain_data[, to_avg]
        } else {
            brain_data[, roi] = rowMeans(brain_data[, to_avg])
        }
        new_brain_vars = c(new_brain_vars, roi)
    }
}
brain_vars = new_brain_vars
iso <- isolationForest$new()
iso$fit(brain_data[, brain_vars])
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(brain_data[, brain_vars], k = round(.5 * nrow(brain_data)))
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if
clean_brain_data = brain_data[idx, ]
keep_me = c()
for (s in unique(clean_brain_data$MRN)) {
    subj_rows = which(clean_brain_data$MRN == s)
    subj_data = clean_brain_data[subj_rows, ]
    min_subj_row = which.min(subj_data$age_at_scan)
    if (abs(subj_data[min_subj_row, 'base_age'] -
            subj_data[min_subj_row, 'age_at_scan'])<1) {
        keep_me = c(keep_me, subj_rows[min_subj_row])
    }
}
anat_data = clean_brain_data[keep_me, ]
data = merge(data, anat_data[, c('MRN', qc_vars, brain_vars)], by='MRN',
             all.x=T, all.y=F)
```

And we do the same with DTI:

```r
qc = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
brain_demo = read.csv('~/data/heritability_change/ready_1020.csv')
m = merge(brain_demo, qc, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)
brain_meta = merge(m, df, by.x="Medical.Record...MRN...Subjects", by.y='MRN',
                   all.x=F, all.y=F)
colnames(brain_meta)[1] = 'MRN'
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
brain_meta = brain_meta[brain_meta$"age_at_scan...Scan...Subjects" < 18, ]
brain_meta = brain_meta[brain_meta$"goodVolumes" <= 61, ]
brain_meta = brain_meta[brain_meta$"numVolumes" < 80, ]
qtile=.95
library(solitude)
iso <- isolationForest$new()
iso$fit(brain_meta[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]
library(dbscan)
scores_lof = lof(brain_meta[, qc_vars], k = round(.5 * nrow(brain_meta)))
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if
all_brain_data = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
x = duplicated(all_brain_data$id)
brain_data = merge(brain_meta[idx,], all_brain_data[!x, ], by.x='Mask.ID...Scan',
                   by.y='id')
all_brain_data = read.csv('~/data/baseline_prediction/jhu_tracts_mode.csv')
x = duplicated(all_brain_data$id)
brain_data = merge(brain_data, all_brain_data[!x, ], by.x='Mask.ID...Scan',
                   by.y='id')
for (p in c('fa', 'ad', 'rd', 'mode')) {
    brain_data[, sprintf('ATR_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_1', p),
                                                                 sprintf('%s_2', p))])
    brain_data[, sprintf('CST_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_3', p),
                                                                 sprintf('%s_4', p))])
    brain_data[, sprintf('CIN_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_5', p),
                                                                 sprintf('%s_6', p),
                                                                 sprintf('%s_7', p),
                                                                 sprintf('%s_8', p))])
    brain_data[, sprintf('CC_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_9', p),
                                                                 sprintf('%s_10', p))])
    brain_data[, sprintf('IFO_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_11', p),
                                                                 sprintf('%s_12', p))])
    brain_data[, sprintf('ILF_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_13', p),
                                                                 sprintf('%s_14', p))])
    brain_data[, sprintf('SLF_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_15', p),
                                                                 sprintf('%s_16', p),
                                                                 sprintf('%s_19', p),
                                                                 sprintf('%s_20', p))])
    brain_data[, sprintf('UNC_%s', p)] = rowMeans(brain_data[, c(sprintf('%s_17', p),
                                                                 sprintf('%s_18', p))])
}
brain_vars = colnames(brain_data)[grepl(colnames(brain_data), pattern="_fa$") &
                                  !grepl(colnames(brain_data), pattern="^mean")]
iso <- isolationForest$new()
iso$fit(brain_data[, brain_vars])
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(brain_data[, brain_vars], k = round(.5 * nrow(brain_data)))
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if
clean_brain_data = brain_data[idx, ]
keep_me = c()
for (s in unique(clean_brain_data$MRN)) {
    subj_rows = which(clean_brain_data$MRN == s)
    subj_data = clean_brain_data[subj_rows, ]
    min_subj_row = which.min(subj_data$age_at_scan...Scan...Subjects)
    if (abs(subj_data[min_subj_row, 'base_age'] -
            subj_data[min_subj_row, 'age_at_scan...Scan...Subjects'])<1) {
        keep_me = c(keep_me, subj_rows[min_subj_row])
    }
}
dti_data = clean_brain_data[keep_me, ]
data = merge(data, dti_data[, c('MRN', qc_vars, brain_vars)], by='MRN',
             all.x=T, all.y=F)
```

Finally, let's merge in the neuropsych and the clinical remaining data:

```r
iq = read.csv('~/data/baseline_prediction/basics.csv')
# no need to curb to 1 year difference
data = merge(data, iq, by='MRN', all.x=T, all.y=F)

library(gdata)
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
beery = read.xls('~/data/baseline_prediction/prs_start/Subjects_Beery_clean.xlsx')
colnames(beery) = c('MRN', 'DOA.beery', 'VMI.beery')
neuropsych = mergeOnClosestDate(df, beery, unique(df$MRN), x.date='base_DOA',
                                y.date='DOA.beery')
colnames(neuropsych)[ncol(neuropsych)] = 'dateDiff.beery'
wisc = read.xls('~/data/baseline_prediction/prs_start/Subjects_WISC_clean.xlsx')
colnames(wisc) = c('MRN', 'DOA.wisc', 'DS.wisc', 'DSB.wisc', 'DSF.wisc', 'SS.wisc',
                   'SSB.wisc', 'SSF.wisc')
neuropsych = mergeOnClosestDate(neuropsych, wisc, unique(df$MRN),
                                x.date='base_DOA', y.date='DOA.wisc')
colnames(neuropsych)[ncol(neuropsych)] = 'dateDiff.wisc'
wj = read.xls('~/data/baseline_prediction/prs_start/Subjects_Woodcock_Johnson_clean.xlsx')
colnames(wj) = c('MRN', 'DOA.wj', 'PS.wj', 'DS.wj', 'VM.wj')
neuropsych = mergeOnClosestDate(neuropsych, wj, unique(df$MRN),
                                x.date='base_DOA', y.date='DOA.wj')
colnames(neuropsych)[ncol(neuropsych)] = 'dateDiff.wj'
for (suf in c('.beery', '.wisc', '.wj')) {
    doa_col = sprintf('DOA%s', suf)
    date_diff = abs(as.Date(neuropsych[, 'base_DOA'], tryFormats='%m/%d/%Y') -
                    as.Date(neuropsych[, doa_col], tryFormats='%m/%d/%Y'))
    mycols = colnames(neuropsych)[grepl(colnames(neuropsych), pattern=sprintf('%s$', suf))]
    idx = which(date_diff > 365)
    neuropsych[idx, mycols] = NA
}
brain_vars = c('VMI.beery',
               'DS.wisc', 'DSB.wisc', 'DSF.wisc', 'SS.wisc', 'SSB.wisc', 'SSF.wisc',
               'PS.wj', 'DS.wj', 'VM.wj')
data = merge(data, neuropsych[, c('MRN', brain_vars)], by='MRN', all.x=T, all.y=F)
saveRDS(data, file='~/data/baseline_prediction/prs_start/complete_massaged_data_02032020.rds', compress=T)
```

Now it's a matter of figuring out which variables we'll be using:

```r
library(MASS)

data = readRDS('~/data/baseline_prediction/prs_start/complete_massaged_data_02032020.rds')
data$externalizing = as.factor(data$externalizing)

# Q < .1, including PRS
inatt_vars = c('VMI.beery', 'DSF.wisc', 'PS.wj', 'DS.wj', 'VM.wj', 'IFO_fa',
               'ADHD_PRS0.000100', 'ADHD_PRS0.001000', 'ADHD_PRS0.000500',
               'FSIQ', 'externalizing')
hi_vars = c('VMI.beery', 'PS.wj', 'DS.wj', 'VM.wj', 'ATR_fa', 'CST_fa',
            'IFO_fa', 'ADHD_PRS0.001000', 'ADHD_PRS0.000500', 'FSIQ',
            'externalizing', 'OFC', 'cingulate')
covars = c('base_age', sapply(1:10, function(x) sprintf('PC%02d', x)),
           "meanX.trans", "meanY.trans", "meanZ.trans", "meanX.rot",
           "meanY.rot", "meanZ.rot", "goodVolumes",
           "mprage_score", "ext_avg", "int_avg", 'sex')
predictors = inatt_vars
sx = 'inatt'
min_sx = 6

if (sx == 'inatt') {
    thresh = 0
} else if (sx == 'hi') {
    thresh = -.5
}
phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data[, phen] = NA
data[which(data[, phen_slope] < thresh), phen] = 'imp'
data[which(data[, phen_slope] >= thresh), phen] = 'nonimp'
data[, phen] = factor(data[, phen], ordered=F)
data[, phen] = relevel(data[, phen], ref='imp')

predictors_str = paste(predictors, collapse='+')
fm_str = paste(phen, " ~ ", predictors_str, ' + ',
               paste(covars, collapse='+'),
               sep="")

use_me = !is.na(data[, phen])
this_data = data[use_me, c(phen, predictors, covars)]
fit_all = glm(as.formula(fm_str), data=this_data, family=binomial(link='logit'), na.action=na.exclude)

use_me = !is.na(data[, phen]) & data$bestInFamily
this_data = data[use_me, c(phen, predictors, covars)]
fit_bif = glm(as.formula(fm_str), data=this_data, family=binomial(link='logit'))

# make sure to step to remove only covariates!
step = stepAIC(fit, direction='both', trace=F,
                scope = list(lower = as.formula(sprintf('~ %s', predictors_str))))
```

Neither model is converging... maybe I'll need to impute afterall.

```r
> summary(this_data[,phen])
   imp nonimp 
    74     59 
> summary(this_data[, predictors])
   VMI.beery         DSF.wisc          PS.wj            DS.wj       
 Min.   : 49.00   Min.   : 3.000   Min.   : 41.00   Min.   : 55.00  
 1st Qu.: 83.50   1st Qu.: 8.000   1st Qu.: 87.25   1st Qu.: 91.00  
 Median : 92.00   Median :10.000   Median : 96.00   Median :100.00  
 Mean   : 92.75   Mean   : 9.821   Mean   : 96.04   Mean   : 98.98  
 3rd Qu.:100.00   3rd Qu.:12.000   3rd Qu.:106.00   3rd Qu.:109.00  
 Max.   :145.00   Max.   :16.000   Max.   :134.00   Max.   :134.00  
 NA's   :30       NA's   :55       NA's   :27       NA's   :24      
     VM.wj            IFO_fa       ADHD_PRS0.000100    ADHD_PRS0.001000   
 Min.   : 36.00   Min.   :0.3731   Min.   :-0.012671   Min.   :-0.007106  
 1st Qu.: 86.00   1st Qu.:0.3948   1st Qu.:-0.010625   1st Qu.:-0.006158  
 Median : 95.00   Median :0.4040   Median :-0.009605   Median :-0.005847  
 Mean   : 94.32   Mean   :0.4036   Mean   :-0.009610   Mean   :-0.005764  
 3rd Qu.:105.00   3rd Qu.:0.4115   3rd Qu.:-0.008492   3rd Qu.:-0.005346  
 Max.   :129.00   Max.   :0.4336   Max.   :-0.006035   Max.   :-0.004222  
 NA's   :24       NA's   :77                                              
 ADHD_PRS0.000500         FSIQ        externalizing
 Min.   :-0.008341   Min.   : 78.00   0:118        
 1st Qu.:-0.007437   1st Qu.: 94.75   1: 15        
 Median :-0.006846   Median :107.00                
 Mean   :-0.006829   Mean   :105.96                
 3rd Qu.:-0.006344   3rd Qu.:116.00                
 Max.   :-0.004740   Max.   :143.00                
                     NA's   :1                     
```

As expected, DTI is the one with the most data missing. But that looks very odd,
as aparently we have way more NAs than expected, even for neuropsych (55???).

The issue is cropping the scans to the closest to baseline date. If we use
collapsed FA for example, we start with 859 scans after qc_var cleaning. After
cleaning using the data, we get we get 804. Of those, we have 277 out of the 284
original subjects. Then, when after removing any scans longer than 1 year, we
only keep 180 subjects, which include NVs. After removing those, we only have
79!

I guess that was the original reason to look at the NVs as well... can we do a
multinomial logistic regression then?

```r
library(nnet)
library(MASS)

data = readRDS('~/data/baseline_prediction/prs_start/complete_massaged_data_02032020.rds')
data$externalizing = as.factor(data$externalizing)

# Q < .1, including PRS
inatt_vars = c('VMI.beery', 'DSF.wisc', 'PS.wj', 'DS.wj', 'VM.wj', 'IFO_fa',
               'ADHD_PRS0.000100', 'ADHD_PRS0.001000', 'ADHD_PRS0.000500',
               'FSIQ', 'externalizing')
hi_vars = c('VMI.beery', 'PS.wj', 'DS.wj', 'VM.wj', 'ATR_fa', 'CST_fa',
            'IFO_fa', 'ADHD_PRS0.001000', 'ADHD_PRS0.000500', 'FSIQ',
            'externalizing', 'OFC', 'cingulate')
covars = c('base_age', sapply(1:10, function(x) sprintf('PC%02d', x)),
           "meanX.trans", "meanY.trans", "meanZ.trans", "meanX.rot",
           "meanY.rot", "meanZ.rot", "goodVolumes",
           "mprage_score", "ext_avg", "int_avg", 'sex')
predictors = inatt_vars
sx = 'inatt'
min_sx = 6

if (sx == 'inatt') {
    thresh = 0
} else if (sx == 'hi') {
    thresh = -.5
}
phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data[, phen_slope]))
idx = data[my_nvs, 'base_inatt'] <= 2 & data[my_nvs, 'base_hi'] <= 2
data[my_nvs[idx], phen] = 'nv012'
data[which(data[, phen_slope] < thresh), phen] = 'imp'
data[which(data[, phen_slope] >= thresh), phen] = 'nonimp'
data[, phen] = factor(data[, phen], ordered=F)
data[, phen] = relevel(data[, phen], ref='nv012')

predictors_str = paste(predictors, collapse='+')
fm_str = paste(phen, " ~ ", predictors_str, ' + ',
               paste(covars, collapse='+'),
               sep="")

use_me = !is.na(data[, phen])
this_data = data[use_me, c(phen, predictors, covars)]
scale_me = c()
for (v in )
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
fit_all = multinom(as.formula(fm_str), data=this_data)

use_me = !is.na(data[, phen]) & data$bestInFamily
this_data = data[use_me, c(phen, predictors, covars)]
fit_bif = glm(as.formula(fm_str), data=this_data, family=binomial(link='logit'))

# # make sure to step to remove only covariates!
# step = stepAIC(fit_all, direction='both', trace=F,
#                 scope = list(lower = as.formula(sprintf('~ %s', predictors_str))))

z <- summary(fit_all)$coefficients/summary(fit_all)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
```

But even by doing something like that I'm running into trouble in the ps...

After chatting with Philip, a few more things we can try. First, does ordered
logistic regression help at all?

```r
library(nnet)
library(MASS)

data = readRDS('~/data/baseline_prediction/prs_start/complete_massaged_data_02032020.rds')
data$externalizing = as.factor(data$externalizing)

# Q < .1, including PRS but after cleaning correlated variables
hi_vars = c('VMI.beery', 'VM.wj', 'FSIQ', 'externalizing', 'IFO_fa', 'DS.wj',
            'ADHD_PRS0.001000', 'OFC', 'ATR_fa', 'CST_fa', 'cingulate',
            'DSF.wisc')
inatt_vars = c('FSIQ', 'VMI.beery', 'VM.wj', 'externalizing',
               'ADHD_PRS0.000500', 'DSF.wisc', 'IFO_fa', 'DS.wj')
covars = c('base_age', sapply(1:10, function(x) sprintf('PC%02d', x)),
           "meanX.trans", "meanY.trans", "meanZ.trans", "meanX.rot",
           "meanY.rot", "meanZ.rot", "goodVolumes",
           "mprage_score", "ext_avg", "int_avg", 'sex')
predictors = inatt_vars
sx = 'inatt'
min_sx = 6

if (sx == 'inatt') {
    thresh = 0
} else if (sx == 'hi') {
    thresh = -.5
}
phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data[, phen_slope]))
idx = data[my_nvs, 'base_inatt'] <= 2 & data[my_nvs, 'base_hi'] <= 2
data[my_nvs[idx], phen] = 'nv012'
data[which(data[, phen_slope] < thresh), phen] = 'imp'
data[which(data[, phen_slope] >= thresh), phen] = 'nonimp'
data[, phen] = factor(data[, phen], ordered=F)
data[, phen] = relevel(data[, phen], ref='nv012')

predictors_str = paste(predictors, collapse='+')
fm_str = paste(phen, " ~ ", predictors_str, ' + ',
               paste(covars, collapse='+'),
               sep="")

use_me = !is.na(data[, phen])
this_data = data[use_me, c(phen, predictors, covars)]
scale_me = c()
for (v in c(predictors, covars)) {
    if (!is.factor(this_data[, v])) {
        scale_me = c(scale_me, v)
    }
}
this_data[, scale_me] = scale(this_data[, scale_me])
fit_all = multinom(as.formula(fm_str), data=this_data, maxit=2000)

use_me = !is.na(data[, phen]) & data$bestInFamily
this_data = data[use_me, c(phen, predictors, covars)]
scale_me = c()
for (v in c(predictors, covars)) {
    if (!is.factor(this_data[, v])) {
        scale_me = c(scale_me, v)
    }
}
this_data[, scale_me] = scale(this_data[, scale_me])
fit_bif = multinom(as.formula(fm_str), data=this_data, maxit=2000)

# make sure to step to remove only covariates!
step = stepAIC(fit_all, direction='both', trace=F,
                scope = list(lower = as.formula(sprintf('~ %s', predictors_str))))

z <- summary(fit_all)$coefficients/summary(fit_all)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2
```

Now I need to assess variable importance here. Let's see if dominance analysis
works with multinomial logistic regression. If not, I could find a different way
to quantify variable importance, or run several pairwise logistic regressions,
having pie-charts for each comparison. That might actually be a bit more
informative, but I don't think it'll converge... no, it didn't. 

OK, so next steps are to try dominance analysis in multinomial, or find other
ways to quantify the variables.





# TODO
* impute data?
* does it work with lme?
* does it work with multinomial?
* add base_sx to see what happens to the variable importance? 