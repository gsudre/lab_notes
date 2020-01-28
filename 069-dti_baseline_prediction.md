# 2020-01-28 10:50:26

I'll follow the same idea as before (068). But now I have an idea of what works
best for PRS in tems of thresholding, so I'll keep those. But the approach is to
start with an ANOVA mixed model, and then calculate where the differences are.
Then, see where the linear model works best.

First we'll need to grab the DTI data:

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

b = read.csv('/Volumes/Shaw/MasterQC/master_qc_20190314.csv')
a = read.csv('~/data/heritability_change/ready_1020.csv')
# m has all 1020 processed scans
m = merge(a, b, by.y='Mask.ID', by.x='Mask.ID...Scan', all.x=F)

# dti_meta has all scans for everyone who has PRS data (954 scans)
dti_meta = merge(m, df, by.x="Medical.Record...MRN...Subjects", by.y='MRN',
          all.x=F, all.y=F)

# restrict based on QC
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
dti_meta = dti_meta[dti_meta$"age_at_scan...Scan...Subjects" < 18, ]
dti_meta = dti_meta[dti_meta$"goodVolumes" <= 61, ]
dti_meta = dti_meta[dti_meta$"numVolumes" < 80, ]

# down to 928 scans that obey criteria and have PRS
library(solitude)
iso <- isolationForest$new()
iso$fit(dti_meta[, qc_vars])
scores_if = as.matrix(iso$scores)[,3]
library(dbscan)
# here I set the number of neighbors to a percentage of the total data
scores_lof = lof(dti_meta[, qc_vars], k = round(.5 * nrow(dti_meta)))

qtile=.95
thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

tracts = read.csv('~/data/heritability_change/jhu_tracts_1020.csv')
# somehow I have two entries for 1418?
x = duplicated(tracts$id)
jhu_data = merge(dti_meta[idx,], tracts[!x, ], by.x='Mask.ID...Scan', by.y='id')
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])

iso <- isolationForest$new()
iso$fit(jhu_data[, tract_names])
scores_if = as.matrix(iso$scores)[,3]
scores_lof = lof(jhu_data[, tract_names], k = round(.5 * nrow(jhu_data)))

thresh_lof = quantile(scores_lof, qtile)
thresh_if = quantile(scores_if, qtile)
idx = scores_lof < thresh_lof & scores_if < thresh_if

clean_jhu_data = jhu_data[idx, ]

# down to 800 scans when only scans at .95 in both criteria are used

# selecting earliest scan for each subject, regardless of score, as we're assu ing every scan now is good
keep_me = c()
for (s in unique(clean_jhu_data$Medical.Record...MRN...Subjects)) {
    subj_rows = which(clean_jhu_data$Medical.Record...MRN...Subjects == s)
    subj_data = clean_jhu_data[subj_rows, ]
    min_subj_row = which.min(subj_data$age_at_scan...Scan...Subjects)
    keep_me = c(keep_me, subj_rows[min_subj_row])
}
data_dti = clean_jhu_data[keep_me, ]
# finished with 277 scans when using baseline for each subject
demo = read.csv('prs_demo.csv')
# just to get FAMID
data_dti = merge(data_dti, demo, by.x='Medical.Record...MRN...Subjects',
                 by.y='MRN')
```

Before we add FA or MO, which will need some extra work and will also change the
scans we're using, let's run our usual models, with the thresholds we got from
PRS:

```r
library(nlme)
library(car)

hold = c()
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
covars = c(qc_vars, 'age_at_scan...Scan...Subjects')
out_fname = '~/data/baseline_prediction/prs_start/univar_JHUtractsADRD_all_PCsAgeSex_lme.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    } else {
        thresholds = -1
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data_dti[, phen] = 'notGE6adhd'
    my_nvs = which(is.na(data_dti[, phen_slope]))
    idx = data_dti[my_nvs, 'base_inatt'] <= 2 & data_dti[my_nvs, 'base_hi'] <= 2
    data_dti[my_nvs[idx], phen] = 'nv012'
    data_dti[which(data_dti[, phen_slope] < thresh), phen] = 'imp'
    data_dti[which(data_dti[, phen_slope] >= thresh), phen] = 'nonimp'
    data_dti[, phen] = factor(data_dti[, phen], ordered=F)
    data_dti[, phen] = relevel(data_dti[, phen], ref='nv012')
    use_me = T#data_dti$isWNH

    this_data = data_dti[use_me, c(phen, 'FAMID', tract_names, covars)]
    this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
    this_data$sex = data_dti[use_me, 'sex.x']
    tmp_covars = c(covars, 'sex')
    phen_res = c()
    for (tract in tract_names) {
        fm_str = paste(tract, " ~ ", phen, " +",
                           paste(tmp_covars, collapse='+'),
                           sep="")
        fit = lme(as.formula(fm_str), ~1|FAMID, data=this_data)
        p = Anova(fit)
        temp = c(p[1,3], summary(fit)$logLik, summary(fit)$AIC,
                 summary(fit)$BIC)
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
    }
    phen_res = data.frame(phen_res)
    phen_res$formula = rownames(phen_res)
    phen_res$predictor = tract_names
    phen_res$outcome = phen
    hold = rbind(hold, phen_res)
}
colnames(hold)[1:4] = c('pval', 'logLik', 'AIC', 'BIC')
write.csv(hold, file=out_fname, row.names=F)
```

![](images/2020-01-28-11-26-48.png)

We do have some results, mostly hi. Maybe it'd survive some sort of comparison
correction... let's plot the results to see where the group differences lie:

```r
library(multcomp)
library(ggplot2)

sx = 'hi'
min_sx = 6
thresh = -.5
use_me = T
my_tracts = c('ad_10', 'ad_7', 'rd_18', 'ad_16', 'rd_10',
              'rd_3', 'rd_16', 'ad_8', 'ad_6', 'rd_17', 'rd_11')
# my_tracts = c('ad_12', 'ad_10', 'ad_6')

phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data_dti[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data_dti[, phen_slope]))
idx = data_dti[my_nvs, 'base_inatt'] <= 2 & data_dti[my_nvs, 'base_hi'] <= 2
data_dti[my_nvs[idx], phen] = 'nv012'
data_dti[which(data_dti[, phen_slope] < thresh), phen] = 'imp'
data_dti[which(data_dti[, phen_slope] >= thresh), phen] = 'nonimp'
data_dti[, phen] = factor(data_dti[, phen], ordered=F)
data_dti[, phen] = relevel(data_dti[, phen], ref='nv012')
use_me = T

this_data = data_dti[use_me, c(phen, 'FAMID', tract_names, covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$sex = data_dti[use_me, 'sex.x']
tmp_covars = c(covars, 'sex')
for (tract in my_tracts) {
    fm_str = paste(tract, " ~ ", phen, " +",
                        paste(tmp_covars, collapse='+'),
                        sep="")
    fit = lme(as.formula(fm_str), ~1|FAMID, data=this_data)
    p = Anova(fit)
    posthoc = glht(fit, linfct=mcp(thresh0.00_inatt_GE6_wp05 = "Tukey"))
    # posthoc = glht(fit, linfct=mcp(thresh0.50_hi_GE6_wp05 = "Tukey"))
    sig_idx = summary(posthoc)$test$pvalues < .05
    print(sprintf('%s, %s', phen, tract))
    print(names(coef(posthoc))[sig_idx])
}
```

This is what I get:

```
[1] "thresh0.50_hi_GE6_wp05, ad_10"
[1] "imp - nv012"
[1] "thresh0.50_hi_GE6_wp05, ad_7"
[1] "imp - nv012"  "nonimp - imp"
[1] "thresh0.50_hi_GE6_wp05, rd_18"
[1] "nonimp - nv012"      "nonimp - imp"        "notGE6adhd - nonimp"
[1] "thresh0.50_hi_GE6_wp05, ad_16"
character(0)
[1] "thresh0.50_hi_GE6_wp05, rd_10"
[1] "imp - nv012"  "nonimp - imp"
[1] "thresh0.50_hi_GE6_wp05, rd_3"
[1] "nonimp - nv012"
[1] "thresh0.50_hi_GE6_wp05, rd_16"
character(0)
[1] "thresh0.50_hi_GE6_wp05, ad_8"
[1] "imp - nv012"
[1] "thresh0.50_hi_GE6_wp05, ad_6"
[1] "notGE6adhd - nv012"
[1] "thresh0.50_hi_GE6_wp05, rd_17"
[1] "nonimp - imp"
[1] "thresh0.50_hi_GE6_wp05, rd_11"
[1] "nonimp - nv012"

[1] "thresh0.00_inatt_GE6_wp05, ad_12"
[1] "notGE6adhd - imp"
[1] "thresh0.00_inatt_GE6_wp05, ad_10"
[1] "nonimp - nv012"
[1] "thresh0.00_inatt_GE6_wp05, ad_6"
[1] "notGE6adhd - nv012"
```

Very few have the direction we'd expect. Maybe it's just easier to test for
linearity, like before:

```r
hold = c()
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
covars = c(qc_vars, 'age_at_scan...Scan...Subjects')
out_fname = '~/data/baseline_prediction/prs_start/univar_JHUtractsADRD_all_PCsAgeSex_4groupOrdered_lme.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    } else {
        thresholds = -1
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data_dti[, phen] = 'notGE6adhd'
    my_nvs = which(is.na(data_dti[, phen_slope]))
    idx = data_dti[my_nvs, 'base_inatt'] <= 2 & data_dti[my_nvs, 'base_hi'] <= 2
    data_dti[my_nvs[idx], phen] = 'nv012'
    data_dti[which(data_dti[, phen_slope] < thresh), phen] = 'imp'
    data_dti[which(data_dti[, phen_slope] >= thresh), phen] = 'nonimp'
    data_dti[, phen] = factor(data_dti[, phen], ordered=F)
    data_dti[, phen] = relevel(data_dti[, phen], ref='nv012')
    use_me = T#data_dti$isWNH

    this_data = data_dti[use_me, c(phen, 'FAMID', tract_names, covars)]
    this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
    this_data$sex = data_dti[use_me, 'sex.x']
    tmp_covars = c(covars, 'sex')
    this_data$ordered = factor(this_data[, phen],
                           levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                           ordered=T)
    phen_res = c()
    for (tract in tract_names) {
        fm_str = paste(tract, " ~ ordered +",
                           paste(tmp_covars, collapse='+'),
                           sep="")
        fit = lme(as.formula(fm_str), ~1|FAMID, data=this_data)
        temp = c(summary(fit)$tTable['ordered.L', ],
                     summary(fit)$logLik, summary(fit)$AIC, summary(fit)$BIC)
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
    }
    phen_res = data.frame(phen_res)
    phen_res$formula = rownames(phen_res)
    phen_res$predictor = tract_names
    phen_res$outcome = phen
    hold = rbind(hold, phen_res)
}
colnames(hold)[6:8] = c('logLik', 'AIC', 'BIC')
write.csv(hold, file=out_fname, row.names=F)
```

![](images/2020-01-28-12-41-12.png)

We get the usual suspects, such as 18 and 19, which are uncinate-R and SLF-L,
respectively. For reference, _3 is CST-L and _11 is IFO-L.

Does it work the other way around as well, as we predict those classes using the
tracts?

```r
library(ordinal)
hold = c()
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
covars = c(qc_vars, 'age_at_scan...Scan...Subjects')
out_fname = '~/data/baseline_prediction/prs_start/univar_JHUtractsADRD_all_PCsAgeSex_4groupOrdered_clmm2.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    } else {
        thresholds = -1
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data_dti[, phen] = 'notGE6adhd'
    my_nvs = which(is.na(data_dti[, phen_slope]))
    idx = data_dti[my_nvs, 'base_inatt'] <= 2 & data_dti[my_nvs, 'base_hi'] <= 2
    data_dti[my_nvs[idx], phen] = 'nv012'
    data_dti[which(data_dti[, phen_slope] < thresh), phen] = 'imp'
    data_dti[which(data_dti[, phen_slope] >= thresh), phen] = 'nonimp'
    data_dti[, phen] = factor(data_dti[, phen], ordered=F)
    data_dti[, phen] = relevel(data_dti[, phen], ref='nv012')
    use_me = T#data_dti$isWNH

    this_data = data_dti[use_me, c(phen, 'FAMID', tract_names, covars)]
    this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
    this_data$sex = data_dti[use_me, 'sex.x']
    tmp_covars = c(covars, 'sex')
    this_data$ordered = factor(this_data[, phen],
                           levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                           ordered=T)
    this_data$FAMID = factor(this_data$FAMID)
    phen_res = c()
    for (tract in tract_names) {
        fm_str = paste("ordered ~", tract, '+',
                        paste(tmp_covars, collapse='+'),
                        sep="")
        fit = clmm2(as.formula(fm_str), random=FAMID, data=this_data, Hess=T)
        temp = c(summary(fit)$coefficients[tract, 'Pr(>|z|)'],
                    summary(fit)$logLik, summary(fit)$condHess)
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
    }
    phen_res = data.frame(phen_res)
    phen_res$formula = rownames(phen_res)
    phen_res$predictor = tract_names
    phen_res$outcome = phen
    hold = rbind(hold, phen_res)
}
colnames(hold)[1:3] = c('pval', 'loglik', 'hessian')
write.csv(hold, file=out_fname, row.names=F)
```

![](images/2020-01-28-13-41-41.png)

Yes, we get some similar hits, which is nice. Now, I wonder if I'd get better
results if I did just the multinomial logistic regression, without fixing the
group order?

My difficulty here is that we could potentially have linear, quadratic or even
cubic relationships, depending on the brain region. We could come up with
stories for any of the relationships. Let's see if some regions are captured
better than others:

```r
hold = c()
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])
qc_vars = c("meanX.trans", "meanY.trans", "meanZ.trans",
            "meanX.rot", "meanY.rot", "meanZ.rot",
            "goodVolumes")
covars = c(qc_vars, 'age_at_scan...Scan...Subjects')
out_fname = '~/data/baseline_prediction/prs_start/univar_JHUtractsADRD_all_PCsAgeSex_4groupOrdered_lme_allFits.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    } else {
        thresholds = -1
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data_dti[, phen] = 'notGE6adhd'
    my_nvs = which(is.na(data_dti[, phen_slope]))
    idx = data_dti[my_nvs, 'base_inatt'] <= 2 & data_dti[my_nvs, 'base_hi'] <= 2
    data_dti[my_nvs[idx], phen] = 'nv012'
    data_dti[which(data_dti[, phen_slope] < thresh), phen] = 'imp'
    data_dti[which(data_dti[, phen_slope] >= thresh), phen] = 'nonimp'
    data_dti[, phen] = factor(data_dti[, phen], ordered=F)
    data_dti[, phen] = relevel(data_dti[, phen], ref='nv012')
    use_me = T#data_dti$isWNH

    this_data = data_dti[use_me, c(phen, 'FAMID', tract_names, covars)]
    this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
    this_data$sex = data_dti[use_me, 'sex.x']
    tmp_covars = c(covars, 'sex')
    this_data$ordered = factor(this_data[, phen],
                           levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                           ordered=T)
    phen_res = c()
    for (tract in tract_names) {
        fm_str = paste(tract, " ~ ordered +",
                           paste(tmp_covars, collapse='+'),
                           sep="")
        fit = lme(as.formula(fm_str), ~1|FAMID, data=this_data)
        temp = c(summary(fit)$tTable['ordered.L', ],
                     summary(fit)$logLik, summary(fit)$AIC, summary(fit)$BIC,
                     tract, 'linear')
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
        temp = c(summary(fit)$tTable['ordered.Q', ],
                     summary(fit)$logLik, summary(fit)$AIC, summary(fit)$BIC,
                     tract, 'quadratic')
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
        temp = c(summary(fit)$tTable['ordered.C', ],
                     summary(fit)$logLik, summary(fit)$AIC, summary(fit)$BIC,
                     tract, 'cubic')
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
    }
    phen_res = data.frame(phen_res)
    phen_res$formula = rownames(phen_res)
    phen_res$outcome = phen
    hold = rbind(hold, phen_res)
}
colnames(hold)[6:10] = c('logLik', 'AIC', 'BIC', 'tract', 'modtype')
write.csv(hold, file=out_fname, row.names=F)
```

We actually get more fits for quadratic and cubic... we will eventually get into
an issue of multiple comparisons here. 

Here's another idea. Philip's daisywheel had an inspiring title:

![](images/2020-01-28-15-07-11.png)

So, why not first filter our variables on how well they predict outcome of ADHD
(simple, not multinomial) logistic regression, and then run other models to
include the NVs, just to see where the NVs are placed in those interesting
regions?

If this approach works, I'll reproduce it for PRS...

```r
hold = c()
tract_names = c(colnames(tracts)[grepl(colnames(tracts), pattern="^ad")],
                colnames(tracts)[grepl(colnames(tracts), pattern="^rd")])
qc_vars = c("norm.trans", "norm.rot", "goodVolumes")
covars = c(qc_vars, 'age_at_scan...Scan...Subjects')
out_fname = '~/data/baseline_prediction/prs_start/univar_JHUtractsADRD_all_GE6_norm_outcomeOnly.csv'
for (sx in c('inatt', 'hi')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresh = 0
    } else if (sx == 'hi') {
        thresh = -.5
    } else {
        thresh = -.48
    }
    phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
    phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
    data_dti[, phen] = NA
    data_dti[which(data_dti[, phen_slope] < thresh), phen] = 'imp'
    data_dti[which(data_dti[, phen_slope] >= thresh), phen] = 'nonimp'
    data_dti[, phen] = factor(data_dti[, phen], ordered=F)
    data_dti[, phen] = relevel(data_dti[, phen], ref='imp')
    use_me = !is.na(data_dti[, phen]) #data_dti$isWNH

    this_data = data_dti[use_me, c(phen, 'FAMID', tract_names, covars)]
    this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
    this_data$sex = data_dti[use_me, 'sex.x']
    tmp_covars = c(covars, 'sex')
    phen_res = c()
    for (tract in tract_names) {
        fm_str = paste(phen, " ~ ", tract, ' + ',
                           paste(tmp_covars, collapse='+'), ' + (1|FAMID)',
                           sep="")
        fit = glmer(as.formula(fm_str), data=this_data,
                    family=binomial(link='logit'))
        if (isSingular(fit)) {
            temp = c(summary(fit)$coefficients[tract, ], summary(fit)$AIC[1:3])
        } else {
            temp = rep(NA, 7)
        }
        phen_res = rbind(phen_res, temp)
        rownames(phen_res)[nrow(phen_res)] = fm_str
    }
    phen_res = data.frame(phen_res)
    phen_res$formula = rownames(phen_res)
    phen_res$predictor = tract_names
    phen_res$outcome = phen
    colnames(phen_res) =  c("Estimate", "Std..Error", "z.value", "Pr...z..",
                            "AIC", "BIC", "logLik", "formula", "predictor",
                            "outcome")
    hold = rbind(hold, phen_res)
}
write.csv(hold, file=out_fname, row.names=F)
```

The models are not converging... but that's because we're leaving a whole bunch
of data on the floor by removing NVs. Using the norms for covariates didn't help
much either. I'm going to try to just use the best kid in each family to make
the models converge. If they do, we could potentially use the left out kids for
validation? Only 19 out of the 133 scans we're using (ADHD only) will be thrown
away... maybe this will work.

<!-- data_prs$bestInFamily = F
nvisits = table(clin_long$MRN)
data_prs = merge(data_prs, as.matrix(nvisits), by.x='MRN', by.y=0)
colnames(data_prs)[ncol(data_prs)] = 'nvisits'
for (f in unique(data_prs$FAMID)) {
    fam_rows = which(data_prs$FAMID == f)
    fam_data = data_prs[fam_rows,]
    if (nrow(fam_data) == 1) {
        data_prs[fam_rows,]$bestInFamily = T
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
                data_prs[fam_rows[svisits$ix[1]], ]$bestInFamily = T
            }
        } else {
            data_prs[fam_rows[stotal$ix[1]], ]$bestInFamily = T
        }
    }
}
```

There are only 6 ties, so I can select them manually.

```r
data_prs[data_prs$MRN==4585574, ]$bestInFamily = T
data_prs[data_prs$MRN==4925051, ]$bestInFamily = T
data_prs[data_prs$MRN==7079035, ]$bestInFamily = T
data_prs[data_prs$MRN==7378993, ]$bestInFamily = T
# chosen because of overall best MPRAGE QC
data_prs[data_prs$MRN==4640378, ]$bestInFamily = T
# chosen because of overall best MPRAGE QC
data_prs[data_prs$MRN==7218965, ]$bestInFamily = T
``` -->



# TODO
* check what are the group differences
* add mode and FA
* could we do probability function easier if our logistic regression weren't
  ordered? 
# links

https://strengejacke.github.io/ggeffects/articles/logisticmixedmodel.html
https://www.researchgate.net/post/How_to_get_P-value_associated_to_explanatory_from_binomial_glmer