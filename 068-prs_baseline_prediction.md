# 2020-01-24 13:02:29

Let's continue what we were doing in 067, but I'm not going to pick up from
persistence PRS. After talking to Philip, it's a better approach to see where
NVs lie on the PRS spectrum. In other ways, can we choose the best outcome
formulation by including NVs in the mix. Where are they with respect to
improvers and non-improvers?

```r
data_prs = read.csv('/Volumes/Shaw/tmp/for_philip/working_gf.csv')
phen = 'bin0.33_hi_GE6_wp05'
data_prs[, phen] = as.character(data_prs[, phen])
imnv = is.na(data_prs[, phen])
data_prs[imnv, phen] = 'nv'
use_me = data_prs$bestInFamily
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
library("nnet")
fit = multinom(bin0.33_hi_GE6_wp05~ADHD_PRS0.000100+ PC01 + PC02 + PC03 + PC04 + PC05 + PC06 + PC07 + PC08 + PC09 + PC10 + base_age + sex.x, data=data_prs[use_me,])
z <- summary(fit)$coefficients/summary(fit)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1))*2
```

The code above seems to work. So, let's loop through our possible thresholds:

```r
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = '+ PC01 + PC02 + PC03 + PC04 + PC05 + PC06 + PC07 + PC08 + PC09 + PC10 + base_age + sex.x'
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_all_PCsAgeSex_multinom.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(3, 4, 6)) {
        for (qtile in c(.2, .25, .33, .5)) {
            phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
            phen = sprintf('bin%.2f_%s_GE%d_wp05', qtile, sx, min_sx)
            thresh = quantile(data_prs[, phen_slope], qtile, na.rm=T)
            data_prs[, phen] = 'nv'
            data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
            data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
            data_prs[, phen] = factor(data_prs[, phen], ordered=F)
            data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
            use_me = data_prs$bestInFamily #& data_prs$isWNH

            phen_res = c()
            for (prs in prs_var_names) {
                fm_str = paste(phen, "~", prs, covars, sep="")
                fit = multinom(as.formula(fm_str), data=data_prs[use_me,])
                z <- summary(fit)$coefficients/summary(fit)$standard.errors
                p <- (1 - pnorm(abs(z), 0, 1))*2
                temp = c(p[, 2], summary(fit)$AIC, summary(fit)$deviance)
                phen_res = rbind(phen_res, temp)
                rownames(phen_res)[nrow(phen_res)] = fm_str
            }
            phen_res = data.frame(phen_res)
            phen_res$formula = rownames(phen_res)
            phen_res$predictor = prs_var_names
            phen_res$outcome = phen
            hold = rbind(hold, phen_res)
        }
    }
}
colnames(hold)[3:4] = c('AIC', 'deviance')
write.csv(hold, file=out_fname, row.names=F)
```

This code is working, but way too many results have p-values too good to be
true. It's not an issue with the code either... it's something funky going on.
But what if we invert the model? Isn't the idea just to show that PRS is
different for the 3 different groups? Then, why not just do an ANOVA?

```r
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
covars = ''#'PC01 + PC02 + PC03 + PC04 + PC05 + PC06 + PC07 + PC08 + PC09 + PC10 + base_age + sex.x + '
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_WNH_noCovs_aov.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(3, 4, 6)) {
        for (qtile in c(.2, .25, .33, .5)) {
            phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
            phen = sprintf('bin%.2f_%s_GE%d_wp05', qtile, sx, min_sx)
            thresh = quantile(data_prs[, phen_slope], qtile, na.rm=T)
            data_prs[, phen] = 'nv'
            data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
            data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
            data_prs[, phen] = factor(data_prs[, phen], ordered=F)
            data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
            use_me = data_prs$bestInFamily & data_prs$isWNH

            phen_res = c()
            for (prs in prs_var_names) {
                fm_str = paste(prs, "~", covars, phen, sep="")
                fit = aov(as.formula(fm_str), data=data_prs[use_me,])
                p = summary(fit)[[1]]
                temp = c(p[nrow(p)-1, 5], p[nrow(p), 3])
                phen_res = rbind(phen_res, temp)
                rownames(phen_res)[nrow(phen_res)] = fm_str
            }
            phen_res = data.frame(phen_res)
            phen_res$formula = rownames(phen_res)
            phen_res$predictor = prs_var_names
            phen_res$outcome = phen
            hold = rbind(hold, phen_res)
        }
    }
}
colnames(hold)[1:2] = c('pval', 'residMeanSq')
write.csv(hold, file=out_fname, row.names=F)
```

Accumulated results in univariate_aov_results.xlsx. If we want to include inatt
and hi results, it seems like the best thing is to use bin.5 and GE6. That's
fine as well. We just need to make sure our ANCOVA assumptions hold in that
threshold, and that the posthoc comparisons look fine (this is a great tutorial,
http://faculty.missouri.edu/huangf/data/quantf/ancova_in_r_handout.pdf)

```r
sx = 'inatt'
min_sx = 6
qtile = .5
prs = 'ADHD_PRS0.000100'
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age', 'sex.x')
phen = sprintf('bin%.2f_%s_GE%d_wp05', qtile, sx, min_sx)
thresh = quantile(data_prs[, phen_slope], qtile, na.rm=T)
data_prs[, phen] = 'nv'
data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
use_me = data_prs$bestInFamily

# interaction between groups and covariates
for (covar in covars) {
    fm_str = sprintf('%s ~ %s + %s + %s:%s', prs, phen, covar, phen, covar)
    fit = aov(as.formula(fm_str), data=data_prs[use_me,])
    ps = summary(fit)[[1]]
    pval = ps[nrow(ps)-1, 5]
    if (pval < .05) {
        print(sprintf('pval for %s is significant', covar))
    }
}
```

PC04 violates it for ADHD_PRS0.001000, 



# TODO
* try persistence PRS
* continue work on glmer model just for robustness 