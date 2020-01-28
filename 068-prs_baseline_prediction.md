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

# 2020-01-27 09:01:48

I was chatting with Philip and there are a few things to try:

* Focus on mixed model ANOVA so we can include everyone. Then we can do it just
one per family to confirm the results
* Use more meaningful breaks in the data, such as 0 or -.5.
* Maybe stratify NVs into squeaky clean (0 and 1), then 2, 4, 5, and then GE6.

Let's run the model as is first.

```r
library(lme4)
library(car)

hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
add_sex = T
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_all_PCsAgeSex_lmer.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(3, 4, 6)) {
        for (qtile in c(.33, .5)) {
            phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
            phen = sprintf('bin%.2f_%s_GE%d_wp05', qtile, sx, min_sx)
            thresh = quantile(data_prs[, phen_slope], qtile, na.rm=T)
            data_prs[, phen] = 'nv'
            data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
            data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
            data_prs[, phen] = factor(data_prs[, phen], ordered=F)
            data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
            use_me = T#data_prs$isWNH
            print(sprintf('%s GE%d at %.2f = %.2f', sx, min_sx, qtile, thresh))

            this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                                           covars)]
            this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
            if (add_sex) {
                this_data$sex = data_prs[use_me, 'sex.x']
                tmp_covars = c(covars, 'sex')
            } else {
                tmp_covars = covars
            }
            phen_res = c()
            for (prs in prs_var_names) {
                fm_str = paste(prs, "~", phen, '+',
                               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
                               sep="")
                fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
                p = Anova(fit)
                temp = c(p[1,3], summary(fit)$AIC[1:4])
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
colnames(hold)[1] = c('pval')
write.csv(hold, file=out_fname, row.names=F)
```

The thresholds don't depend on all or WNH. I get:

```
[1] "inatt GE3 at 0.33 = -0.33"
[1] "inatt GE3 at 0.50 = -0.08"
[1] "inatt GE4 at 0.33 = -0.33"
[1] "inatt GE4 at 0.50 = -0.08"
[1] "inatt GE6 at 0.33 = -0.33"
[1] "inatt GE6 at 0.50 = -0.10"
[1] "hi GE3 at 0.33 = -0.57"
[1] "hi GE3 at 0.50 = -0.36"
[1] "hi GE4 at 0.33 = -0.57"
[1] "hi GE4 at 0.50 = -0.36"
[1] "hi GE6 at 0.33 = -0.57"
[1] "hi GE6 at 0.50 = -0.38"
[1] "total GE3 at 0.33 = -0.73"
[1] "total GE3 at 0.50 = -0.44"
[1] "total GE4 at 0.33 = -0.74"
[1] "total GE4 at 0.50 = -0.43"
[1] "total GE6 at 0.33 = -0.76"
[1] "total GE6 at 0.50 = -0.48"
```

So it seems like natural thresholds for inatt would be 0 and -.3. For hi I have
-.5 and -.3. For total I have -.4 and -.8. I'll round them -.5 and -1. Let's
try those thresholds then:

```r
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
add_sex = T
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_WNH_PCsAgeSex_lmer.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(3, 4, 6)) {
        if (sx == 'inatt') {
            thresholds = c(0, -.3)
        } else if (sx == 'hi') {
            thresholds = c(-.3, -.5)
        } else {
            thresholds = c(-.5, -1)
        }
        for (thresh in thresholds) {
            phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
            phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
            data_prs[, phen] = 'nv'
            data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
            data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
            data_prs[, phen] = factor(data_prs[, phen], ordered=F)
            data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
            use_me = data_prs$isWNH

            this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                                           covars)]
            this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
            if (add_sex) {
                this_data$sex = data_prs[use_me, 'sex.x']
                tmp_covars = c(covars, 'sex')
            } else {
                tmp_covars = covars
            }
            phen_res = c()
            for (prs in prs_var_names) {
                fm_str = paste(prs, "~", phen, '+',
                               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
                               sep="")
                fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
                p = Anova(fit)
                temp = c(p[1,3], summary(fit)$AIC[1:4])
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
colnames(hold)[1] = c('pval')
write.csv(hold, file=out_fname, row.names=F)
```

And let's check on how it works if we split the NV groups as well:

```r
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
add_sex = T
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_all_PCsAgeSex_4_groups_lmer.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresholds = c(0, -.3)
    } else if (sx == 'hi') {
        thresholds = c(-.3, -.5)
    } else {
        thresholds = c(-.5, -1)
    }
    for (thresh in thresholds) {
        phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
        data_prs[, phen] = 'notGE6adhd'
        my_nvs = which(is.na(data_prs[, phen_slope]))
        idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
        data_prs[my_nvs[idx], phen] = 'nv012'
        data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
        data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
        data_prs[, phen] = factor(data_prs[, phen], ordered=F)
        data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
        use_me = T#data_prs$isWNH

        this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                                       covars)]
        this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
        if (add_sex) {
            this_data$sex = data_prs[use_me, 'sex.x']
            tmp_covars = c(covars, 'sex')
        } else {
            tmp_covars = covars
        }
        phen_res = c()
        for (prs in prs_var_names) {
            fm_str = paste(prs, "~", phen, '+',
                           paste(tmp_covars, collapse='+'), '+(1|FAMID)',
                           sep="")
            fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
            p = Anova(fit)
            temp = c(p[1,3], summary(fit)$AIC[1:4])
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
colnames(hold)[1] = c('pval')
write.csv(hold, file=out_fname, row.names=F)
```

I compiled the results under univariate_lmer_results.xlsx. There were a few hits
that worked for the 4group cases, and also for all and WNH. That being the
criteria, I'l keep HI at -.5 and inatt at 0. Let's make a few plots to see what
they mean:

```r
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
sx = 'hi'
prs = 'ADHD_PRS0.000500'
min_sx = 6
thresh = -.5
use_me = T#data_prs$isWNH

phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data_prs[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data_prs[, phen_slope]))
idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
data_prs[my_nvs[idx], phen] = 'nv012'
data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                               covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$sex = data_prs[use_me, 'sex.x']
tmp_covars = c(covars, 'sex')
fm_str = paste(prs, "~", phen, '+',
               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
               sep="")
fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
print(Anova(fit))

# posthoc = glht(fit, linfct=mcp(thresh0.00_inatt_GE6_wp05 = "Tukey"))
posthoc = glht(fit, linfct=mcp(thresh0.50_hi_GE6_wp05 = "Tukey"))
print(summary(posthoc))
tmp <- as.data.frame(confint(posthoc)$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
    geom_errorbar() + geom_point() + ggtitle(sprintf('%s, %s', phen, prs)) +
    theme(plot.title = element_text(hjust = 0.5))
```

![](images/2020-01-27-11-40-45.png)

So, this is a good result in that the meaningful differences are between nonimp
and imp, and nonimp and clean nvs. But we do have many comparisons here... I
could just keep nvs as nv012, or try the 2 group result. Let me see what some of
the other 4 group results look like:

![](images/2020-01-27-11-45-28.png)

![](images/2020-01-27-11-49-05.png)

![](images/2020-01-27-15-08-18.png)

![](images/2020-01-27-15-09-39.png)

So, it looks like most of the differences are in the nonimp to clean nv group.
That's make sense, especially since our PRS comes from ADHD GWAS. The first
result had a difference between nonimp and imp, which looked interestingly. It's
up to us whether we want this comparison to be more evident though. What if we
focus on WNH, or just the 3 gorup differences?

Philip also asked me to explore ordered differences. Let me see if it makes any
differences in the current result:

```r
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
sx = 'hi'
prs = 'ADHD_PRS0.000500'
min_sx = 6
thresh = -0.5
use_me = T#data_prs$isWNH

phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data_prs[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data_prs[, phen_slope]))
idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
data_prs[my_nvs[idx], phen] = 'nv012'
data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                               covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$sex = data_prs[use_me, 'sex.x']
tmp_covars = c(covars, 'sex')
fm_str = paste(prs, "~", phen, '+',
               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
               sep="")
fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
print(Anova(fit))
# posthoc = glht(fit, linfct=mcp(thresh0.00_inatt_GE6_wp05 = "Tukey"))
posthoc = glht(fit, linfct=mcp(thresh0.50_inatt_GE6_wp05 = "Tukey"))
print(summary(posthoc))

data_prs$ordered = factor(data_prs[, phen],
                          levels=c('nv012', 'notGE6adhd', 'nonimp', 'imp'),
                          ordered=T)
fmo_str = paste(prs, "~ ordered +",
               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
               sep="")
fito = lmer(as.formula(fmo_str), data=this_data, REML = FALSE)
print(Anova(fito))
posthoco = glht(fito, linfct=mcp(thresh0.00_inatt_GE6_wp05 = "Tukey"))
print(summary(posthoco))
```

lmer doesn't support that. We'd have to run an ordered logistic regression, like
note 009-ordered_logistc.Rmd. Let me focus on these for now. I can always just
plot the coefficients, if I want to check where they lie in the spectrum:

```r
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
sx = 'hi'
prs = 'ADHD_PRS0.000500'
min_sx = 6
thresh = 0
use_me = T#data_prs$isWNH

phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data_prs[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data_prs[, phen_slope]))
idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
data_prs[my_nvs[idx], phen] = 'nv012'
data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                               covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$sex = data_prs[use_me, 'sex.x']
tmp_covars = c(covars, 'sex')
fm_str = paste(prs, "~", phen, '+',
               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
               sep="")
fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
print(Anova(fit))
posthoc = glht(fit, linfct=mcp(thresh0.00_inatt_GE6_wp05 = "Tukey"))
print(summary(posthoc))

data_prs$ordered = factor(data_prs[, phen],
                          levels=c('nv012', 'notGE6adhd', 'nonimp', 'imp'),
                          ordered=T)
fmo_str = paste(prs, "~ ordered +",
               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
               sep="")
fito = lmer(as.formula(fmo_str), data=this_data, REML = FALSE)
print(Anova(fito))
posthoco = glht(fito, linfct=mcp(thresh0.00_inatt_GE6_wp05 = "Tukey"))
print(summary(posthoco))
```

Note that the actual model doesn't have the exact same stats, but I think it
should be fine if all I'm doing is plotting the coefficients.

```r
fm0_str = paste(prs, "~ 0 +", phen, '+',
               paste(tmp_covars, collapse='+'), '+(1|FAMID)',
               sep="")
fit0 = lmer(as.formula(fm0_str), data=this_data, REML = FALSE)
print(Anova(fit0))
posthoc = glht(fit, linfct=mcp(thresh0.50_hi_GE6_wp05 = "Tukey"))
print(summary(posthoc)) + (1|experiment), data = dataset)
tmp <- as.data.frame(confint(glht(model))$confint)
tmp$Comparison <- rownames(tmp)
ggplot(tmp, aes(x = Comparison, y = Estimate, ymin = lwr, ymax = upr)) +
  geom_errorbar() + geom_point()
```

## Ordinal logistic regression

OK, let's go back to the ordinal case then. We'll need to invert the model
again, and follow
https://cran.r-project.org/web/packages/ordinal/vignettes/clmm2_tutorial.pdf

Let's run the whole shabang:

```r
library(ordinal)

hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_all_PCsAgeSex_4groups_clmm2.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresholds = c(0, -.3)
    } else if (sx == 'hi') {
        thresholds = c(-.3, -.5)
    } else {
        thresholds = c(-.5, -1)
    }
    for (thresh in thresholds) {
        phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
        data_prs[, phen] = 'notGE6adhd'
        my_nvs = which(is.na(data_prs[, phen_slope]))
        idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
        data_prs[my_nvs[idx], phen] = 'nv012'
        data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
        data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
        data_prs[, phen] = factor(data_prs[, phen], ordered=F)
        data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
        use_me = T#data_prs$isWNH

        this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                                       covars)]
        this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
        this_data$ordered = factor(this_data[, phen],
                                   levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                                   ordered=T)
        this_data$sex = data_prs[use_me, 'sex.x']
        tmp_covars = c(covars, 'sex')
        this_data$FAMID = factor(this_data$FAMID)
        phen_res = c()
        for (prs in prs_var_names) {
            fm_str = paste("ordered ~", prs, '+',
                           paste(tmp_covars, collapse='+'),
                           sep="")
            fit = clmm2(as.formula(fm_str), random=FAMID, data=this_data, Hess=T)
            temp = c(summary(fit)$coefficients[prs, 'Pr(>|z|)'],
                     summary(fit)$logLik, summary(fit)$condHess)
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
colnames(hold)[1:3] = c('pval', 'loglik', 'hessian')
write.csv(hold, file=out_fname, row.names=F)
```

Results are not very consistent here. I get two hits for inatt at 0 threshold
but only for WNH. If looking at all, there are only some mixed hits for hi and
total. Let me see what heppens if I go down to 3 levels.

```r
library(ordinal)

hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_all_PCsAgeSex_3groups_clmm2.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresholds = c(0, -.3)
    } else if (sx == 'hi') {
        thresholds = c(-.3, -.5)
    } else {
        thresholds = c(-.5, -1)
    }
    for (thresh in thresholds) {
        phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
        data_prs[, phen] = 'nv'
        data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
        data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
        data_prs[, phen] = factor(data_prs[, phen], ordered=F)
        data_prs[, phen] = relevel(data_prs[, phen], ref='nv')
        use_me = T#data_prs$isWNH

        this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                                       covars)]
        this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
        this_data$ordered = factor(this_data[, phen],
                                   levels=c('nv', 'imp', 'nonimp'),
                                   ordered=T)
        this_data$sex = data_prs[use_me, 'sex.x']
        tmp_covars = c(covars, 'sex')
        this_data$FAMID = factor(this_data$FAMID)
        phen_res = c()
        for (prs in prs_var_names) {
            fm_str = paste("ordered ~", prs, '+',
                           paste(tmp_covars, collapse='+'),
                           sep="")
            fit = clmm2(as.formula(fm_str), random=FAMID, data=this_data, Hess=T)
            temp = c(summary(fit)$coefficients[prs, 'Pr(>|z|)'],
                     summary(fit)$logLik, summary(fit)$condHess)
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
colnames(hold)[1:3] = c('pval', 'loglik', 'hessian')
write.csv(hold, file=out_fname, row.names=F)
```

We get some interesting results here, so we need to start plotting them. The
results with 4 groups were actually a bit better than 3 groups. But we need to
make sure it all makes sense.

# 2020-01-28 08:20:44

Philip asked me to run the lme version of this as well, with the model inverted.
Let's see what we get:

```r
library(nlme)

hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_WNH_PCsAgeSex_4groups_lme.csv'
for (sx in c('inatt', 'hi', 'total')) {
    min_sx = 6
    if (sx == 'inatt') {
        thresholds = c(0, -.3)
    } else if (sx == 'hi') {
        thresholds = c(-.3, -.5)
    } else {
        thresholds = c(-.5, -1)
    }
    for (thresh in thresholds) {
        phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
        data_prs[, phen] = 'notGE6adhd'
        my_nvs = which(is.na(data_prs[, phen_slope]))
        idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
        data_prs[my_nvs[idx], phen] = 'nv012'
        data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
        data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
        data_prs[, phen] = factor(data_prs[, phen], ordered=F)
        data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
        use_me = data_prs$isWNH

        this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                                       covars)]
        this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
        this_data$ordered = factor(this_data[, phen],
                                   levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                                   ordered=T)
        this_data$sex = data_prs[use_me, 'sex.x']
        tmp_covars = c(covars, 'sex')
        this_data$FAMID = factor(this_data$FAMID)
        phen_res = c()
        for (prs in prs_var_names) {
            fm_str = paste(prs, " ~ ordered +",
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
        phen_res$predictor = prs_var_names
        phen_res$outcome = phen
        hold = rbind(hold, phen_res)
    }
}
colnames(hold)[6:8] = c('loglik', 'AIC', 'BIC')
write.csv(hold, file=out_fname, row.names=F)
```

We're still doing well. It's just a matter of figuring out which model looks
nicer, how we want to plot it, etc.

## Plotting clmm2 model

I'll focus on one of the good results for clmm2 and plot it. Then we can make
plots for all good results. In the next section I'll do the same for the lme
results.

```r
library(ordinal)
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
sx = 'hi'
min_sx = 6
thresh = -.5
prs = 'ADHD_PRS0.001000'

phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data_prs[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data_prs[, phen_slope]))
idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
data_prs[my_nvs[idx], phen] = 'nv012'
data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
use_me = T#data_prs$isWNH
this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                               covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$ordered = factor(this_data[, phen],
                           levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                           ordered=T)
this_data$sex = data_prs[use_me, 'sex.x']
tmp_covars = c(covars, 'sex')
this_data$FAMID = factor(this_data$FAMID)
fm_str = paste("ordered ~", prs, '+',
               paste(tmp_covars, collapse='+'),
               sep="")
fit = clmm2(as.formula(fm_str), random=FAMID, data=this_data, Hess=T)
```

Cannot really find a good way to do this... apparently the function to predict
the probability for different classes is failing... let's play with lme then.

```r
# newdat <- expand.grid(sex=unique(this_data$sex),
#                       ADHD_PRS0.001000=c(min(this_data$ADHD_PRS0.001000),
#                                          max(this_data$ADHD_PRS0.001000)),
#                       base_age=c(min(this_data$base_age), max(this_data$base_age)),
#                       PC01=c(min(this_data$PC01), max(this_data$PC01)),
#                       PC02=c(min(this_data$PC02), max(this_data$PC02)),
#                       PC03=c(min(this_data$PC03), max(this_data$PC03)),
#                       PC04=c(min(this_data$PC04), max(this_data$PC04)),
#                       PC05=c(min(this_data$PC05), max(this_data$PC05)),
#                       PC06=c(min(this_data$PC06), max(this_data$PC06)),
#                       PC07=c(min(this_data$PC07), max(this_data$PC07)),
#                       PC08=c(min(this_data$PC08), max(this_data$PC08)),
#                       PC09=c(min(this_data$PC09), max(this_data$PC09)),
#                       PC10=c(min(this_data$PC10), max(this_data$PC10)))
newdat <- expand.grid(sex=unique(this_data$sex),
                      ordered=unique(this_data$ordered),
                      FAMID=unique(this_data$FAMID),
                      base_age=c(min(this_data$base_age), max(this_data$base_age)),
                      PC01=c(min(this_data$PC01), max(this_data$PC01)),
                      PC02=c(min(this_data$PC02), max(this_data$PC02)),
                      PC03=c(min(this_data$PC03), max(this_data$PC03)),
                      PC04=c(min(this_data$PC04), max(this_data$PC04)),
                      PC05=c(min(this_data$PC05), max(this_data$PC05)),
                      PC06=c(min(this_data$PC06), max(this_data$PC06)),
                      PC07=c(min(this_data$PC07), max(this_data$PC07)),
                      PC08=c(min(this_data$PC08), max(this_data$PC08)),
                      PC09=c(min(this_data$PC09), max(this_data$PC09)),
                      PC10=c(min(this_data$PC10), max(this_data$PC10)))
library(nlme)
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHD_')]
covars = c(sapply(1:10, function(x) sprintf('PC%02d', x)), 'base_age')
sx = 'hi'
min_sx = 6
thresh = -.5
prs = 'ADHD_PRS0.001000'

phen_slope = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
phen = sprintf('thresh%.2f_%s_GE%d_wp05', abs(thresh), sx, min_sx)
data_prs[, phen] = 'notGE6adhd'
my_nvs = which(is.na(data_prs[, phen_slope]))
idx = data_prs[my_nvs, 'base_inatt'] <= 2 & data_prs[my_nvs, 'base_hi'] <= 2
data_prs[my_nvs[idx], phen] = 'nv012'
data_prs[which(data_prs[, phen_slope] < thresh), phen] = 'imp'
data_prs[which(data_prs[, phen_slope] >= thresh), phen] = 'nonimp'
data_prs[, phen] = factor(data_prs[, phen], ordered=F)
data_prs[, phen] = relevel(data_prs[, phen], ref='nv012')
use_me = T#data_prs$isWNH
this_data = data_prs[use_me, c(phen, 'FAMID', prs_var_names,
                               covars)]
this_data[, 3:ncol(this_data)] = scale(this_data[, 3:ncol(this_data)])
this_data$ordered = factor(this_data[, phen],
                           levels=c('nv012', 'notGE6adhd', 'imp', 'nonimp'),
                           ordered=T)
this_data$sex = data_prs[use_me, 'sex.x']
tmp_covars = c(covars, 'sex')
this_data$FAMID = factor(this_data$FAMID)
fm_str = paste(prs, " ~ ordered +",
               paste(tmp_covars, collapse='+'), sep="")
fit = lme(as.formula(fm_str), ~1|FAMID, data=this_data)
```

Let's not spend much time here. We'll worry about the plots later. Let's go for
imaging now.

Chatting with Philip, he suggested I should look into the cubic fits as well,
using the idea that remitted might have protective alleles. For the brain
analysis, we shouldn't really strict ourselves to the ordered model, as we are
not sure of what's going show. But we could say it's based on our previous
framework. So, maybe start with the regular ANOVA, then try to look for it in
ordered fashion, and finally look at the ordinal classes as the outcome.

Philip also suggested I should look at even more strict PRS thresholds, as it
sounds like our results are in that lower side.


# TODO
* which ones can we keep? which ones don't violate the assumptions and have good
  pairwise comparisons?
* try persistence PRS
* continue work on glmer model just for robustness 