# 2021-06-08 05:58:17

Let's see what we get running the univariate analysis on the same data we sent
to Gang. First, no outlier removal. I'll even use the same code he did all the
way up until I can't anymore.

## All data, univariate, random slope and intercept

```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
dat$famID <- as.factor(dat$famID) # 522 families
dat$SID <- as.factor(dat$SID) # 925 subjects
dd <- dat[complete.cases(dat), ] # 20053 observations remaining
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction

fm_str = 'FA ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (age_scan|famID)+(age_scan|SID)'

library(lme4)
tracts = unique(dd$tract)

this_data = dd[dd$tract == tracts[1], ]
fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)
```

I got an error running a single tract:

```
r$> fit = lmer(as.formula(fm_str), data=this_data, REML = FALSE)                          
Error: number of observations (=1823) <= number of random effects (=1834) for term (age_scan | SID); the random-effects parameters and the residual variance (or scale parameter) are probably unidentifiable
```

I wonder if adding this term in the Bayesian model could be causing those
results? Worth re-running it with and without it.

## All data, univariate, random intercept only

```r
library(car)

fm_str = '%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)'

tracts = unique(dd$tract)

res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_all_rndIntercept.csv',
          row.names=F)
```

```
r$> # quick assessment 
    res.df[res.df$pval_grp < .05, 1:3]
        prop     tract             pval_grp
temp.6    FA right_ifo   0.0128237473046893
temp.8    FA right_slf   0.0243673966451264
temp.23   RD  left_ifo  0.00694181538202031
temp.28   RD right_ifo   0.0014938787235355
temp.30   RD right_slf 0.000256452532878001
temp.31   RD right_unc  0.00163695177969821
```

## Clean data, univariate, random intercept only

Let's do some cleaning and see how that affects the univariate results. I want
to bring in the master spreadsheet so we can keep it only to 60 volumes scans.

```r
ncr = read.csv('~/data/all_DTI_tracts_03222021.csv')
ncr = ncr[ncr$maskid != 1090, ]  # remove one scan we didn't finish processing
ncr$seq = 'adult'
ncr[which(ncr$numVolumes < 80), 'seq'] = 'child'
kid_scans = ncr[ncr$seq == 'child', 'maskid']

dd.clean = dd[dd$maskid %in% kid_scans, ] # 14784 observations remaining

# some histograms to remove QC outliers
quartz()
par(mfrow = c(1,3))
hist(dd.clean$missingVolumes, breaks=50)
hist(dd.clean$norm.trans, breaks=50)
hist(dd.clean$norm.rot, breaks=50)
```

![](images/2021-06-08-06-50-53.png) 

A bit crude, but for now this will do:

```r
keep_idx = (dd.clean$missingVolumes <= 4 & dd.clean$norm.trans <= 2.5 &
            dd.clean$norm.rot <= .04)
dd.clean = dd.clean[keep_idx,]  # 10824 observations remaining
```

Let's remove based on the data quality as well:

```r
par(mfrow=c(3, 4))
for (v in tracts) {
    hist(dd.clean[dd.clean$tract == v, 'FA'], breaks=25, main=v)
}
```

![](images/2021-06-08-06-58-05.png)

```r
keep_scans = dd.clean[dd.clean$tract=='left_unc' & dd.clean$FA > .25, 'maskid']
dd.clean = dd.clean[dd.clean$maskid %in% keep_scans,]  # 10703 observations

par(mfrow=c(3, 4))
for (v in tracts) {
    hist(dd.clean[dd.clean$tract == v, 'FA'], breaks=25, main=v)
}
```

![](images/2021-06-08-07-03-13.png)

Maybe just another one:

```r
keep_scans = dd.clean[dd.clean$tract=='right_unc' & dd.clean$FA > .22, 'maskid']
dd.clean = dd.clean[dd.clean$maskid %in% keep_scans,]  # 10692 observations
keep_scans = dd.clean[dd.clean$tract=='left_cst' & dd.clean$FA > .29, 'maskid']
dd.clean = dd.clean[dd.clean$maskid %in% keep_scans,]  # 10670 observations

par(mfrow=c(3, 4))
for (v in tracts) {
    hist(dd.clean[dd.clean$tract == v, 'FA'], breaks=25, main=v)
}
```

![](images/2021-06-08-07-06-07.png)

That looks better. Let's save it and re-run all analyses:

```r
write.csv(dd.clean, file='~/data/long_data_for_gang_with_FAMID_clean.csv',
          row.names=F)
```

```r
library(lme4)
library(car)

dat <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
dat$famID <- as.factor(dat$famID) # 277 families
dat$SID <- as.factor(dat$SID) # 415 subjects
dd <- dat[complete.cases(dat), ] # 10670 observations remaining
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction

tracts = unique(dd$tract)

fm_str = '%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)'

res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_clean_rndIntercept.csv',
          row.names=F)
```

```
r$> # quick assessment 
    res.df[res.df$pval_grp < .05, 1:3]
        prop     tract            pval_grp
temp      FA  left_cst  0.0146252351616683
temp.2    FA  left_ilf  0.0133132661423609
temp.7    FA right_ilf  0.0141072577552538
temp.13   AD  left_ilf 0.00051669816065722
temp.29   RD right_ilf  0.0435892170760677
temp.31   RD right_unc  0.0495423982674764
```

Results change a bit from using all data points. 

## Clean data, univariate, random intercept and slope

However, now we can actually run the random slopes as well (lots of singular fits, though):

```r
fm_str = '%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (age_scan|famID/SID)'

res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_clean_rndInterceptSlope.csv',
          row.names=F)
print(res.df[res.df$pval_grp < .05, 1:3])
```

```
        prop    tract           pval_grp
temp      FA left_cst 0.0215027658358196
temp.2    FA left_ilf 0.0199250826036848
temp.13   AD left_ilf 0.0103914545675075
```


# All data, Bayesian model, random intercept only

Let's then use Gang's code to run the Bayesian model. First, testing how
estimating the random slope affects the results:

```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
dat$famID <- as.factor(dat$famID) # 522 families
dat$SID <- as.factor(dat$SID) # 925 subjects
dd <- dat[complete.cases(dat), ] # 20053 observations remaining
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction
 
mds = c('FA', 'AD', 'RD')
fixed_vars = colnames(dd)[c(1:3, 5:9, 13, 15, 16)]
# write out 3 tables for Bayesian modeling
for (m in mds) {
    fname = sprintf('~/data/bayesian/%s_all.tbl', m)
    write.table(dd[, c(fixed_vars, m)], file = fname,
            quote = FALSE, row.names = FALSE)
}
```

Copy the new tables to BW, and make sure the right packages are installed:

```bash
# bw
module load R
conda acitvate radian
radian
```

```r
# bw
install.packages('brms')
install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/',
                 getOption('repos')))
# had to load gcc module first
install_cmdstan(dir='~/', cores=2)
```

I had to install my own version of CmdStan because the HPC module did not expose
the Makefile for cmdstanr to extract the version. 

Then, run the R code for the Bayesian model:

```r
m = 'FA'
fm_str = sprintf('%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)+ (grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes|tract)', m)
out_fname = sprintf('~/data/bayesian/%s_all_rndIntercept.RData', m)
in_fname = sprintf('~/data/bayesian/%s_all.tbl', m)

dat <- read.table(in_fname, header=TRUE)
dat$famID <- as.factor(dat$famID)
library('brms')
ncpu = future::availableCores()
options(mc.cores = ncpu)
require('cmdstanr')
set_cmdstan_path('~/cmdstan-2.27.0')
nit <- 1000
nch <- 4 # number of chains
nthreads = floor(ncpu / nch)
set.seed(1234)
fm <- brm(as.formula(fm_str),
          data=dat, family='student', chains = nch, iter=nit,
          control = list(adapt_delta = 0.99, max_treedepth = 20),
          backend = "cmdstanr", threads = threading(nthreads))
save.image(file=out_fname)
```

So, as Gang mentioned, this is shaping up to take forever. Now, I need to figure
out some ways to optimize it. For now, let's leave it running for several
different models, and we can try to optimize it later.

Just put this into a swarm file:

```
Rscript FA_all_int.R | tee -a FA_all_int.log
Rscript AD_all_int.R | tee -a AD_all_int.log
Rscript RD_all_int.R | tee -a RD_all_int.log
```

And swarm it:

```bash
cd ~/data/bayesian
swarm --job-name all_int -f swarm.all_int -t 32 -g 10 -m R --time=10-00:00:00
```

# Clean data, Bayesian model, random intercept only

```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
dat$famID <- as.factor(dat$famID) # 277 families
dat$SID <- as.factor(dat$SID) # 415 subjects
dd <- dat[complete.cases(dat), ] # 10670 observations remaining
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction
 
mds = c('FA', 'AD', 'RD')
fixed_vars = colnames(dd)[c(1:3, 5:9, 13, 15, 16)]
# write out 3 tables for Bayesian modeling
for (m in mds) {
    fname = sprintf('~/data/bayesian/%s_clean.tbl', m)
    write.table(dd[, c(fixed_vars, m)], file = fname,
            quote = FALSE, row.names = FALSE)
}
```

```r
m = 'FA'
fm_str = sprintf('%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)+ (grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes|tract)', m)
out_fname = sprintf('~/data/bayesian/%s_clean_rndIntercept.RData', m)
in_fname = sprintf('~/data/bayesian/%s_clean.tbl', m)

dat <- read.table(in_fname, header=TRUE)
dat$famID <- as.factor(dat$famID)
library('brms')
ncpu = future::availableCores()
options(mc.cores = ncpu)
require('cmdstanr')
set_cmdstan_path('~/cmdstan-2.27.0')
nit <- 1000
nch <- 4 # number of chains
nthreads = floor(ncpu / nch)
print(nthreads)
set.seed(1234)
fm <- brm(as.formula(fm_str),
          data=dat, family='student', chains = nch, iter=nit,
          control = list(adapt_delta = 0.99, max_treedepth = 20),
          backend = "cmdstanr", threads = threading(nthreads))
save.image(file=out_fname)
```

# Clean data, Bayesian model, random intercept and slope

```r
m = 'FA'
fm_str = sprintf('%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (age_scan|famID/SID)+ (grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes|tract)', m)
out_fname = sprintf('~/data/bayesian/%s_clean_rndInterceptSlope.RData', m)
in_fname = sprintf('~/data/bayesian/%s_clean.tbl', m)

dat <- read.table(in_fname, header=TRUE)
dat$famID <- as.factor(dat$famID)
library('brms')
ncpu = future::availableCores()
options(mc.cores = ncpu)
require('cmdstanr')
set_cmdstan_path('~/cmdstan-2.27.0')
nit <- 1000
nch <- 4 # number of chains
nthreads = floor(ncpu / nch)
print(nthreads)
set.seed(1234)
fm <- brm(as.formula(fm_str),
          data=dat, family='student', chains = nch, iter=nit,
          control = list(adapt_delta = 0.99, max_treedepth = 20),
          backend = "cmdstanr", threads = threading(nthreads))
save.image(file=out_fname)
```

## Optimizing?

Let's try some new things to try to speed this up:

```r
m = 'FA'
fm_str = sprintf('%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)+ (grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes|tract)', m)
out_fname = sprintf('~/data/bayesian/%s_all_rndIntercept_opt.RData', m)
in_fname = sprintf('~/data/bayesian/%s_all.tbl', m)

dat <- read.table(in_fname, header=TRUE)
dat$famID <- as.factor(dat$famID)
library('brms')
ncpu = future::availableCores()
options(mc.cores = ncpu)
require('cmdstanr')
set_cmdstan_path('~/cmdstan-2.27.0')
nit <- 1000
nch <- 4 # number of chains
nthreads = floor(ncpu / nch)
set.seed(1234)
fm <- brm(as.formula(fm_str),
          data=dat, family='student', chains = nch, iter=nit,
          control = list(adapt_delta = 0.99, max_treedepth = 20),
          backend = "rstan")
save.image(file=out_fname)
```

It's taking forever in the warmup step. But that's not just in rstan backedn...
it's happening all the time.


# 2021-06-09 09:28:58

Actually, the clean runs finished between and 6 and 12h. The runs with all data
are still going though. In the meanwhile, let's collect the clean run results.

```r
#local
require(brms)
require(ggplot2)
source('~/research_code/stackOrdered.R')

for (it in c('', 'Slope')) {
    for (m in c('FA', 'AD', 'RD')) {
        froot = sprintf('%s_cleanLong_rndIntercept%s', m, it)
        cat(froot, '\n')
        load(sprintf('~/data/bayesian/%s.RData', froot))
        ns <- 2000
        pe <- fixef(fm, summary = FALSE) # Population-Level Estimates
        ge <- ranef(fm, summary = FALSE)
        ps <- apply(ge[['tract']][,,'GA'], 2, '+', pe[,'GA'])
        p = stackPlot(ps, range(ps), froot)
        ggsave(file = sprintf('~/data/bayesian/%s.png', froot),
               width=8, height=10, units = 'in', dpi = 300)
    }
}
```

As I'm putting the presentation together, I noticed we could run only the actual
longitudinal subjects too. We can do it for all and clean.

## All longitudinal data, univariate, random slope and intercept

```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
dat$SID <- as.factor(dat$SID) # 925 subjects
dat$famID <- as.factor(dat$famID) # 522 families
dd <- dat[complete.cases(dat), ] # 20053 observations remaining
long_subjs = names(table(dd$SID))[table(dd$SID)>=22]
dd = dd[dd$SID %in% long_subjs, ]
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction

fm_str = 'FA ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (age_scan|famID)+(age_scan|SID)'

library(lme4)
library(car)
res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_allLong_rndInterceptSlope.csv',
          row.names=F)
print(res.df[res.df$pval_grp < .05, 1:3])
```

Got several singular issues, but it ran without crashing errors. Nothing with
nominal p < .05.


## All longitudinal data, univariate, random intercept only

```r
library(car)

fm_str = '%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)'

tracts = unique(dd$tract)

res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_allLong_rndIntercept.csv',
          row.names=F)
print(res.df[res.df$pval_grp < .05, 1:3])
```

A few warnings that the model didn't converge, but lots of results p < .05:

```
        prop     tract             pval_grp
temp      FA  left_cst   0.0316321241337544
temp.1    FA  left_ifo   0.0337826201828238
temp.2    FA  left_ilf   0.0433848055434377
temp.3    FA  left_slf   0.0231147443920181
temp.4    FA  left_unc   0.0326683982040515
temp.5    FA right_cst   0.0130615344944267
temp.6    FA right_ifo   0.0251524701096188
temp.7    FA right_ilf  0.00549799020703181
temp.8    FA right_slf   0.0238289894786366
temp.9    FA right_unc   0.0178932267213103
temp.10   FA        cc   0.0396473454287343
temp.23   RD  left_ifo   0.0291829587974558
temp.25   RD  left_slf  0.00448023186489197
temp.26   RD  left_unc   0.0123828097616283
temp.28   RD right_ifo   0.0114010177064767
temp.29   RD right_ilf 0.000127576266537002
temp.30   RD right_slf  0.00516034212602674
temp.31   RD right_unc 0.000502365259487459
```

## Clean data, longitudinal only, univariate, random intercept only

```r
library(lme4)
library(car)

dat <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
dat$famID <- as.factor(dat$famID) # 277 families
dat$SID <- as.factor(dat$SID) # 415 subjects
dd <- dat[complete.cases(dat), ] # 10670 observations remaining
long_subjs = names(table(dd$SID))[table(dd$SID)>=22]
dd = dd[dd$SID %in% long_subjs, ]
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction

tracts = unique(dd$tract)

fm_str = '%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (1|famID/SID)'

res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_cleanLong_rndIntercept.csv',
          row.names=F)
print(res.df[res.df$pval_grp < .05, 1:3])
```

```
        prop     tract            pval_grp
temp.13   AD  left_ilf 0.00579010277909848
temp.31   RD right_unc  0.0200698648576646
```

## Clean data, univariate, random intercept and slope

```r
fm_str = '%s ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+ (age_scan|famID/SID)'

res = c()
for (md in c('FA', 'AD', 'RD')) {
    fm = as.formula(sprintf(fm_str, md))
    for (tr in tracts) {
        this_data = dd[dd$tract == tr, ]
        fit = lmer(fm, data=this_data, REML = FALSE)
        p = Anova(fit)
        temp = c(md, tr, p[1,3], summary(fit)$AIC[1:4], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_grp', 'AIC', 'BIC', 'logLik',
                     'deviance', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_cleanLong_rndInterceptSlope.csv',
          row.names=F)
print(res.df[res.df$pval_grp < .05, 1:3])
```

Lots of singular warnings and model not converging. Still...

```
        prop    tract           pval_grp
temp      FA left_cst 0.0430661303801234
temp.13   AD left_ilf 0.0403563164388045
```

# All longitudinal data, Bayesian model, random intercept only


```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
dat$famID <- as.factor(dat$famID) # 522 families
dat$SID <- as.factor(dat$SID) # 925 subjects
dd <- dat[complete.cases(dat), ] # 20053 observations remaining
long_subjs = names(table(dd$SID))[table(dd$SID)>=22]
dd = dd[dd$SID %in% long_subjs, ]
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction
 
mds = c('FA', 'AD', 'RD')
fixed_vars = colnames(dd)[c(1:3, 5:9, 13, 15, 16)]
# write out 3 tables for Bayesian modeling
for (m in mds) {
    fname = sprintf('~/data/bayesian/%s_allLong.tbl', m)
    write.table(dd[, c(fixed_vars, m)], file = fname,
            quote = FALSE, row.names = FALSE)
}
```

# Clean longitudinal data, Bayesian model, random intercept only

```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
dat$famID <- as.factor(dat$famID) # 277 families
dat$SID <- as.factor(dat$SID) # 415 subjects
dd <- dat[complete.cases(dat), ] # 10670 observations remaining
long_subjs = names(table(dd$SID))[table(dd$SID)>=22]
dd = dd[dd$SID %in% long_subjs, ]
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction
 
mds = c('FA', 'AD', 'RD')
fixed_vars = colnames(dd)[c(1:3, 5:9, 13, 15, 16)]
# write out 3 tables for Bayesian modeling
for (m in mds) {
    fname = sprintf('~/data/bayesian/%s_cleanLong.tbl', m)
    write.table(dd[, c(fixed_vars, m)], file = fname,
            quote = FALSE, row.names = FALSE)
}
```

These shouldn't take very long, because it's way less data than before.




# TODO
 * I saw this: "Ensure that the data is randomly sorted such that consecutive
   subsets of the data are roughly of the same computational effort.", in
   https://cran.r-project.org/web/packages/brms/vignettes/brms_threading.html.
   Maybe something to keep in mind?





# Gang's original code:

```r
dat <- read.csv('long_data_for_gang_with_FAMID.csv', header=T)  # read in table Gustavo provided
dat$famID <- as.factor(dat$famID) # 522 families
dat$SID <- as.factor(dat$SID) # 925 subjects
#dd <- dat[!is.na(dat),] # remove rows with NAs
dd <- dat[complete.cases(dat), ] # 20053 observations remaining
dd$grp <- ifelse(dd$group=='NV', -0.5, 0.5)  # dummy code group
dd$GA <- dd$grp*dd$age_scan   # group-by-age interaction
 
# write out 3 tables for Bayesian modeling
write.table(dd[,c(1:3,5:9,13,15,16,10)], file = "FA.tbl", quote = FALSE, row.names = FALSE)
write.table(dd[,c(1:3,5:9,13,15,16,11)], file = "AD.tbl", quote = FALSE, row.names = FALSE)
write.table(dd[,c(1:3,5:9,13,15,16,12)], file = "RD.tbl", quote = FALSE, row.names = FALSE)
 
###########  Perform Bayesian modeling for RD ###########
## required installations of R libraries
## 1) install the R package 'brms'
## 2) install the R package 'cmdstan' using the following command in R:
install.packages('cmdstanr', repos = c('https://mc-stan.org/r-packages/', getOption('repos')))
## 3) install the R package  'cmdstanr' in your home directory by following the instruction here:
    https://mc-stan.org/cmdstanr/articles/cmdstanr.html
```

A. Save the following R code as a text file and call it, for example, RDr.R (assuming 24 CPUs available on the computer)

```r 
dat <- read.table('RD.tbl', header=TRUE)
dat$famID <- as.factor(dat$famID)
library('brms')
options(mc.cores = parallel::detectCores())
require('cmdstanr')
set_cmdstan_path('~/cmdstan')
nit <- 1000
nch <- 4 # number of chains
set.seed(1234)
fm <- brm(RD ~ grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes+
               (age_scan|famID)+(age_scan|SID)+
               (grp+age_scan+GA+scanner_update+sex+norm.trans+norm.rot+missingVolumes|tract),
          data=dat, family='student', chains = nch, iter=nit,
          control = list(adapt_delta = 0.99, max_treedepth = 20),
          backend = "cmdstanr", threads = threading(6))
save.image(file='RDr.RData')
```

B. Run the script RDr.R by typing the following at the terminal (may take a few days):

```bash 
nohup R CMD BATCH –no-save RDr.R RDr.diary &
```

C. Download the following plotting function and put it somewhere (e.g., home directory):
 
https://afni.nimh.nih.gov/sscc/staff/gangc/pub/stackOrdered.R
 
D. Once B is done, obtain the group-by-age effect by running the following in R:

```r 
load('RDr.RData')
require(brms)
require(ggplot2)
ns <- 2000
pe <- fixef(fm, summary = FALSE) # Population-Level Estimates
ge <- ranef(fm, summary = FALSE)
ps <- apply(ge[['tract']][,,'GA'], 2, '+', pe[,'GA'])
source('~/stackOrdered.R')
stackPlot(ps, range(ps), 'RD-GA')
```