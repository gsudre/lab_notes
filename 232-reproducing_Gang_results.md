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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_all_rndIntercept.csv',
          row.names=F)
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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_clean_rndIntercept.csv',
          row.names=F)
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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_clean_rndInterceptSlope.csv',
          row.names=F)
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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_allLong_rndInterceptSlope.csv',
          row.names=F)
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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_allLong_rndIntercept.csv',
          row.names=F)
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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_cleanLong_rndIntercept.csv',
          row.names=F)
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
        temp = c(md, tr, p['GA',3], summary(fit)$AIC[1:4],
                 summary(fit)$coefficients['GA',1:2], fm_str)
        res = rbind(res, temp)
    }
}
res.df = data.frame(res)
colnames(res.df) = c('prop', 'tract', 'pval_GA', 'AIC', 'BIC', 'logLik',
                     'deviance', 'Estimate', 'StdError', 'formula')
write.csv(res.df, file='~/data/bayesian/univar_cleanLong_rndInterceptSlope.csv',
          row.names=F)
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

# 2021-06-14 17:21:40

Let's also create some data density plots to show how it changes as we cleaned
the data.

```r
dat <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
dat$status = 'all'
datc <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
datc$status = 'clean'
datc$grp = NULL
datc$GA=NULL
df = rbind(dat, datc)

quartz()
myplots = list()
cnt = 1
for (m in c('FA', 'AD', 'RD')) {
    df$x = df[, m]
    p = ggplot(df, aes(x=x, fill=status)) + geom_density(alpha=0.4) +
        labs(x=m)
    myplots[[cnt]] = p
    cnt = cnt + 1
}

dat <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
dat$status = 'all'
datl <- read.csv('~/data/long_data_for_gang_with_FAMID.csv', header=T)
datl$status = 'all_long'
long_subjs = names(table(datl$SID))[table(datl$SID)>=22]
datl = datl[datl$SID %in% long_subjs, ]
df = rbind(dat, datl)
for (m in c('FA', 'AD', 'RD')) {
    df$x = df[, m]
    p = ggplot(df, aes(x=x, fill=status)) + geom_density(alpha=0.4) +
        labs(x=m)
    myplots[[cnt]] = p
    cnt = cnt + 1
}

dat <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
dat$status = 'clean'
datl <- read.csv('~/data/long_data_for_gang_with_FAMID_clean.csv', header=T)
datl$status = 'clean_long'
long_subjs = names(table(datl$SID))[table(datl$SID)>=22]
datl = datl[datl$SID %in% long_subjs, ]
df = rbind(dat, datl)
for (m in c('FA', 'AD', 'RD')) {
    df$x = df[, m]
    p = ggplot(df, aes(x=x, fill=status)) + geom_density(alpha=0.4) +
        labs(x=m)
    myplots[[cnt]] = p
    cnt = cnt + 1
}

library(ggpubr)
ggarrange(plotlist=myplots, nrow=3, ncol=3)
```

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

