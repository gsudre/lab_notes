# 2020-03-14 05:31:40

Philip sent me guidelines for the project on Slack yesterday. Let's start
running it. I first removed the data matrix to an RDS file, and we'll go from
there.

The approach will be to run lm() within ACC and within caudate first, filtering
the covariates at p< .1, p<.05, and using stepAIC. That gives me 6 tabs of
results. Then, I'll do the same using lme() and the brain region as the random
term, for another 3 tabs. Then we can check for consistencies across result
tabs.

```r
myregion = 'ACC'
pthresh = .05
keep_str = 'Diagnosis*Age'

data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data$substance_group = as.factor(data$substance_group)
data$batch = as.factor(data$batch)
# no column names as numbers!
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names
data = data[data$Region==myregion, ]

# dependent
dep_vars = colnames(data)[34:ncol(data)]
# keep these regardless of significance
keep_vars = c(keep_str)
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'pH', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex')
# spit out the results
out_fname = sprintf('~/tmp/res_%s_pLT%.02f_%s.csv', myregion, pthresh,
                    gsub(pattern='\\*',replacement='',x=keep_str))

hold = c()
for (dp in 1:length(dep_vars)) {
    if (dp %% 50 == 0) {
        print(sprintf('%d of %d (%s)', dp, length(dep_vars), out_fname))
    }

    dep_var = dep_vars[dp]
    fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                   paste(test_vars, collapse='+'), sep="")
    fit = lm(as.formula(fm_str), data=data)
    res = summary(fit)$coefficients
    # filtering variables
    sig_vars = c()
    for (v in 1:length(test_vars)) {
        # rows in results table that correspond to the screened variable
        var_rows = which(grepl(rownames(res),
                         pattern=sprintf('^%s', test_vars[v])))
        for (r in var_rows) {
            if (res[r, 'Pr(>|t|)'] < pthresh) {
                sig_vars = c(sig_vars, test_vars[v])
            }
        }
    }
    # factors might get added several times, so here we clean it up
    sig_vars = unique(sig_vars)
    if (length(sig_vars) > 0) {
        clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                       paste(sig_vars, collapse='+'), sep="")
    } else {
        clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), sep="")
    }
    # new model
    clean_fit = lm(as.formula(clean_fm_str), data=data)
    res = data.frame(summary(clean_fit)$coefficients)
    # remove intercept
    res = res[2:nrow(res),]
    res$dep_var = dep_var
    res$formula = clean_fm_str
    res$orig_formula = fm_str
    res$predictor = rownames(res)
    hold = rbind(hold, res)
}
write.csv(hold, file=out_fname, row.names=F)
```

I did some tests running stepAIC but the funciton kept dying. Need to explore
that a bit further. For now, let's focus on the p-value exclusion first.

For the model that adds in region we'll end up with two measurements per
subject. So, we'll run lme():

```r
library(nlme)
pthresh = .05
keep_str = 'Diagnosis*Region'

data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data$substance_group = as.factor(data$substance_group)
data$batch = as.factor(data$batch)
data$hbcc_brain_id = as.factor(data$hbcc_brain_id)

# no column names as numbers!
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names
# dependent
dep_vars = colnames(data)[34:ncol(data)]
# keep these regardless of significance
keep_vars = c(keep_str)
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'pH', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'Age')
# spit out the results
out_fname = sprintf('~/tmp/res_pLT%.02f_%s.csv', pthresh,
                    gsub(pattern='\\*',replacement='',x=keep_str))

hold = c()
for (dp in 1:length(dep_vars)) {
    if (dp %% 50 == 0) {
        print(sprintf('%d of %d (%s)', dp, length(dep_vars), out_fname))
    }

    dep_var = dep_vars[dp]
    fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                   paste(test_vars, collapse='+'), sep="")
    fit = try(lme(as.formula(fm_str), ~1|hbcc_brain_id, data=data, na.action=na.omit))
    if (length(fit) > 1) {
        res = summary(fit)$tTable
        # filtering variables
        sig_vars = c()
        for (v in 1:length(test_vars)) {
            # rows in results table that correspond to the screened variable
            var_rows = which(grepl(rownames(res),
                            pattern=sprintf('^%s', test_vars[v])))
            for (r in var_rows) {
                if (res[r, 'p-value'] < pthresh) {
                    sig_vars = c(sig_vars, test_vars[v])
                }
            }
        }
        # factors might get added several times, so here we clean it up
        sig_vars = unique(sig_vars)
        if (length(sig_vars) > 0) {
            clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                        paste(sig_vars, collapse='+'), sep="")
        } else {
            clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), sep="")
        }
        # new model
        clean_fit = try(lme(as.formula(clean_fm_str), ~1|hbcc_brain_id, data=data,
                        na.action=na.omit))
        if (length(clean_fit) > 1) {
            res = data.frame(summary(clean_fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = clean_fm_str
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        } else {
            res = data.frame(summary(fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = NA
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        }
        hold = rbind(hold, res)
    }
}
write.csv(hold, file=out_fname, row.names=F)
```

So, in the end I had to change from the original approach. Now, we have only p <
.05 and p<.1 to select the variables. I'll play with stepAIC later if necessary.
I then tried lme models if using both regions, but only lm if doing it within
region. I also played with Diagnosis or Diagnosis*Age for lm (keeping Age as
covariate candidate in the former), and Diagnosis*Region, Diagnosis*Age for lme,
but keeping Region as fixed covariate and Age as fitereable when appropriate.

I'm compiling them into 2 different Excel sheets, with different tabs
each.

Maybe it's easier to do the filtering in R:

```r
res = read.csv('res_ACC_pLT0.05_Diagnosis.csv')
res = res[res$predictor=='DiagnosisControl',]
p = res[, 'Pr...t..']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests p < .05: %d', sum(p<.05)))
print(sprintf('Tests p < .01: %d', sum(p<.01)))
print(sprintf('Tests q < .05: %d', sum(p2<.05)))
print(sprintf('Tests q < .1: %d', sum(p2<.1)))
```

```
> res = read.csv('res_ACC_pLT0.05_Diagnosis.csv')
[1] "Tests p < .05: 2215"
[1] "Tests p < .01: 581"
[1] "Tests q < .05: 2"
[1] "Tests q < .1: 2"
> res = read.csv('res_ACC_pLT0.10_Diagnosis.csv')
[1] "Tests p < .05: 2590"
[1] "Tests p < .01: 696"
[1] "Tests q < .05: 3"
[1] "Tests q < .1: 3"
> res = read.csv('res_Caudate_pLT0.05_Diagnosis.csv')
[1] "Tests p < .05: 2619"
[1] "Tests p < .01: 673"
[1] "Tests q < .05: 0"
[1] "Tests q < .1: 0"
> res = read.csv('res_Caudate_pLT0.10_Diagnosis.csv')
[1] "Tests p < .05: 2947"
[1] "Tests p < .01: 779"
[1] "Tests q < .05: 0"
[1] "Tests q < .1: 1"
```

Caudate had more nominally significant tests in general than ACC where the gene
expression differed by diagnosis. However, a few more in ACC were significant
using FDR. Looking at the Excel spreadsheet, they have to be the 3 lowest
p-values:

![](images/2020-03-14-18-22-50.png)

The top 2 are the same one in both p thresholds.

Also, filtering at the less conservative p<.1 seemed to show more results.

Maybe a different approach here would be to take some sort ofintersection list
of nominally significant genes across the different model we ran, and so sort of
gene set analysis?

But for now let's look at some of the other interactions:

```r
res0 = read.csv('res_ACC_pLT0.05_DiagnosisAge.csv')
res = res0[res0$predictor=='DiagnosisControl',]
p = res[, 'Pr...t..']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with DX q < .1: %d', sum(p2<.1)))
p1 = p
res = res0[res0$predictor=='DiagnosisControl:Age',]
p = res[, 'Pr...t..']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with Age:DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with Age:DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with Age:DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with Age:DX q < .1: %d', sum(p2<.1)))
print(sprintf('Tests with both DX and Age:DX p < .05: %d', sum(p<.05 & p1<.05)))
print(sprintf('Tests with both DX and Age:DX p < .01: %d', sum(p<.01 & p1<.01)))
```

```
> res0 = read.csv('res_ACC_pLT0.05_DiagnosisAge.csv')
[1] "Tests with DX p < .05: 1867"
[1] "Tests with DX p < .01: 378"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Age:DX p < .05: 1861"
[1] "Tests with Age:DX p < .01: 353"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 0"
[1] "Tests with both DX and Age:DX p < .05: 1215"
[1] "Tests with both DX and Age:DX p < .01: 222"
> res0 = read.csv('res_ACC_pLT0.10_DiagnosisAge.csv')
[1] "Tests with DX p < .05: 2160"
[1] "Tests with DX p < .01: 451"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Age:DX p < .05: 2143"
[1] "Tests with Age:DX p < .01: 439"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 0"
[1] "Tests with both DX and Age:DX p < .05: 1393"
[1] "Tests with both DX and Age:DX p < .01: 261"
> res0 = read.csv('res_Caudate_pLT0.05_DiagnosisAge.csv')
[1] "Tests with DX p < .05: 2355"
[1] "Tests with DX p < .01: 498"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Age:DX p < .05: 2834"
[1] "Tests with Age:DX p < .01: 724"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 1"
[1] "Tests with both DX and Age:DX p < .05: 1833"
[1] "Tests with both DX and Age:DX p < .01: 381"
> res0 = read.csv('res_Caudate_pLT0.10_DiagnosisAge.csv')
[1] "Tests with DX p < .05: 2689"
[1] "Tests with DX p < .01: 611"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Age:DX p < .05: 3260"
[1] "Tests with Age:DX p < .01: 858"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 1"
[1] "Tests with both DX and Age:DX p < .05: 2099"
[1] "Tests with both DX and Age:DX p < .01: 465"
```

Like before, filtering covariates at p < .1 seemed to have more results than p <
.05. In general, there were more gene expression variables significant for
Age:DX than DX by itself. Again, the caudate seems to be a bit better in this
model, especially because this no FDR adjusted terms came out of the ACC
regressions. For reference, the single FDR result is:

![](images/2020-03-14-18-35-37.png)

The same in both p thresholds. 

Let's take a look at the LME models:

```r
res0 = read.csv('res_pLT0.05_DiagnosisRegion.csv')
res = res0[res0$predictor=='DiagnosisControl',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with DX q < .1: %d', sum(p2<.1)))
p1 = p
res = res0[res0$predictor=='DiagnosisControl:RegionCaudate',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with Region:DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with Region:DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with Region:DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with Region:DX q < .1: %d', sum(p2<.1)))
print(sprintf('Tests with both DX and Region:DX p < .05: %d', sum(p<.05 & p1<.05)))
print(sprintf('Tests with both DX and Region:DX p < .01: %d', sum(p<.01 & p1<.01)))
```

```
> res0 = read.csv('res_pLT0.05_DiagnosisRegion.csv')
[1] "Tests with DX p < .05: 1681"
[1] "Tests with DX p < .01: 373"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Region:DX p < .05: 1396"
[1] "Tests with Region:DX p < .01: 269"
[1] "Tests with Region:DX q < .05: 0"
[1] "Tests with Region:DX q < .1: 0"
[1] "Tests with both DX and Region:DX p < .05: 290"
[1] "Tests with both DX and Region:DX p < .01: 36"
> res0 = read.csv('res_pLT0.10_DiagnosisRegion.csv')
[1] "Tests with DX p < .05: 1933"
[1] "Tests with DX p < .01: 481"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Region:DX p < .05: 1425"
[1] "Tests with Region:DX p < .01: 260"
[1] "Tests with Region:DX q < .05: 0"
[1] "Tests with Region:DX q < .1: 0"
[1] "Tests with both DX and Region:DX p < .05: 292"
[1] "Tests with both DX and Region:DX p < .01: 37"
```

Same patterns we were seeing before seem to hold here. Nothing for FDR though,
but I'm not too worried as we might have other methods to analyze this.

```r
res0 = read.csv('res_pLT0.05_Diagnosis.csv')
res = res0[res0$predictor=='DiagnosisControl',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with DX q < .1: %d', sum(p2<.1)))
```

```
> res0 = read.csv('res_pLT0.05_Diagnosis.csv')
[1] "Tests with DX p < .05: 2365"
[1] "Tests with DX p < .01: 595"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
> res0 = read.csv('res_pLT0.10_Diagnosis.csv')
[1] "Tests with DX p < .05: 2719"
[1] "Tests with DX p < .01: 699"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
```

Those results are for the model that holds Region fixed as an additive
covariate. I didn't report how many tests had that as significant because, as
expected, there were many of them. Same pattern for DX significant holds. Also,
good to remember here that there are always 35917 tests.

```r
res0 = read.csv('res_pLT0.05_DiagnosisAge.csv')
res = res0[res0$predictor=='DiagnosisControl',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with DX q < .1: %d', sum(p2<.1)))
p1 = p
res = res0[res0$predictor=='DiagnosisControl:Age',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with Age:DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with Age:DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with Age:DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with Age:DX q < .1: %d', sum(p2<.1)))
print(sprintf('Tests with both DX and Age:DX p < .05: %d', sum(p<.05 & p1<.05)))
print(sprintf('Tests with both DX and Age:DX p < .01: %d', sum(p<.01 & p1<.01)))
```

```
> res0 = read.csv('res_pLT0.05_DiagnosisAge.csv')
[1] "Tests with DX p < .05: 1964"
[1] "Tests with DX p < .01: 350"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Age:DX p < .05: 2028"
[1] "Tests with Age:DX p < .01: 421"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 0"
[1] "Tests with both DX and Age:DX p < .05: 1311"
[1] "Tests with both DX and Age:DX p < .01: 229"
> res0 = read.csv('res_pLT0.10_DiagnosisAge.csv')
[1] "Tests with DX p < .05: 2082"
[1] "Tests with DX p < .01: 388"
[1] "Tests with DX q < .05: 0"
[1] "Tests with DX q < .1: 0"
[1] "Tests with Age:DX p < .05: 2180"
[1] "Tests with Age:DX p < .01: 453"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 0"
[1] "Tests with both DX and Age:DX p < .05: 1381"
[1] "Tests with both DX and Age:DX p < .01: 239"
```

Those results were also holding Region as a fixed additive covariate. Same usual
patterns we're seeing.

```r
res0 = read.csv('res_pLT0.05_DiagnosisAgeRegion.csv')
res = res0[res0$predictor=='DiagnosisControl',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with DX q < .1: %d', sum(p2<.1)))
p1 = p
res = res0[res0$predictor=='DiagnosisControl:Age',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with Age:DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with Age:DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with Age:DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with Age:DX q < .1: %d', sum(p2<.1)))
p3 = p
res = res0[res0$predictor=='DiagnosisControl:RegionCaudate',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with Region:DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with Region:DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with Region:DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with Region:DX q < .1: %d', sum(p2<.1)))
p4 = p
res = res0[res0$predictor=='DiagnosisControl:Age:RegionCaudate',]
p = res[, 'p.value']
p2 = p.adjust(p, method='fdr')
print(sprintf('Tests with Region:Age:DX p < .05: %d', sum(p<.05)))
print(sprintf('Tests with Region:Age:DX p < .01: %d', sum(p<.01)))
print(sprintf('Tests with Region:Age:DX q < .05: %d', sum(p2<.05)))
print(sprintf('Tests with Region:Age:DX q < .1: %d', sum(p2<.1)))
p5 = p
```

```
> res0 = read.csv('res_pLT0.05_DiagnosisAgeRegion.csv')
[1] "Tests with DX p < .05: 2130"
[1] "Tests with DX p < .01: 425"
[1] "Tests with DX q < .05: 1"
[1] "Tests with DX q < .1: 1"
[1] "Tests with Age:DX p < .05: 2068"
[1] "Tests with Age:DX p < .01: 408"
[1] "Tests with Age:DX q < .05: 1"
[1] "Tests with Age:DX q < .1: 1"
[1] "Tests with Region:DX p < .05: 2608"
[1] "Tests with Region:DX p < .01: 671"
[1] "Tests with Region:DX q < .05: 1"
[1] "Tests with Region:DX q < .1: 1"
[1] "Tests with Region:Age:DX p < .05: 2432"
[1] "Tests with Region:Age:DX p < .01: 566"
[1] "Tests with Region:Age:DX q < .05: 0"
[1] "Tests with Region:Age:DX q < .1: 0"
> res0 = read.csv('res_pLT0.10_DiagnosisAgeRegion.csv')
[1] "Tests with DX p < .05: 2167"
[1] "Tests with DX p < .01: 424"
[1] "Tests with DX q < .05: 1"
[1] "Tests with DX q < .1: 1"
[1] "Tests with Age:DX p < .05: 2131"
[1] "Tests with Age:DX p < .01: 431"
[1] "Tests with Age:DX q < .05: 0"
[1] "Tests with Age:DX q < .1: 0"
[1] "Tests with Region:DX p < .05: 2606"
[1] "Tests with Region:DX p < .01: 648"
[1] "Tests with Region:DX q < .05: 0"
[1] "Tests with Region:DX q < .1: 0"
[1] "Tests with Region:Age:DX p < .05: 2439"
[1] "Tests with Region:Age:DX p < .01: 561"
[1] "Tests with Region:Age:DX q < .05: 0"
[1] "Tests with Region:Age:DX q < .1: 0"
```

Just for kicks, I ran Meff to see where we'd be. I did it withinACC, Caudate,
and then using both:

```r
data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data = data[data$Region=='ACC', 34:ncol(data)
# some variables have zero sd after removing one of the regions
sds = apply(data, 2, sd)
keep_me = which(sds>0)
data = data[, keep_me]
cc = cor(data)
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

So, if I include both regions the divider is 77.70. For ACC it's 44.26, and for Caudate
it's 47.54. Even for the biggest one, we're looing at 0.00064, which is not too
bad. 

Philip also liked that the p<.1 results gave a few more significant gene
expression. So, let's stick with that for now, add in the population code (while
we decide whether or not to use the population PCs from Chandra), and re-run
everything. The other option would be to add the different domains sequentially,
in blocks, to see what works best. Add them to the fixed variables, I mean.

```r
library(nlme)
pthresh = .1
keep_str = 'Diagnosis*Age'

data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data$substance_group = as.factor(data$substance_group)
data$batch = as.factor(data$batch)
data$hbcc_brain_id = as.factor(data$hbcc_brain_id)
# no column names as numbers!
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names

pop_code = read.csv('~/data/rnaseq_derek/file_pop.csv')
data = merge(data, pop_code, by='hbcc_brain_id')

# dependent
dep_vars = colnames(data)[grepl(colnames(data), pattern='^grex')]
# keep these regardless of significance
keep_vars = c(keep_str, 'Region')
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'pH', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'POP_CODE')
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/resPOP_pLT%.02f_%s.csv', pthresh,
                    gsub(pattern='\\*',replacement='',x=keep_str))

hold = c()
for (dp in 1:length(dep_vars)) {
    if (dp %% 50 == 0) {
        print(sprintf('%d of %d (%s)', dp, length(dep_vars), out_fname))
    }

    dep_var = dep_vars[dp]
    fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                   paste(test_vars, collapse='+'), sep="")
    fit = try(lme(as.formula(fm_str), ~1|hbcc_brain_id, data=data, na.action=na.omit))
    if (length(fit) > 1) {
        res = summary(fit)$tTable
        # filtering variables
        sig_vars = c()
        for (v in 1:length(test_vars)) {
            # rows in results table that correspond to the screened variable
            var_rows = which(grepl(rownames(res),
                            pattern=sprintf('^%s', test_vars[v])))
            for (r in var_rows) {
                if (res[r, 'p-value'] < pthresh) {
                    sig_vars = c(sig_vars, test_vars[v])
                }
            }
        }
        # factors might get added several times, so here we clean it up
        sig_vars = unique(sig_vars)
        if (length(sig_vars) > 0) {
            clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                        paste(sig_vars, collapse='+'), sep="")
        } else {
            clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), sep="")
        }
        # new model
        clean_fit = try(lme(as.formula(clean_fm_str), ~1|hbcc_brain_id, data=data,
                        na.action=na.omit))
        if (length(clean_fit) > 1) {
            res = data.frame(summary(clean_fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = clean_fm_str
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        } else {
            res = data.frame(summary(fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = NA
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        }
        hold = rbind(hold, res)
    }
}
write.csv(hold, file=out_fname, row.names=F)
```

I'll also try running WNH only, just in case:

```r
library(nlme)
pthresh = .1
keep_str = 'Diagnosis*Region'

data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data$substance_group = as.factor(data$substance_group)
data$batch = as.factor(data$batch)
data$hbcc_brain_id = as.factor(data$hbcc_brain_id)
# no column names as numbers!
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names

pop_code = read.csv('~/data/rnaseq_derek/file_pop.csv')
data = merge(data, pop_code, by='hbcc_brain_id')
data = data[data$POP_CODE=='WNH', ]

# dependent
dep_vars = colnames(data)[grepl(colnames(data), pattern='^grex')]
# keep these regardless of significance
keep_vars = c(keep_str)
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'pH', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'Age')
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/resWNH_pLT%.02f_%s.csv', pthresh,
                    gsub(pattern='\\*',replacement='',x=keep_str))

hold = c()
for (dp in 1:length(dep_vars)) {
    if (dp %% 50 == 0) {
        print(sprintf('%d of %d (%s)', dp, length(dep_vars), out_fname))
    }

    dep_var = dep_vars[dp]
    fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                   paste(test_vars, collapse='+'), sep="")
    fit = try(lme(as.formula(fm_str), ~1|hbcc_brain_id, data=data, na.action=na.omit))
    if (length(fit) > 1) {
        res = summary(fit)$tTable
        # filtering variables
        sig_vars = c()
        for (v in 1:length(test_vars)) {
            # rows in results table that correspond to the screened variable
            var_rows = which(grepl(rownames(res),
                            pattern=sprintf('^%s', test_vars[v])))
            for (r in var_rows) {
                if (res[r, 'p-value'] < pthresh) {
                    sig_vars = c(sig_vars, test_vars[v])
                }
            }
        }
        # factors might get added several times, so here we clean it up
        sig_vars = unique(sig_vars)
        if (length(sig_vars) > 0) {
            clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                        paste(sig_vars, collapse='+'), sep="")
        } else {
            clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), sep="")
        }
        # new model
        clean_fit = try(lme(as.formula(clean_fm_str), ~1|hbcc_brain_id, data=data,
                        na.action=na.omit))
        if (length(clean_fit) > 1) {
            res = data.frame(summary(clean_fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = clean_fm_str
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        } else {
            res = data.frame(summary(fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = NA
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        }
        hold = rbind(hold, res)
    }
}
write.csv(hold, file=out_fname, row.names=F)
```

And we can start playing with the lm models as well:

```r
myregion = 'ACC'
pthresh = .1
keep_str = 'Diagnosis'

data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
data$substance_group = as.factor(data$substance_group)
data$batch = as.factor(data$batch)
# no column names as numbers!
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names

pop_code = read.csv('~/data/rnaseq_derek/file_pop.csv')
data = merge(data, pop_code, by='hbcc_brain_id')
# data = data[data$POP_CODE=='WNH', ]

data = data[data$Region==myregion, ]

# dependent
dep_vars = colnames(data)[grepl(colnames(data), pattern='^grex')]
# keep these regardless of significance
keep_vars = c(keep_str)
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'pH', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'POP_CODE', 'Age')
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/resPOP_%s_pLT%.02f_%s.csv', myregion, pthresh,
                    gsub(pattern='\\*',replacement='',x=keep_str))

hold = c()
for (dp in 1:length(dep_vars)) {
    if (dp %% 50 == 0) {
        print(sprintf('%d of %d (%s)', dp, length(dep_vars), out_fname))
    }

    dep_var = dep_vars[dp]
    fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                   paste(test_vars, collapse='+'), sep="")
    fit = lm(as.formula(fm_str), data=data)
    res = summary(fit)$coefficients
    # filtering variables
    sig_vars = c()
    for (v in 1:length(test_vars)) {
        # rows in results table that correspond to the screened variable
        var_rows = which(grepl(rownames(res),
                         pattern=sprintf('^%s', test_vars[v])))
        for (r in var_rows) {
            if (res[r, 'Pr(>|t|)'] < pthresh) {
                sig_vars = c(sig_vars, test_vars[v])
            }
        }
    }
    # factors might get added several times, so here we clean it up
    sig_vars = unique(sig_vars)
    if (length(sig_vars) > 0) {
        clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), ' + ',
                       paste(sig_vars, collapse='+'), sep="")
    } else {
        clean_fm_str = paste(dep_var, ' ~ ', paste(keep_vars, collapse='+'), sep="")
    }
    # new model
    clean_fit = lm(as.formula(clean_fm_str), data=data)
    res = data.frame(summary(clean_fit)$coefficients)
    # remove intercept
    res = res[2:nrow(res),]
    res$dep_var = dep_var
    res$formula = clean_fm_str
    res$orig_formula = fm_str
    res$predictor = rownames(res)
    hold = rbind(hold, res)
}
write.csv(hold, file=out_fname, row.names=F)
```

# 2020-03-18 19:38:05

Let's gather the results again. I had to create a script because it was getting
too cumbersome to cut and paste code. So, now we run:

```bash
Rscript ~/research_code/compile_rnaseq_results.R resWNH_ACC_pLT0.10_DiagnosisAge.csv
```

```
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_ACC_pLT0.10_DiagnosisAge.csv
[1] "resPOP_ACC_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2172"
[1] "Tests with DiagnosisControl p < .01: 459"
[1] "Tests with DiagnosisControl p < Meff: 64"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 2"
[1] "Tests with DiagnosisControl:Age p < .05: 2219"
[1] "Tests with DiagnosisControl:Age p < .01: 473"
[1] "Tests with DiagnosisControl:Age p < Meff: 66"
[1] "Tests with DiagnosisControl:Age q < .05: 1"
[1] "Tests with DiagnosisControl:Age q < .1: 2"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.05: 1437"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.01: 280"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_Caudate_pLT0.10_DiagnosisAge.csv
[1] "resPOP_Caudate_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2709"
[1] "Tests with DiagnosisControl p < .01: 609"
[1] "Tests with DiagnosisControl p < Meff: 89"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 3309"
[1] "Tests with DiagnosisControl:Age p < .01: 895"
[1] "Tests with DiagnosisControl:Age p < Meff: 116"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 0"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.05: 2119"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.01: 469"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_ACC_pLT0.10_Diagnosis.csv
[1] "resPOP_ACC_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 2626"
[1] "Tests with DiagnosisControl p < .01: 726"
[1] "Tests with DiagnosisControl p < Meff: 132"
[1] "Tests with DiagnosisControl q < .05: 2"
[1] "Tests with DiagnosisControl q < .1: 3"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_Caudate_pLT0.10_Diagnosis.csv
[1] "resPOP_Caudate_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 2966"
[1] "Tests with DiagnosisControl p < .01: 795"
[1] "Tests with DiagnosisControl p < Meff: 125"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_ACC_pLT0.10_DiagnosisAge.csv
[1] "resWNH_ACC_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 4181"
[1] "Tests with DiagnosisControl p < .01: 1231"
[1] "Tests with DiagnosisControl p < Meff: 258"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 74"
[1] "Tests with DiagnosisControl:Age p < .05: 4198"
[1] "Tests with DiagnosisControl:Age p < .01: 1354"
[1] "Tests with DiagnosisControl:Age p < Meff: 299"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 110"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.05: 3207"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.01: 948"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_Caudate_pLT0.10_DiagnosisAge.csv
[1] "resWNH_Caudate_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2602"
[1] "Tests with DiagnosisControl p < .01: 612"
[1] "Tests with DiagnosisControl p < Meff: 85"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 2168"
[1] "Tests with DiagnosisControl:Age p < .01: 498"
[1] "Tests with DiagnosisControl:Age p < Meff: 90"
[1] "Tests with DiagnosisControl:Age q < .05: 1"
[1] "Tests with DiagnosisControl:Age q < .1: 1"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.05: 1586"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.01: 364"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_ACC_pLT0.10_Diagnosis.csv
[1] "resWNH_ACC_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3739"
[1] "Tests with DiagnosisControl p < .01: 1492"
[1] "Tests with DiagnosisControl p < Meff: 415"
[1] "Tests with DiagnosisControl q < .05: 110"
[1] "Tests with DiagnosisControl q < .1: 442"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_Caudate_pLT0.10_Diagnosis.csv
[1] "resWNH_Caudate_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3357"
[1] "Tests with DiagnosisControl p < .01: 879"
[1] "Tests with DiagnosisControl p < Meff: 144"
[1] "Tests with DiagnosisControl q < .05: 5"
[1] "Tests with DiagnosisControl q < .1: 5"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_pLT0.10_Diagnosis.csv
[1] "resPOP_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 2756"
[1] "Tests with DiagnosisControl p < .01: 740"
[1] "Tests with DiagnosisControl p < Meff: 72"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_pLT0.10_Diagnosis.csv
[1] "resWNH_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 2901"
[1] "Tests with DiagnosisControl p < .01: 698"
[1] "Tests with DiagnosisControl p < Meff: 55"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_pLT0.10_DiagnosisAge.csv
[1] "resPOP_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2132"
[1] "Tests with DiagnosisControl p < .01: 406"
[1] "Tests with DiagnosisControl p < Meff: 22"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 2201"
[1] "Tests with DiagnosisControl:Age p < .01: 449"
[1] "Tests with DiagnosisControl:Age p < Meff: 30"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 0"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.05: 1414"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.01: 243"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_pLT0.10_DiagnosisAge.csv
[1] "resWNH_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2965"
[1] "Tests with DiagnosisControl p < .01: 579"
[1] "Tests with DiagnosisControl p < Meff: 35"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 2198"
[1] "Tests with DiagnosisControl:Age p < .01: 389"
[1] "Tests with DiagnosisControl:Age p < Meff: 28"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 0"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.05: 1732"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age p < 0.01: 276"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_pLT0.10_DiagnosisRegion.csv
[1] "resPOP_pLT0.10_DiagnosisRegion.csv"
[1] "Tests with DiagnosisControl p < .05: 2014"
[1] "Tests with DiagnosisControl p < .01: 503"
[1] "Tests with DiagnosisControl p < Meff: 45"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 1427"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 260"
[1] "Tests with DiagnosisControl:RegionCaudate p < Meff: 16"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 0"
[1] "Tests with DiagnosisControl and DiagnosisControl:RegionCaudate p < 0.05: 276"
[1] "Tests with DiagnosisControl and DiagnosisControl:RegionCaudate p < 0.01: 33"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_pLT0.10_DiagnosisRegion.csv
[1] "resWNH_pLT0.10_DiagnosisRegion.csv"
[1] "Tests with DiagnosisControl p < .05: 2217"
[1] "Tests with DiagnosisControl p < .01: 457"
[1] "Tests with DiagnosisControl p < Meff: 32"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 1243"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 232"
[1] "Tests with DiagnosisControl:RegionCaudate p < Meff: 18"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 0"
[1] "Tests with DiagnosisControl and DiagnosisControl:RegionCaudate p < 0.05: 287"
[1] "Tests with DiagnosisControl and DiagnosisControl:RegionCaudate p < 0.01: 51"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOP_pLT0.10_DiagnosisRegionAge.csv
[1] "resPOP_pLT0.10_DiagnosisRegionAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2143"
[1] "Tests with DiagnosisControl p < .01: 418"
[1] "Tests with DiagnosisControl p < Meff: 22"
[1] "Tests with DiagnosisControl q < .05: 1"
[1] "Tests with DiagnosisControl q < .1: 1"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 2564"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 638"
[1] "Tests with DiagnosisControl:RegionCaudate p < Meff: 59"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 1"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 1"
[1] "Tests with DiagnosisControl:Age p < .05: 2148"
[1] "Tests with DiagnosisControl:Age p < .01: 421"
[1] "Tests with DiagnosisControl:Age p < Meff: 22"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < .05: 2407"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < .01: 546"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < Meff: 43"
[1] "Tests with DiagnosisControl:RegionCaudate:Age q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate:Age q < .1: 0"
[1] "Tests with DiagnosisControl and DiagnosisControl:RegionCaudate and DiagnosisControl:Age and DiagnosisControl:RegionCaudate:Age p < 0.05: 391"
[1] "Tests with DiagnosisControl and DiagnosisControl:RegionCaudate and DiagnosisControl:Age and DiagnosisControl:RegionCaudate:Age p < 0.01: 43"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNH_pLT0.10_DiagnosisAgeRegion.csv
[1] "resWNH_pLT0.10_DiagnosisAgeRegion.csv"
[1] "Tests with DiagnosisControl p < .05: 3659"
[1] "Tests with DiagnosisControl p < .01: 699"
[1] "Tests with DiagnosisControl p < Meff: 39"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 2757"
[1] "Tests with DiagnosisControl:Age p < .01: 483"
[1] "Tests with DiagnosisControl:Age p < Meff: 35"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 3150"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 821"
[1] "Tests with DiagnosisControl:RegionCaudate p < Meff: 59"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 4"
[1] "Tests with DiagnosisControl:Age:RegionCaudate p < .05: 2618"
[1] "Tests with DiagnosisControl:Age:RegionCaudate p < .01: 584"
[1] "Tests with DiagnosisControl:Age:RegionCaudate p < Meff: 38"
[1] "Tests with DiagnosisControl:Age:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:Age:RegionCaudate q < .1: 2"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age and DiagnosisControl:RegionCaudate and DiagnosisControl:Age:RegionCaudate p < 0.05: 541"
[1] "Tests with DiagnosisControl and DiagnosisControl:Age and DiagnosisControl:RegionCaudate and DiagnosisControl:Age:RegionCaudate p < 0.01: 62"
```

* play with adding the different covariate domains sequentially
* can Meff be used here safely?