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
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names
data = data[data$Region=='ACC', grex_names]
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
data = data[data$POP_CODE=='WNH', ]

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
              'Sex', 'Age')  # POP_CODE
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

Then we can look for some intersections in the lists:

```bash
sort file1 file2 | uniq -d > out8
```

For example, the intersection of the main result and WNH only:

```
(base) HG-02035307-LM3:rnaseq_derek sudregp$ sort grexlist_resPOP_pLT0.10_DiagnosisRegionAge_DiagnosisControl:RegionCaudate:Age_Meff.txt grexlist_resWNH_pLT0.10_DiagnosisAgeRegion_DiagnosisControl:Age:RegionCaudate_Meff.txt | uniq -d
grex13678
grex26766
grex4894
grex9272
(base) HG-02035307-LM3:rnaseq_derek sudregp$ sort grexlist_resPOP_pLT0.10_Diagnosis_DiagnosisControl_Meff.txt grexlist_resWNH_pLT0.10_Diagnosis_DiagnosisControl_Meff.txt | uniq -d
grex13047
grex1325
grex21764
grex2308
grex3593
grex6827
```

# 2020-03-20 20:56:08

I figured out what happened to the RNAseq for WNH ACC. When we subset it to WNH
only, we are left with 35 subjects in ACC (37 for caudate). Depending on how
many variables are left after screening, the model completely overfits the data
(R2=.9999), because the pH variable has lots of NAs. Specifically, 23 in
subjects in ACC (18 WNH) don't have the pH data. That only leaves 17 subjects to
run the model with. That seems to be the only variable with so many NAs, so I'll
try the analysis without adding it as a covariate. We could always test the good
hits later to see if pH makes any difference to it.

So, let's re-run the code without pH.

```r
library(nlme)
pthresh = .1
keep_str = 'Diagnosis'

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
              "bainbank", 'PMI', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'Age', 'POP_CODE')
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/resPOPnoPH_pLT%.02f_%s.csv', pthresh,
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
keep_str = 'Diagnosis'

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
keep_vars = c(keep_str, 'Region')
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'Age')
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/resWNHnoPH_pLT%.02f_%s.csv', pthresh,
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
myregion = 'Caudate'
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
data = data[data$POP_CODE=='WNH', ]
data = data[data$Region==myregion, ]
# dependent
dep_vars = colnames(data)[grepl(colnames(data), pattern='^grex')]
# keep these regardless of significance
keep_vars = c(keep_str)
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank", 'PMI', 'Manner.of.Death',
              # technical
              'batch', 'RINe',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex', 'Age')#, 'POP_CODE')
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/resWNHnoPH_%s_pLT%.02f_%s.csv',
                    myregion, pthresh,
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

# 2020-03-21 08:44:03

I'm trying to calculate Meff in Python instead, since it's not quite working in
R. It's giving an awnfully low number there.

```python
# bw
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
readRDS = robjects.r['readRDS']
df = readRDS('/home/sudregp/data/rnaseq_derek/data_from_philip.rds')
a = df.iloc[:, 33:len(df.columns)].corr()
import numpy as np
ev = np.linalg.eigvals(a)
```

Can Meff be used here safely? The actual article for Meff is here:
  https://onlinelibrary.wiley.com/doi/pdf/10.1002/gepi.20310 and its
  approximation (Keff) is here:
  https://onlinelibrary.wiley.com/doi/pdf/10.1002/gepi.20331
  
They have an even simpler formula for Meff there:

![](images/2020-03-21-15-04-39.png)

So, let's use that, because then we avoid calculating the SVD, which is quite
expensive:
  
```r
data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
grex_names = sapply(colnames(data)[34:ncol(data)],
                    function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names
data = data[data$Region=='ACC', grex_names]
# some variables have zero sd after removing one of the regions
sds = apply(data, 2, sd)
keep_me = which(sds>0)
data = data[, keep_me]
cc = cor(data)
M = nrow(cc)
cnt = 0
for (j in 1:M) {
    print(j)
    for (k in 1:M) {
        cnt = cnt + (1 - cc[j, k]**2)
    }
}
meff = 1 + cnt / M
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

OK, now our Meff is around 34K for all 3 cases, so it's definitely more
realistic. Also, it doesn't help anything.

But I do have to re-run the summary scripts without the pH variable, so here it
is:

```
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_ACC_pLT0.10_Diagnosis.csv 
[1] "resPOPnoPH_ACC_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3014"
[1] "Tests with DiagnosisControl p < .01: 829"
[1] "Tests with DiagnosisControl q < .05: 3"
[1] "Tests with DiagnosisControl q < .1: 3"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_ACC_pLT0.10_Diagnosis.csv 
[1] "resWNHnoPH_ACC_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3369"
[1] "Tests with DiagnosisControl p < .01: 905"
[1] "Tests with DiagnosisControl q < .05: 1"
[1] "Tests with DiagnosisControl q < .1: 1"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_Caudate_pLT0.10_Diagnosis.csv 
[1] "resPOPnoPH_Caudate_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3600"
[1] "Tests with DiagnosisControl p < .01: 974"
[1] "Tests with DiagnosisControl q < .05: 1"
[1] "Tests with DiagnosisControl q < .1: 1"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_Caudate_pLT0.10_Diagnosis.csv 
[1] "resWNHnoPH_Caudate_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3319"
[1] "Tests with DiagnosisControl p < .01: 874"
[1] "Tests with DiagnosisControl q < .05: 5"
[1] "Tests with DiagnosisControl q < .1: 10"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_ACC_pLT0.10_DiagnosisAge.csv 
[1] "resPOPnoPH_ACC_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2492"
[1] "Tests with DiagnosisControl p < .01: 526"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 1"
[1] "Tests with DiagnosisControl:Age p < .05: 2714"
[1] "Tests with DiagnosisControl:Age p < .01: 661"
[1] "Tests with DiagnosisControl:Age q < .05: 1"
[1] "Tests with DiagnosisControl:Age q < .1: 1"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_ACC_pLT0.10_DiagnosisAge.csv 
[1] "resWNHnoPH_ACC_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 3349"
[1] "Tests with DiagnosisControl p < .01: 709"
[1] "Tests with DiagnosisControl q < .05: 1"
[1] "Tests with DiagnosisControl q < .1: 2"
[1] "Tests with DiagnosisControl:Age p < .05: 3345"
[1] "Tests with DiagnosisControl:Age p < .01: 838"
[1] "Tests with DiagnosisControl:Age q < .05: 2"
[1] "Tests with DiagnosisControl:Age q < .1: 2"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_Caudate_pLT0.10_DiagnosisAge.csv 
[1] "resPOPnoPH_Caudate_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2642"
[1] "Tests with DiagnosisControl p < .01: 575"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 3270"
[1] "Tests with DiagnosisControl:Age p < .01: 798"
[1] "Tests with DiagnosisControl:Age q < .05: 1"
[1] "Tests with DiagnosisControl:Age q < .1: 1"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_Caudate_pLT0.10_DiagnosisAge.csv 
[1] "resWNHnoPH_Caudate_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2811"
[1] "Tests with DiagnosisControl p < .01: 622"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 2"
[1] "Tests with DiagnosisControl:Age p < .05: 2596"
[1] "Tests with DiagnosisControl:Age p < .01: 586"
[1] "Tests with DiagnosisControl:Age q < .05: 2"
[1] "Tests with DiagnosisControl:Age q < .1: 3"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_pLT0.10_Diagnosis.csv 
[1] "resPOPnoPH_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 3092"
[1] "Tests with DiagnosisControl p < .01: 755"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_pLT0.10_Diagnosis.csv 
[1] "resWNHnoPH_pLT0.10_Diagnosis.csv"
[1] "Tests with DiagnosisControl p < .05: 2924"
[1] "Tests with DiagnosisControl p < .01: 667"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_pLT0.10_DiagnosisAge.csv 
[1] "resPOPnoPH_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 1988"
[1] "Tests with DiagnosisControl p < .01: 393"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:Age p < .05: 2387"
[1] "Tests with DiagnosisControl:Age p < .01: 474"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 1"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_pLT0.10_DiagnosisAge.csv 
[1] "resWNHnoPH_pLT0.10_DiagnosisAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2865"
[1] "Tests with DiagnosisControl p < .01: 595"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 1"
[1] "Tests with DiagnosisControl:Age p < .05: 2629"
[1] "Tests with DiagnosisControl:Age p < .01: 597"
[1] "Tests with DiagnosisControl:Age q < .05: 1"
[1] "Tests with DiagnosisControl:Age q < .1: 1"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_pLT0.10_DiagnosisRegion.csv 
[1] "resPOPnoPH_pLT0.10_DiagnosisRegion.csv"
[1] "Tests with DiagnosisControl p < .05: 2139"
[1] "Tests with DiagnosisControl p < .01: 473"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 1578"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 296"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_pLT0.10_DiagnosisRegion.csv 
[1] "resWNHnoPH_pLT0.10_DiagnosisRegion.csv"
[1] "Tests with DiagnosisControl p < .05: 2219"
[1] "Tests with DiagnosisControl p < .01: 428"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 1438"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 270"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resPOPnoPH_pLT0.10_DiagnosisRegionAge.csv 
[1] "resPOPnoPH_pLT0.10_DiagnosisRegionAge.csv"
[1] "Tests with DiagnosisControl p < .05: 2167"
[1] "Tests with DiagnosisControl p < .01: 415"
[1] "Tests with DiagnosisControl q < .05: 1"
[1] "Tests with DiagnosisControl q < .1: 1"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 2961"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 778"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 1"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 1"
[1] "Tests with DiagnosisControl:Age p < .05: 2216"
[1] "Tests with DiagnosisControl:Age p < .01: 426"
[1] "Tests with DiagnosisControl:Age q < .05: 0"
[1] "Tests with DiagnosisControl:Age q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < .05: 2589"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < .01: 623"
[1] "Tests with DiagnosisControl:RegionCaudate:Age q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate:Age q < .1: 0"
(base) HG-02035307-LM3:rnaseq_derek sudregp$ Rscript ~/research_code/compile_rnaseq_results.R resWNHnoPH_pLT0.10_DiagnosisRegionAge.csv 
[1] "resWNHnoPH_pLT0.10_DiagnosisRegionAge.csv"
[1] "Tests with DiagnosisControl p < .05: 3389"
[1] "Tests with DiagnosisControl p < .01: 600"
[1] "Tests with DiagnosisControl q < .05: 0"
[1] "Tests with DiagnosisControl q < .1: 0"
[1] "Tests with DiagnosisControl:RegionCaudate p < .05: 3463"
[1] "Tests with DiagnosisControl:RegionCaudate p < .01: 891"
[1] "Tests with DiagnosisControl:RegionCaudate q < .05: 0"
[1] "Tests with DiagnosisControl:RegionCaudate q < .1: 5"
[1] "Tests with DiagnosisControl:Age p < .05: 2891"
[1] "Tests with DiagnosisControl:Age p < .01: 572"
[1] "Tests with DiagnosisControl:Age q < .05: 1"
[1] "Tests with DiagnosisControl:Age q < .1: 2"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < .05: 2743"
[1] "Tests with DiagnosisControl:RegionCaudate:Age p < .01: 605"
[1] "Tests with DiagnosisControl:RegionCaudate:Age q < .05: 2"
[1] "Tests with DiagnosisControl:RegionCaudate:Age q < .1: 2"
```

Let's do some filtering now. If we look at resPOPnoPH_ACC_pLT0.10_Diagnosis.csv
we have:

grex28267
grex29351
grex9100

as interesting candidates. In the WNH population, our top candidate is
grex20601, but those 3 genes are not doing poorly. 

0.00401195	grex28267
9.16E-05	grex29351
0.002668081	grex9100

For reference, grex20601 in the entire sample has p-value of 0.000369942. Now,
let's look at Caudate. The entire population only has one interesting candidate:
grex7973. But, if we look at the WNH population, we have 10 at q<.1:

grex5426
grex21133
grex18580
grex29858
grex29704
grex22959
grex31626
grex3223
grex14858
grex29583

For reference:

0.002067775	grex7973

## Population PCS

Let's calculate some of the population PCs. If anything, it'll be a good way to
check if the self-described population is working out fine.

* play with adding the different covariate domains sequentially
* maybe add PRS?
* add population covariates?

# 2020-04-14 12:38:06

Philip asked for a summary sheet:

```
could you send me a summary sheet that has these results all aligned by gene (with the p and q values).  And could you add in (2) the p value (and beta) of the region (ACC/caudate)* diagnosis interaction term- WNH and ALL; (3) the p value (and beta) for the region*diagnosis*age term - WNH and ALL; (4) the pvalue (and beta) for region*age; (5) the pvalues for ACC*age; (6) the pvalue for caudate*age- WNH/ALL.   In the file could you also include a column that indicates which GREX have the unusual distribution?
```

So, in R:

```r
a = read.csv('~/data/rnaseq_derek/resPOPnoPH_ACC_pLT0.10_Diagnosis.csv')
a = a[a$predictor=='DiagnosisControl',]
a$q = p.adjust(a[,'Pr...t..'], method='fdr')
a$predictor=NULL;
colnames(a) = sapply(colnames(a), function(x) sprintf('%s.%s', x, 'ACC_all_DX'))
b = read.csv('~/data/rnaseq_derek/resWNHnoPH_ACC_pLT0.10_Diagnosis.csv')
b = b[b$predictor=='DiagnosisControl',]
b$q = p.adjust(b[,'Pr...t..'], method='fdr')
b$predictor=NULL;
colnames(b) = sapply(colnames(b), function(x) sprintf('%s.%s', x, 'ACC_WNH_DX'))
m = merge(a, b, by=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_Caudate_pLT0.10_Diagnosis.csv')
d = d[d$predictor=='DiagnosisControl',]
d$q = p.adjust(d[,'Pr...t..'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'Caudate_all_DX'))
m = merge(m, d, by.x=1, by.y=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_Caudate_pLT0.10_Diagnosis.csv')
d = d[d$predictor=='DiagnosisControl',]
d$q = p.adjust(d[,'Pr...t..'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'Caudate_WNH_DX'))
m = merge(m, d, by.x=1, by.y=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_ACC_pLT0.10_DiagnosisAge.csv')
d = d[d$predictor=='DiagnosisControl:Age',]
d$q = p.adjust(d[,'Pr...t..'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'ACC_all_DX:Age'))
m = merge(m, d, by.x=1, by.y=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_ACC_pLT0.10_DiagnosisAge.csv')
d = d[d$predictor=='DiagnosisControl:Age',]
d$q = p.adjust(d[,'Pr...t..'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'ACC_WNH_DX:Age'))
m = merge(m, d, by.x=1, by.y=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_Caudate_pLT0.10_DiagnosisAge.csv')
d = d[d$predictor=='DiagnosisControl:Age',]
d$q = p.adjust(d[,'Pr...t..'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'Caudate_all_DX:Age'))
m = merge(m, d, by.x=1, by.y=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_Caudate_pLT0.10_DiagnosisAge.csv')
d = d[d$predictor=='DiagnosisControl:Age',]
d$q = p.adjust(d[,'Pr...t..'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'Caudate_WNH_DX:Age'))
m = merge(m, d, by.x=1, by.y=5, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_pLT0.10_Diagnosis.csv')
d = d[d$predictor=='DiagnosisControl',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'all_DX'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_pLT0.10_Diagnosis.csv')
d = d[d$predictor=='DiagnosisControl',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'WNH_DX'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_pLT0.10_DiagnosisAge.csv')
d = d[d$predictor=='DiagnosisControl:Age',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'all_DX:Age'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_pLT0.10_DiagnosisAge.csv')
d = d[d$predictor=='DiagnosisControl:Age',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'WNH_DX:Age'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_pLT0.10_DiagnosisRegion.csv')
d = d[d$predictor=='DiagnosisControl:RegionCaudate',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'all_DX:Region'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_pLT0.10_DiagnosisRegion.csv')
d = d[d$predictor=='DiagnosisControl:RegionCaudate',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'WNH_DX:Region'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resPOPnoPH_pLT0.10_DiagnosisRegionAge.csv')
d = d[d$predictor=='DiagnosisControl:RegionCaudate:Age',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'all_DX:Region:Age'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

d = read.csv('~/data/rnaseq_derek/resWNHnoPH_pLT0.10_DiagnosisRegionAge.csv')
d = d[d$predictor=='DiagnosisControl:RegionCaudate:Age',]
d$q = p.adjust(d[,'p.value'], method='fdr')
d$predictor=NULL;
colnames(d) = sapply(colnames(d), function(x) sprintf('%s.%s', x, 'WNH_DX:Region:Age'))
m = merge(m, d, by.x=1, by.y=6, all.x=T, all.y=F)

# just mark which variables are weird
m$imweird = F
data = readRDS('~/data/rnaseq_derek/data_from_philip.rds')
grex_names = sapply(colnames(data)[34:ncol(data)], function(x) sprintf('grex%s', x))
colnames(data)[34:ncol(data)] = grex_names
library(caret)
pp = preProcess(data[, grex_names], method=c('range'), rangeBounds=c(0,1))
a = predict(pp, data[, grex_names])
n0 = colSums(a==0)
imbad = names(n0)[n0>1]
m[m[, 1] %in% imbad, 'imweird'] = T

write.csv(m, file='~/data/rnaseq_derek/compiled_univariate.csv', row.names=F)
```
