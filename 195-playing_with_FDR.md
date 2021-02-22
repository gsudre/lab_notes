# 2021-02-19 11:38:12

In working with the methylation data I ran into a few methods to increase the
power in FDR:

 * IHW
 * CAMT
 * adaptMT

Let's see if any of them make a difference in our results.

## CAMT

I'll start with this one just because of the function that tests possible
covariates.

```r
library(CAMT)
library(minfi)
load('~/data/methylation_post_mortem/filt_ACC_02182021.RData')
load('~/data/methylation_post_mortem/res_ACC_02182021.RData')
#vals <- getM(mSetSqFlt)
vals <- getBeta(mSetSqFlt)

bad_probes = rownames(which(abs(vals)==Inf, arr.ind = T))
mSetSqFlt = mSetSqFlt[!(rownames(vals) %in% bad_probes), ]
vals = vals[!(rownames(vals) %in% bad_probes), ]
pvals = res_acc[['all']]$DMPs[, 'P.Value']
names(pvals) = rownames(res_acc[['all']]$DMPs)

covars = rowMeans(vals)
# covars = apply(vals, 1, sd)

# sorting so the variables correspond to each other
covars2 = covars[match(names(pvals), names(covars))]
set.seed(42)
stats = cov.test.n(pvals, covars2, silence=F, perm.no = 1000)
```

I also tested meanB, and sd for M and B. Here are the results:

```
# beta, sd
$p.value
[1] 0.000999001

$x.cut.optim
[1] 8

$p.cut.optim
      10% 
0.1050184 

# beta, mean
$p.value
[1] 0.000999001

$x.cut.optim
[1] 16

$p.cut.optim
      10% 
0.1050184 

# M, sd
$p.value
[1] 0.000999001

$x.cut.optim
[1] 16

$p.cut.optim
     20% 
0.204451 

# M, mean
$p.value
[1] 0.000999001

$x.cut.optim
[1] 16

$p.cut.optim
      10% 
0.1050184 
```

These values are pretty much the same. The test statistic is different for all
of them (not copied), so there's no issues with the procedure itself. I need to
read the method a bit more to see how to use it, though. It's not as simple as
just using one covariate.

## IHW

I've been using this for a while now:

```r
library("IHW")
ihwRes <- ihw(pvals ~ covars2,  alpha = 0.05)
padjBH <- p.adjust(pvals, method = "BH")
cat('IHW:', sum(adj_pvalues(ihwRes) <= 0.05, na.rm = TRUE))
cat('FDR:', sum(padjBH <= 0.05, na.rm = TRUE))
```

IHW didn't give any new rejections at .05 for either of the 4 scenarios.

## adaptMT

Here I'll try not only the single variable model, but maybe a combination of
them too.

```r
library("adaptMT")
library("splines")
x <- data.frame(x = covars2)
formulas <- paste0("ns(x, df = ", 6:10, ")")
res_glm <- adapt_glm(x = x, pvals = pvals, pi_formulas = formulas,
                    mu_formulas = formulas)
```

```
# SD M
alpha = 0.13: FDPhat 0.1111, Number of Rej. 9
alpha = 0.12: FDPhat 0.1111, Number of Rej. 9

# mean M
alpha = 0.6: FDPhat 0.6, Number of Rej. 5

# mean B
alpha = 0.38: FDPhat 0.375, Number of Rej. 8
alpha = 0.37: FDPhat 0.3333, Number of Rej. 3
alpha = 0.36: FDPhat 0.3333, Number of Rej. 3
alpha = 0.35: FDPhat 0.3333, Number of Rej. 3
alpha = 0.34: FDPhat 0.3333, Number of Rej. 3

# SD B
alpha = 0.36: FDPhat 0.3571, Number of Rej. 14
alpha = 0.35: FDPhat 0.3077, Number of Rej. 13
alpha = 0.34: FDPhat 0.3077, Number of Rej. 13
alpha = 0.33: FDPhat 0.3077, Number of Rej. 13
alpha = 0.32: FDPhat 0.3077, Number of Rej. 13
alpha = 0.31: FDPhat 0.3077, Number of Rej. 13
```

Not there yet, but let's try combining the covariates.

```r
load('~/data/methylation_post_mortem/filt_ACC_02182021.RData')
load('~/data/methylation_post_mortem/res_ACC_02182021.RData')
vals <- getM(mSetSqFlt)
#vals <- getBeta(mSetSqFlt)

bad_probes = rownames(which(abs(vals)==Inf, arr.ind = T))
mSetSqFlt = mSetSqFlt[!(rownames(vals) %in% bad_probes), ]
vals = vals[!(rownames(vals) %in% bad_probes), ]
pvals = res_acc[['all']]$DMPs[, 'P.Value']
names(pvals) = rownames(res_acc[['all']]$DMPs)

mus = rowMeans(vals)
sds = apply(vals, 1, sd)

# sorting so the variables correspond to each other
mus2 = mus[match(names(pvals), names(mus))]
sds2 = sds[match(names(pvals), names(sds))]

library("mgcv")
set.seed(42)
x <- data.frame(x1 = mus2, x2=sds2)
formula <- "s(x1, x2)"
res_gam <- adapt_gam(x = x, pvals = pvals, pi_formulas = formula,
                     mu_formulas = formula)
```

And if we want to use both means or both sds, we do:

```r
load('~/data/methylation_post_mortem/filt_ACC_02182021.RData')
load('~/data/methylation_post_mortem/res_ACC_02182021.RData')
library(minfi)
bvals <- getM(mSetSqFlt)
mvals <- getBeta(mSetSqFlt)

bad_probes = rownames(which(abs(bvals)==Inf, arr.ind = T))
mSetSqFlt = mSetSqFlt[!(rownames(bvals) %in% bad_probes), ]
bvals = bvals[!(rownames(bvals) %in% bad_probes), ]
mvals = mvals[!(rownames(mvals) %in% bad_probes), ]
pvals = res_acc[['all']]$DMPs[, 'P.Value']
names(pvals) = rownames(res_acc[['all']]$DMPs)

m = rowMeans(mvals)
b = rowMeans(bvals)

# m = apply(mvals, 1, sd)
# b = apply(bvals, 1, sd)

# sorting so the variables correspond to each other
m2 = m[match(names(pvals), names(m))]
b2 = b[match(names(pvals), names(b))]

library("mgcv")
library("adaptMT")
library("splines")
set.seed(42)
x <- data.frame(x1 = m2, x2=b2)
formula <- "s(x1, x2)"
res_gam <- adapt_gam(x = x, pvals = pvals, pi_formulas = formula,
                     mu_formulas = formula)
```

And we can try a version with all 4 variables:

```r
mm = rowMeans(mvals)
mb = rowMeans(bvals)
sm = apply(mvals, 1, sd)
sb = apply(bvals, 1, sd)

# sorting so the variables correspond to each other
mm2 = mm[match(names(pvals), names(mm))]
mb2 = mb[match(names(pvals), names(mb))]
sm2 = sm[match(names(pvals), names(sm))]
sb2 = sb[match(names(pvals), names(sb))]

library("mgcv")
set.seed(42)
x <- data.frame(x1 = mm2, x2=mb2, x3 = sm2, x4=sb2)
formula <- "s(x1, x2, x3, x4)"
res_gam <- adapt_gam(x = x, pvals = pvals, pi_formulas = formula,
                     mu_formulas = formula)
```

Finally, let's see what we can get with glmnet:

```r
res <- adapt_glmnet(as.matrix(x), pvals)
```