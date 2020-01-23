# 2020-01-22 09:23:37

Let's do some further cleaning in the data to push the boundaries for PRS
prediction.

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
prs = read.csv('/Volumes/NCR/reference/merged_NCR_1KG_PRS_12192019.csv')
data_prs = merge(df, prs, by='MRN', all.x=F, all.y=F)
```

Let's merge in self-identified race and ethnicity, as well as family IDs, so
that we can define WNH:

```r
demo = read.csv('prs_demo.csv')
data_prs = merge(data_prs, demo, by='MRN', all.x=F, all.y=F)
library(ggplot2)
sp = ggplot(data_prs, aes(x=PC01, y=PC02, col=population_self)) + geom_point()
sp
```

![](images/2020-01-22-10-32-32.png)

```r
data_prs$isWNH = F
idx = data_prs$PC01 < -.033 & data_prs$PC02 < -.018
data_prs[idx,]$isWNH = T
sp + geom_hline(yintercept=-.018) + geom_vline(xintercept=-.033)
```

![](images/2020-01-22-10-41-45.png)

Let's now select a single kid per family. I was going to do it based on
neuroimage, but I think the criteria should be:

1) Biggest total slope (to increase variance)
2) Most clinical data points (better characterization of the slope)
3) Most scans (higher chances of a good scan)

```r
data_prs$bestInFamily = F
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
```

We're down to 252 subjects if we select only the best in family, 154 if keeping
WNH only (we have 247 out of the 393 as WNH). Now, let's see what PRS can help
predict. I'll start with continuous metrics:

```r
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_') |
                                   grepl(colnames(data_prs), pattern='PC')]
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        pdf(sprintf('~/tmp/%s.pdf', phen))
        data_prs$y = data_prs[, phen]
        use_me = !is.na(data_prs[,]$y) & data_prs$isWNH & data_prs$bestInFamily
        X1 = prcomp(data_prs[use_me, prs_var_names[1:12]])$x
        colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
        X = cbind(X1, as.matrix(data_prs[use_me, prs_var_names[13:22]]))
        X = scale(X)
        Y = data_prs[use_me, ]$y

        # trying to create future classes balanced
        library(caret)
        set.seed(3456)
        nFolds = 10
        trainIndex <- createFolds(Y, k = nFolds)
        # format understood by glmnet
        foldid = vector(length=length(Y))
        for (k in 1:nFolds) {
            foldid[trainIndex[sprintf('Fold%02d', k)][[1]]] = k
        }

        cv1 = cv.glmnet(X, Y, foldid=foldid, alpha=1)
        cv.5 = cv.glmnet(X, Y, foldid=foldid, alpha=.5)
        cv0 = cv.glmnet(X, Y, foldid=foldid,alpha=0)

        mse = mean((Y-mean(Y))**2)
        par(mfrow=c(2,2))
        plot(cv1, main='Lasso'); abline(h=mse, col='green')
        plot(cv.5, main='ENet .5'); abline(h=mse, col='green')
        plot(cv0, main='Ridge'); abline(h=mse, col='green')
        plot(log(cv1$lambda), cv1$cvm, pch=19, col="red", xlab="log(Lambda)",
                 ylab=cv1$name, main=phen)
        points(log(cv.5$lambda), cv.5$cvm, pch=19, col="grey")
        points(log(cv0$lambda), cv0$cvm,pch=19, col="blue")
        legend("topleft", legend=c('Lasso', "Enet .5", "Ridge"), pch=19,
               col=c("red", "grey", "blue"))
        abline(h=mse, col='green')
        dev.off()
    }
}
```

The only situation when we approached the mean error was when predicting HI rate
for GE3 (but GE4 is the same):

![](images/2020-01-22-12-56-42.png)

Let's see what variables were useful there (LASSO):

```r
sx = 'hi'
min_sx = 3
phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
data_prs$y = data_prs[, phen]
use_me = !is.na(data_prs[,]$y) & data_prs$isWNH & data_prs$bestInFamily
X1 = prcomp(data_prs[use_me, prs_var_names[1:12]])$x
colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
X = cbind(X1, as.matrix(data_prs[use_me, prs_var_names[13:22]]))
X = scale(X)
Y = data_prs[use_me, ]$y

# trying to create future classes balanced
library(caret)
set.seed(3456)
nFolds = 10
trainIndex <- createFolds(Y, k = nFolds)
# format understood by glmnet
foldid = vector(length=length(Y))
for (k in 1:nFolds) {
   foldid[trainIndex[sprintf('Fold%02d', k)][[1]]] = k
}

cv1 = cv.glmnet(X, Y, foldid=foldid, alpha=1)
coef(cv1, s = "lambda.min")
```

```
23 x 1 sparse Matrix of class "dgCMatrix"
                       1
(Intercept) -0.325089964
PRSPC01      .          
PRSPC02     -0.045668149
PRSPC03      .          
PRSPC04     -0.010149552
PRSPC05      .          
PRSPC06      0.019601452
PRSPC07      .          
PRSPC08      .          
PRSPC09      .          
PRSPC10      .          
PRSPC11      0.061886196
PRSPC12      .          
PC01         .          
PC02        -0.041012941
PC03         .          
PC04         .          
PC05        -0.108348517
PC06         .          
PC07         .          
PC08         0.008327589
PC09         .          
PC10         .          
```

Nothing earth-shattering. But let me also check that the things I assumed in
previous models are still correct here. First, do the PC of the PRS perform
better? Than, are the population PCs doing anything, since we're focusing on
WNH?

```r
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_') |
                                   grepl(colnames(data_prs), pattern='PC')]
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        pdf(sprintf('~/tmp/%s.pdf', phen))
        data_prs$y = data_prs[, phen]
        use_me = !is.na(data_prs[,]$y) & data_prs$isWNH & data_prs$bestInFamily
        X = scale(data_prs[use_me, prs_var_names])
        Y = data_prs[use_me, ]$y

        # trying to create future classes balanced
        library(caret)
        set.seed(3456)
        nFolds = 10
        trainIndex <- createFolds(Y, k = nFolds)
        # format understood by glmnet
        foldid = vector(length=length(Y))
        for (k in 1:nFolds) {
            foldid[trainIndex[sprintf('Fold%02d', k)][[1]]] = k
        }

        cv1 = cv.glmnet(X, Y, foldid=foldid, alpha=1)
        cv.5 = cv.glmnet(X, Y, foldid=foldid, alpha=.5)
        cv0 = cv.glmnet(X, Y, foldid=foldid,alpha=0)

        mse = mean((Y-mean(Y))**2)
        par(mfrow=c(2,2))
        plot(cv1, main='Lasso'); abline(h=mse, col='green')
        plot(cv.5, main='ENet .5'); abline(h=mse, col='green')
        plot(cv0, main='Ridge'); abline(h=mse, col='green')
        plot(log(cv1$lambda), cv1$cvm, pch=19, col="red", xlab="log(Lambda)",
                 ylab=cv1$name, main=phen)
        points(log(cv.5$lambda), cv.5$cvm, pch=19, col="grey")
        points(log(cv0$lambda), cv0$cvm,pch=19, col="blue")
        legend("topleft", legend=c('Lasso', "Enet .5", "Ridge"), pch=19,
               col=c("red", "grey", "blue"))
        abline(h=mse, col='green')
        dev.off()
    }
}
```

![](images/2020-01-22-13-05-23.png)
![](images/2020-01-22-13-05-49.png)

It doesn't look as good anymore, and there's 2, maybe even just one variable
significant...

```r
sx = 'hi'
min_sx = 3
phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
data_prs$y = data_prs[, phen]
use_me = !is.na(data_prs[,]$y) & data_prs$isWNH & data_prs$bestInFamily
X = scale(data_prs[use_me, prs_var_names])
Y = data_prs[use_me, ]$y
cv1 = cv.glmnet(X, Y, foldid=foldid, alpha=1)
coef(cv1, s = "lambda.min")
```

In fact, for HI it was only the intercept. For inatt we got intercept and
ADHDeur_PRS0.000500, so nothing fancy. Since the population PCs are not helping,
let's see if they're are detriment:

```r
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        pdf(sprintf('~/tmp/%s.pdf', phen))
        data_prs$y = data_prs[, phen]
        use_me = !is.na(data_prs[,]$y) & data_prs$isWNH & data_prs$bestInFamily
        X1 = prcomp(data_prs[use_me, prs_var_names[1:12]])$x
        colnames(X1) = sapply(1:ncol(X1), function(x) sprintf('PRSPC%02d', x))
        X = scale(X1)
        Y = data_prs[use_me, ]$y

        # trying to create future classes balanced
        library(caret)
        set.seed(3456)
        nFolds = 10
        trainIndex <- createFolds(Y, k = nFolds)
        # format understood by glmnet
        foldid = vector(length=length(Y))
        for (k in 1:nFolds) {
            foldid[trainIndex[sprintf('Fold%02d', k)][[1]]] = k
        }

        cv1 = cv.glmnet(X, Y, foldid=foldid, alpha=1)
        cv.5 = cv.glmnet(X, Y, foldid=foldid, alpha=.5)
        cv0 = cv.glmnet(X, Y, foldid=foldid,alpha=0)

        mse = mean((Y-mean(Y))**2)
        par(mfrow=c(2,2))
        plot(cv1, main='Lasso'); abline(h=mse, col='green')
        plot(cv.5, main='ENet .5'); abline(h=mse, col='green')
        plot(cv0, main='Ridge'); abline(h=mse, col='green')
        plot(log(cv1$lambda), cv1$cvm, pch=19, col="red", xlab="log(Lambda)",
                 ylab=cv1$name, main=phen)
        points(log(cv.5$lambda), cv.5$cvm, pch=19, col="grey")
        points(log(cv0$lambda), cv0$cvm,pch=19, col="blue")
        legend("topleft", legend=c('Lasso', "Enet .5", "Ridge"), pch=19,
               col=c("red", "grey", "blue"))
        abline(h=mse, col='green')
        dev.off()
    }
}
```

Nope, nothing there either... let's just make some plots to be sure out slope
variables look fine:

```r
par(mfrow=c(3,2))
hist(df$slope_inatt, breaks=20, main='inatt')
hist(df$slope_hi, breaks=20, main='hi')
idx = df$base_inatt>=3 | df$base_hi>=3
hist(df[idx,]$slope_inatt, breaks=20,
     main=sprintf('inatt (either >=3, n=%d)', sum(idx)))
hist(df[idx,]$slope_hi, breaks=20,
     main=sprintf('hi (either >=3, n=%d)', sum(idx)))
idx = df$base_inatt>=4 | df$base_hi>=4
hist(df[idx,]$slope_inatt, breaks=20,
     main=sprintf('inatt (either >=4, n=%d)', sum(idx)))
hist(df[idx,]$slope_hi, breaks=20,
     main=sprintf('hi (either >=4, n=%d)', sum(idx)))
```

![](images/2020-01-22-13-14-29.png)

They make sense... and here's how they look after winsorizing:

```r
par(mfrow=c(3, 4))
for (sx in c('inatt', 'hi', 'total')) {
    phen = sprintf('slope_%s', sx)
    plot(sort(data_prs[, phen]), main=phen, ylim=c(-4, 2.5))
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        plot(sort(data_prs[, phen]), main=phen, ylim=c(-4, 2.5))
    }
}
```

![](images/2020-01-22-13-19-44.png)

So, winsorizing seems to be working fine too.

## Univariate

Philip suggested I should try a simple model y ~ PRS + PCs, and then maybe add
age and sex. But do each PRS individually, and see which ones are significant
related to the different outcomes. Not use all PRS together. This will make it
easier to use everyone later, if not restricting everyone to one in the family.

```r
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
covars = '+ PC01 + PC02 + PC03 + PC04 + PC05 + PC06 + PC07 + PC08 + PC09 + PC10 + base_age + sex.x'
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_WNH_PCsAgeSex_lm.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen_res = c()
        for (prs in prs_var_names) {
            use_me = !is.na(data_prs[, phen]) & data_prs$bestInFamily & data_prs$isWNH
            fm_str = paste(phen, "~", prs, covars, sep="")
            fit = lm(as.formula(fm_str), data=data_prs[use_me, ])
            # assuming interesting variable is always first one
            temp = c(summary(fit)$coefficients[2, ], summary(fit)$r.squared,
                     summary(fit)$adj.r.squared)
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
colnames(hold)[5:6] = c('R2', 'adjR2')
write.csv(hold, file=out_fname, row.names=F)
```

I tried a few variations of the code above, outputing the different CSV files,
and then copied everything into an Excel file called
prs_univariate_results.xslx.

For lme it changed a bit:

```r
library(nlme)
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
covars = ' + PC01 + PC02 + PC03 + PC04 + PC05 + PC06 + PC07 + PC08 + PC09 + PC10 + base_age + sex.x'
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_WNH_PCsAgeSex_lme.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen_res = c()
        for (prs in prs_var_names) {
            use_me = !is.na(data_prs[, phen]) & data_prs$isWNH
            fm_str = paste(phen, "~", prs, covars, sep="")
            fit = lme(as.formula(fm_str), ~1|FAMID, data=data_prs[use_me, ],
                      control=lmeControl(tolerance=1e-100, returnObject=TRUE))
            # assuming interesting variable is always first one
            temp = c(summary(fit)$tTable[2, ], summary(fit)$BIC,
                     summary(fit)$AIC)
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
colnames(hold)[6:7] = c('BIC', 'AIC')
write.csv(hold, file=out_fname, row.names=F)
```

## Binary groups

Let's see if we do any better if we binarize the data. Let's start with the
median value as the cutoff. But we could even optimize that, if the goal here is
really stacking up the deck...

We start with one per family as before:

<!-- ```r
hold = c()
prs_var_names = colnames(data_prs)[grepl(colnames(data_prs), pattern='ADHDeur_')]
covars = '+ PC01 + PC02 + PC03 + PC04 + PC05 + PC06 + PC07 + PC08 + PC09 + PC10 + base_age + sex.x'
out_fname = '~/data/baseline_prediction/prs_start/univar_prs_WNH_PCsAgeSex_lm.csv'
for (sx in c('inatt', 'hi', 'total')) {
    for (min_sx in c(0, 3, 4)) {
        phen = sprintf('slope_%s_GE%d_wp05', sx, min_sx)
        phen_res = c()
        for (prs in prs_var_names) {
            use_me = !is.na(data_prs[, phen]) & data_prs$bestInFamily & data_prs$isWNH
            fm_str = paste(phen, "~", prs, covars, sep="")
            fit = lm(as.formula(fm_str), data=data_prs[use_me, ])
            # assuming interesting variable is always first one
            temp = c(summary(fit)$coefficients[2, ], summary(fit)$r.squared,
                     summary(fit)$adj.r.squared)
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
colnames(hold)[5:6] = c('R2', 'adjR2')
write.csv(hold, file=out_fname, row.names=F)
```


fm = as.formula(sprintf("%s ~ OLS_%s_slope + (1|nuclearFamID)", colnames(df)[i], sx))
            print(fm)
            fit = try(glmer(fm, data=df[idx2, ], na.action=na.omit, family = binomial(link = "logit")))
            print(summary(fit)$coefficients[2,]) -->



# TODO
* binary groups
* try persistence PRS