# 2019-05-02 11:38:00

Philip suggested to collapse the DTI tract results across hemispheres, to reduce
the number of comparisons, and hopefully keep the results. So, what do we get?

```r
a = read.csv('~/data/heritability_change/dti_JHUtracts_residNoSex_OLS_naSlopes133.csv')
b = cbind(b, rowMeans(a[,c('ad_1', 'ad_2')]))
b = cbind(b, rowMeans(a[,c('ad_3', 'ad_4')]))
b = cbind(b, rowMeans(a[,c('ad_5', 'ad_6')]))
b = cbind(b, rowMeans(a[,c('ad_7', 'ad_8')]))
b = cbind(b, rowMeans(a[,c('ad_11', 'ad_12')]))
b = cbind(b, rowMeans(a[,c('ad_13', 'ad_14')]))
b = cbind(b, rowMeans(a[,c('ad_15', 'ad_16')]))
b = cbind(b, rowMeans(a[,c('ad_17', 'ad_18')]))
b = cbind(b, rowMeans(a[,c('ad_19', 'ad_20')]))
b = cbind(b, rowSums(a[,c('ad_9', 'ad_10')]))
b = cbind(b, rowMeans(a[,c('rd_1', 'rd_2')]))
b = cbind(b, rowMeans(a[,c('rd_3', 'rd_4')]))
b = cbind(b, rowMeans(a[,c('rd_5', 'rd_6')]))
b = cbind(b, rowMeans(a[,c('rd_7', 'rd_8')]))
b = cbind(b, rowMeans(a[,c('rd_11', 'rd_12')]))
b = cbind(b, rowMeans(a[,c('rd_13', 'rd_14')]))
b = cbind(b, rowMeans(a[,c('rd_15', 'rd_16')]))
b = cbind(b, rowMeans(a[,c('rd_17', 'rd_18')]))
b = cbind(b, rowMeans(a[,c('rd_19', 'rd_20')]))
b = cbind(b, rowSums(a[,c('rd_9', 'rd_10')]))
colnames(b)[9:28] = c('ad_atr', 'ad_cst', 'ad_cin', 'ad_hip', 'ad_ifo',
                      'ad_ilf', 'ad_slf', 'ad_unc', 'ad_slf2', 'ad_for',
                      'rd_atr', 'rd_cst', 'rd_cin', 'rd_hip', 'rd_ifo',
                      'rd_ilf', 'rd_slf', 'rd_unc', 'rd_slf2', 'rd_for')
write.csv(b, file='~/data/heritability_change/dti_JHUtractsComb_residNoSex_OLS_naSlopes133.csv',
          row.names=F, na='', quote=F)
```

Let's try combining even further:

```r
a = read.csv('~/data/heritability_change/dti_JHUtracts_residNoSex_OLS_naSlopes133.csv')
b = a[,c(1, 2, 63:68)]
for (m in c('ad', 'rd', 'fa')) {
    b = cbind(b, rowMeans(a[,c(sprintf('%s_1', m), sprintf('%s_2', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_atf', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_3', m), sprintf('%s_4', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_cst', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_5', m), sprintf('%s_6', m),
                            sprintf('%s_7', m), sprintf('%s_8', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_cin', m)
    b = cbind(b, rowSums(a[,c(sprintf('%s_9', m), sprintf('%s_10', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_for', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_11', m), sprintf('%s_12', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_ifo', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_13', m), sprintf('%s_14', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_ilf', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_15', m), sprintf('%s_16', m),
                            sprintf('%s_19', m), sprintf('%s_20', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_slf', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_17', m), sprintf('%s_18', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_unc', m)
}
write.csv(b, file='~/data/heritability_change/dti_JHUtractsComb2_residNoSex_OLS_naSlopes133.csv',
          row.names=F, na='', quote=F)
```

OK, so that gives us results for uncinate AD and RD, and RD even survives
Bonferroni. AD might survive Meff, but let's check association first.

```r
library(nlme)
data = read.csv('~/data/heritability_change/dti_JHUtracts_residNoSex_OLS_naSlopesAndBaseline283.csv')
data$sex = as.factor(data$sex)
b = data[, c(1, 2, 63:68, 129)]
a = data
for (m in c('ad', 'rd')) {
    b = cbind(b, rowMeans(a[,c(sprintf('%s_1', m), sprintf('%s_2', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_atf', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_3', m), sprintf('%s_4', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_cst', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_5', m), sprintf('%s_6', m),
                            sprintf('%s_7', m), sprintf('%s_8', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_cin', m)
    b = cbind(b, rowSums(a[,c(sprintf('%s_9', m), sprintf('%s_10', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_for', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_11', m), sprintf('%s_12', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_ifo', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_13', m), sprintf('%s_14', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_ilf', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_15', m), sprintf('%s_16', m),
                            sprintf('%s_19', m), sprintf('%s_20', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_slf', m)
    b = cbind(b, rowMeans(a[,c(sprintf('%s_17', m), sprintf('%s_18', m))]))
    colnames(b)[ncol(b)] = sprintf('%s_unc', m)
}
tract_names = colnames(b)[10:ncol(b)]
data = b
out_fname = '~/data/heritability_change/assoc_LME_JHUtractsComb2_naSlopes283.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline', 'DX', 'DX2')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-try(lme(as.formula(fm_str), data, ~1|famID, na.action=na.omit))
        if (length(model1) > 1) {
            temp<-summary(model1)$tTable
            a<-as.data.frame(temp)
            a$formula<-fm_str
            a$target = i
            a$predictor = j
            a$term = rownames(temp)
            hold=rbind(hold,a)
        } else {
            hold=rbind(hold, NA)
        }
    }
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX=='ADHD', ]
out_fname = '~/data/heritability_change/assoc_LME_JHUtractsComb2_naSlopes283_dx1.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|famID, na.action=na.omit))
        if (length(model1) > 1) {
            temp<-summary(model1)$tTable
            a<-as.data.frame(temp)
            a$formula<-fm_str
            a$target = i
            a$predictor = j
            a$term = rownames(temp)
            hold=rbind(hold,a)
        } else {
            hold=rbind(hold, NA)
        }
    }
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX2=='ADHD', ]
out_fname = '~/data/heritability_change/assoc_LME_JHUtractsComb2_naSlopes283_dx2.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|famID, na.action=na.omit))
        if (length(model1) > 1) {
            temp<-summary(model1)$tTable
            a<-as.data.frame(temp)
            a$formula<-fm_str
            a$target = i
            a$predictor = j
            a$term = rownames(temp)
            hold=rbind(hold,a)
        } else {
            hold=rbind(hold, NA)
        }
    }
}
write.csv(hold, out_fname, row.names=F)
```

No, unfortunately putting the hemispheres together wipes out the association, at
least when using LME. What if I use straight up regression?

```r
out_fname = '~/data/heritability_change/assoc_JHUlabelsComb2_naSlopes283.csv'
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-lm(as.formula(fm_str), data, na.action=na.omit)
        temp<-summary(model1)$coefficients
        a<-as.data.frame(temp)
        a$formula<-fm_str
        a$target = i
        a$predictor = j
        a$term = rownames(temp)
        hold=rbind(hold,a)
    }
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX=='ADHD', ]
out_fname = '~/data/heritability_change/assoc_JHUlabelsComb2_naSlopes283dx1.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-lm(as.formula(fm_str), data2, na.action=na.omit)
        temp<-summary(model1)$coefficients
        a<-as.data.frame(temp)
        a$formula<-fm_str
        a$target = i
        a$predictor = j
        a$term = rownames(temp)
        hold=rbind(hold,a)
    }
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX2=='ADHD', ]
out_fname = '~/data/heritability_change/assoc_JHUlabelsComb2_naSlopes283_dx2.csv'
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline')
targets = tract_names
hold=NULL
for (i in targets) {
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex', i, j)
        model1<-lm(as.formula(fm_str), data2, na.action=na.omit)
        temp<-summary(model1)$coefficients
        a<-as.data.frame(temp)
        a$formula<-fm_str
        a$target = i
        a$predictor = j
        a$term = rownames(temp)
        hold=rbind(hold,a)
    }
}
write.csv(hold, out_fname, row.names=F)
```

That didn't help either... let's throw a hail marry and try to go for voxelwise
SOLAR.
