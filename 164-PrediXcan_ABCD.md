# 2020-12-21 14:15:46

Let's play a bit with the ABCD transcript imputation. First, let's load the data
and save it as RDS, because the files are quite big.

Now that I have the files downloaded, let's see how many subjects have
Freesurfer data, and then we can worry about the Qc parameters.

I'll also need to remove anyone not WNH because of the imputation model.

```r
# start with everyone who was imputed
imp = readRDS('~/data/expression_impute/results/ABCD_v201_ACC_predict_1KG_mashr.rds')
# keep only WNH
mds.cluster = read.table('~/data/expression_impute/HM3_b37mds.mds', head=1)
quartz()
colors=rep("red",length(mds.cluster$C1));
colors[which(mds.cluster$FID == "CEU")] <- "lightblue";
colors[which(mds.cluster$FID == "CHB")] <- "brown";
colors[which(mds.cluster$FID == "YRI")] <- "yellow";
colors[which(mds.cluster$FID == "TSI")] <- "green";
colors[which(mds.cluster$FID == "JPT")] <- "purple";
colors[which(mds.cluster$FID == "CHD")] <- "orange";
colors[which(mds.cluster$FID == "MEX")] <- "grey50";
colors[which(mds.cluster$FID == "GIH")] <- "black";
colors[which(mds.cluster$FID == "ASW")] <- "darkolivegreen";
colors[which(mds.cluster$FID == "LWK")] <- "magenta";
colors[which(mds.cluster$FID == "MKK")] <- "darkblue";
plot(mds.cluster$C2, mds.cluster$C1, col=colors,
         ylab="Dimension 1", xlab="Dimension 2",pch=20)
legend("topleft", c("CEU", "CHB", "YRI", "TSI", "JPT", "CHD",
                     "MEX", "GIH", "ASW","LWK", "MKK", 'ABCD'),
       fill=c("lightblue", "brown", "yellow", "green", "purple",
              "orange", "grey50", "black", "darkolivegreen", "magenta",
              "darkblue", 'red'))
```

![](images/2020-12-21-19-46-12.png)

I don't want to lose that many people, so let's cap it below and above zero for
now to keep it simple.

```r
imwnh = mds.cluster[mds.cluster$C1 < 0 & mds.cluster$C2 > 0, 'IID']
# 6233 out of 11301
imp$IID2 = sapply(imp$IID, function(x) paste(strsplit(x=x,
                                                      split='_')[[1]][2:3],
                                            collapse='_'))
imp.wnh = imp[imp$IID2 %in% imwnh, ]
# 6031 left, let's see who has freesurfer
junk = read.delim('~/data/expression_impute/abcd_smrip101.txt')
fs = read.delim('~/data/expression_impute/abcd_smrip101.txt', skip=2, head=0)
colnames(fs) = colnames(junk)
fs.imp = fs[fs$src_subject_id %in% imp.wnh$IID2, ]
# 5897 imputed WNH left with Freesurfer
junk = read.delim('~/data/expression_impute/freesqc01.txt')
# quick hack to avoid description line
fsqc = read.delim('~/data/expression_impute/freesqc01.txt', skip=2, head=0)
colnames(fsqc) = colnames(junk)
good_subj = fsqc[which(fsqc$fsqc_qc == 1), 'src_subject_id']
# 11076 out of 11556 are marked as good
fs.use = fs.imp[fs.imp$src_subject_id %in% good_subj, ]
# 5714 usable subjects
brain = merge(fs.use, fsqc, by='src_subject_id')

data = merge(brain, imp, by.x='src_subject_id', by.y='IID2', all.x=F, all.y=F)
```

Now it's just a matter of running the regressions:

```r
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(data[, grex_vars], method = pp_order)
X = predict(pp, data[, grex_vars])

library(bestNormalize)
Xnorm = X
for (v in 1:ncol(Xnorm)) {
    if ((v %% 100)==0) {
        print(sprintf('%d / %d', v, ncol(Xnorm)))
    }
    bn = orderNorm(Xnorm[, v])
    Xnorm[, v] = bn$x.t
}

# clean up some of the brain as well
# whole rhACC
data$rh_ACC_vol = data$smri_vol_cdk_cdacaterh + data$smri_vol_cdk_rracaterh
data$rh_ACC_thick = data$smri_thick_cdk_cdacaterh + data$smri_thick_cdk_rracaterh
data$rh_ACC_area = data$smri_area_cdk_cdacaterh + data$smri_area_cdk_rracaterh
# whole caACC
data$caACC_vol = data$smri_vol_cdk_cdacatelh + data$smri_vol_cdk_cdacaterh
data$caACC_area = data$smri_area_cdk_cdacaterh + data$smri_area_cdk_cdacatelh
data$caACC_thick = data$smri_thick_cdk_cdacaterh + data$smri_thick_cdk_cdacatelh
# whole roACC
data$roACC_vol = data$smri_vol_cdk_rracaterh + data$smri_vol_cdk_rracatelh
data$roACC_area = data$smri_area_cdk_rracaterh + data$smri_area_cdk_rracatelh
data$roACC_thick = data$smri_thick_cdk_rracaterh + data$smri_thick_cdk_rracatelh
# whole ACC
data$ACC_vol = data$roACC_vol + data$caACC_vol
data$ACC_area = data$roACC_area + data$caACC_area
data$ACC_thick = data$roACC_thick + data$caACC_thick

brain_vars = c()
for (p in c('area', 'thick', 'vol')) {
    for (bv in c('rh_ACC', 'caACC', 'roACC', 'ACC')) {
        brain_vars = c(brain_vars, sprintf('%s_%s', bv, p))
    }
    for (bv in c('smri_%s_cdk_cdacaterh', 'smri_%s_cdk_rracaterh')) {
        brain_vars = c(brain_vars, sprintf(bv, p))
    }
}

for (v in brain_vars) {
    m = mean(data[, v], na.rm=T)
    s = sd(data[, v], na.rm=T)
    data[which(data[, v] > m + 3*s), v] = NA
    data[which(data[, v] < m - 3*s), v] = NA
}
dataRaw = data
dataRaw[, colnames(X)] = X
dataNorm = data
dataNorm[, colnames(X)] = Xnorm

library(vows)
form = (~ meta$Diagnosis + mydata$PC1 + mydata$PC2 + mydata$PC3 + mydata$PC4
        + mydata$PC5 + mydata$PC6 + mydata$PC7 + mydata$PC9 + mydata$PC11
        + mydata$PC13 + mydata$PC14 + mydata$PC17)
fm_str = ("~ data$rh_ACC_area + data$fsqc_qu_motion + data$fsqc_qu_pialover + data$fsqc_qu_wmunder + data$fsqc_qu_inhomogeneity + data$fsqc_qu_artifact")
res = summary(lm.mp(X, as.formula(fm_str)))
junk = data.frame(P.Value=res$pvalue['meta$DiagnosisCase',],
                  t=res$tstat['meta$DiagnosisCase',],
                  adj.P.Val=p.adjust(res$pvalue['meta$DiagnosisCase',],
                                     method='fdr'))
res_dx = merge(junk, meta_iso, by.x=0, by.y='id1')
res_dx = res_dx[order(res_dx$P.Value), ]
``` 
