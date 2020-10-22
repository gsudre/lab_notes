# 2020-10-21 16:12:52

Let's go back to the imputation results now that we have some solid results for
gene set analysis (138). We can use everyone, and let's only use the people with
good scans. The postmortem data has an age range from 6.69 to 38.83, so let's
make sure our samples and scans are within the same age range. We start from
expression_impute/gf_1119_09092020.csv and add all their anatomical QC to the
file.

```r
gf = read.csv('~/data/expression_impute/gf_1119_09092020.csv')
scans = read.csv('~/data/expression_impute/all_mprage_qc.csv')
m = merge(gf, scans, all.x=T, all.y=F, by=1)
write.csv(m, file='~/data/expression_impute/gfWithMPRAGE_2456_10212020.csv')
```

I then manually remove everyone without a scan, and then selected the best scan
inside that age range above, based on MPRAGE QC. In the case of a tie, I chose
the eldest scan. That created gfWithMPRAGE_773_10212020.csv. To that, we need to
attach the freesurfer QC and the regional values. For caudate we only have
volume, but for ACC we can do all 3. Ideally volume will work best so we can
keep just one metric, but we'll see. I added the QCs to the same file, and the
Freesurfer metrics. I also created bilateral versions just in case.


# 2020-10-22 06:21:22

Let's go ahead and run a very similar analysis as in limma. I won't use the same
package tools, but at least the same statistical model. For that, I'll need to
check for correlated PCs, so we use the same approach of denoising here. The
tricky thing is to make sure the PCs are not related to any of the phenotypes
I'll try.

The question here is whether I should include scan-specific variables ni the
model. For the RNAseq proper, we did a PCA in the counts and checked the PC
correlation with rna-metrics but also biological. Here, our phenotype is the
brain (not DX), and a few variables (such as QC) wouldn't be related to the
dependent variable (imputed scores). Do we add them to the model regardless? Or
check their correlation first? Need to check both ways!

Philip suggested to just use good scans and potentially regress out the QC
variables if needed. Let's go with that.

I also just noticed I didn't remove based on the age range, so I'm recreating
the gf (gfWithMPRAGE_632_10222020.csv)

```r
gf = read.csv('~/data/expression_impute/gfWithMPRAGE_773_10212020.csv')
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
imp_data = merge(imp_data, gf, by.x='IID', by.y='Subject.Code...Subjects',
                 all.x=F, all.y=F)


# TODO
 * run models with and without neuroimaging-related covariates


