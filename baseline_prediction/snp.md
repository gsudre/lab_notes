# 2018-09-13 13:49:14

The idea is to use subsets of SNPs from the PGC 2017 GWAS, which will be similar
to doing PRS, except that we'll have information on individual SNPs, and won't
be using the beta weights.

let's use the same dataset we used for PRS, after cleaning:

```bash
[sudregp@cn3443 geno3]$ pwd
/home/sudregp/data/prs/geno3
plink --bfile merged_noDups_clean_flipped --out merged_noDups_clean_flipped.tab --recodeA

awk '$9 < 1e-04' ~/pgc2017/adhd_jun2017 | cut -f 2 - > ~/tmp/snps_pgcadhd2017_1e04.txt
```

That generates several subsets of SNPs. Then, in R:

```r
library(h2o)
h2o.init()
df = h2o.importFile('/Users/sudregp/data/baseline_prediction/merged_noDups_clean_flipped.tab.raw')
snps = read.table('~/data/baseline_prediction/snps_pgcadhd2017_1e04.txt', colClasses='character')[,1]
cnames = sapply(colnames(df), function(x) { gsub('_[GACT]$','',x) })
colnames(df) = cnames
keep_me = c('IID'); for (i in snps) { if (i %in% cnames) { keep_me = c(keep_me, i) } }
df2 = df[, keep_me]
write.csv(as.data.frame(df2),file='~/data/baseline_prediction/geno3_snps1e04_09132018.csv', row.names=F)
```

And then repeat that for the other thresholds.

```r
for (p in 5:9) {
    print(p)
    snp_fname = sprintf('~/data/baseline_prediction/snps_pgcadhd2017_1e0%d.txt', p)
    snps = read.table(snp_fname, colClasses='character')[,1]
    keep_me = c('IID'); for (i in snps) { if (i %in% cnames) { keep_me = c(keep_me, i) } }
    df2 = df[, keep_me]
    out_fname = sprintf('~/data/baseline_prediction/geno3_snps1e0%d_09132018.csv', p)
    write.csv(as.data.frame(df2),file=out_fname, row.names=F)
}
```

I saved everything as CSV to make it easier to access in the future. Now, the
only processing that makes sense ot use here is raw and pca, as the univariate
was replaced by the GWAS results. 

Then, we just need to create the gf file:

```r
prs = read.csv('/Volumes/Shaw/prs/FINAL_FILES_08242018/REGRESSION/PRS2017_geno3_1KG_noFlip_genop05MAFbtp01rsbtp9.csv')
nsb_long = as.vector(df2$IID)
need_reformat = grepl('-', nsb_long)
clean_nsb = sapply(nsb_long[need_reformat], function(x) strsplit(strsplit(x, '-')[[1]][3], '@')[[1]][1])
nsb = nsb_long
nsb[need_reformat] = clean_nsb
df2 = h2o.cbind(df2, as.h2o(as.numeric(nsb)))
```



for (i in colnames(snps)) { df[, i] = as.factor(df[,i])}
cnames = sapply(colnames(df), function(x) { gsub('_[GACT]$','',x) })

keep_me = c('IID'); for (i in snps) { if (i %in% colnames(df)) { keep_me = c(keep_me, i) } }
df2 = df[, keep_me]
> dim(df2)
[1] 997   4
```