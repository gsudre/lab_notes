# 2021-02-09 13:54:33

First, let's reproduce the previous plotsw ehad using the meta-analysis for the
cortex, from note 184:

```r
do_boot_corrs = function(both_res, log2FC_col, method) {
    corrs = c()
    nperms = 10000
    set.seed(42)
    options(warn=-1)  # remove annoying spearman warnings
    for (p in 1:nperms) {
        idx = sample(nrow(both_res), replace = T)
        corrs = c(corrs, cor.test(both_res[idx, 'log2FoldChange'],
                                  both_res[idx, log2FC_col],
                                  method=method)$estimate)
    }
    return(corrs)
}

meta = readRDS('~/data/post_mortem/aad6469_Gandal_SM_Data-Table-S1_micro.rds')

st = 'protein_coding'
met = 'spearman'
load('~/data/post_mortem/DGE_02082021.RData')
dge = as.data.frame(dge_acc[[st]][['res']])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)

corrs = list()
disorders = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD', 'IBD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.beta_log2FC', d), met)
}
mylim = max(abs(unlist(corrs)))
quartz()
boxplot(corrs, ylim=c(-mylim, mylim), ylab=sprintf('%s correlation', met),
        main=sprintf('ACC %s (n=%d)', st, nrow(both_res)))
abline(h=0, col='red')
# calculating p-value
for (d in disorders) {
    if (median(corrs[[d]]) > 0) {
        pval = sum(corrs[[d]] <= 0) / 10000
    } else {
        pval = sum(corrs[[d]] >= 0) / 10000
    }
    cat(d, 'pval = ', pval, '\n')
}
```

![](images/2021-02-09-14-06-14.png)

Now, let's do something similar, but now using the data from the ACC paper
Kwangmi sent: https://www.nature.com/articles/s41386-020-00949-5

```r
# read.xls was taking forever. Had to crop them into individual CSVs, from the 
# original ~/data/post_mortem/41386_2020_949_MOESM3_ESM.xlsx
library(gdata)
meta1 = read.csv('~/tmp/BD.csv')
meta2 = read.csv('~/tmp/SCZ.csv')
meta3 = read.csv('~/tmp/MDD.csv')
m = merge(meta1[, c('Ensemble.gene.ID', 'log2FoldChange')],
          meta2[, c('Ensemble.gene.ID', 'log2FoldChange')],
          by='Ensemble.gene.ID', all.x=T, all.y=T, suffix=c('.BD', '.SCZ'))
m = merge(m, meta3[, c('Ensemble.gene.ID', 'log2FoldChange')],
          by='Ensemble.gene.ID', all.x=T, all.y=T)
colnames(m)[ncol(m)] = 'log2FoldChange.MDD'
saveRDS(m, file='~/data/post_mortem/ACC_other_disorders.rds')
```

Now we try this again, now with the new data. Note that the new data only has
genes that are nominally significant at p < .05, so this might be harder:

```r
meta = readRDS('~/data/post_mortem/ACC_other_disorders.rds')

st = 'protein_coding'
met = 'spearman'
load('~/data/post_mortem/DGE_02082021.RData')
dge = as.data.frame(dge_acc[[st]][['res']])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by.x='ensembl_gene_id', by.y='Ensemble.gene.ID',
                 all.x=F, all.y=F)

corrs = list()
disorders = c('BD', 'SCZ', 'MDD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('log2FoldChange.%s', d), met)
}
mylim = max(abs(unlist(corrs)))
quartz()
boxplot(corrs, ylim=c(-mylim, mylim), ylab=sprintf('%s correlation', met),
        main=sprintf('ACC %s (n=%d)', st, nrow(both_res)))
abline(h=0, col='red')
# calculating p-value
for (d in disorders) {
    if (median(corrs[[d]]) > 0) {
        pval = sum(corrs[[d]] <= 0) / 10000
    } else {
        pval = sum(corrs[[d]] >= 0) / 10000
    }
    cat(d, 'pval = ', pval, '\n')
}
```

![](images/2021-02-09-18-05-59.png)

So, this looks good. Now we just need to find data for the Caudate.