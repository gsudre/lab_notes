# 2021-01-29 14:00:01

Philip sent me this interesting paper:

Transcriptomic organization of the human brain in post-traumatic stress disorder
https://doi.org/10.1038/s41593-020-00748-7

and it has two interesting analysis I could try to replicate. One of them is
related to how the DGE results correlate with expression analysis from other
studies. The other is related to a TWAS, which would involve revisiting our
FUSION or METAXcan results. We could also try tackling a similar issue by using
the MESC tool for mediation...

Let's explore the correlation idea first.

```r
library(gdata)
# meta = read.xls('~/data/post_mortem/aad6469_Gandal_SM_Data-Table-S1.xlsx',
                # 'Microarray MetaAnalysis DGE')
# saved as RDS because it took a while to load
meta = readRDS('~/data/post_mortem/aad6469_Gandal_SM_Data-Table-S1_micro.rds')
load('~/data/post_mortem/DGE_01272021.RData')
dge = as.data.frame(dge_acc[['protein_coding']])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)
```

```
r$> dim(both_res)                                                                             
[1] 14852    32

r$> cor.test(both_res$log2FoldChange, both_res$ASD.beta_log2FC, method='spearman')            

        Spearman's rank correlation rho

data:  both_res$log2FoldChange and both_res$ASD.beta_log2FC
S = 2.71e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.2602367 


r$> cor.test(both_res$log2FoldChange, both_res$SCZ.beta_log2FC, method='spearman')            

        Spearman's rank correlation rho

data:  both_res$log2FoldChange and both_res$SCZ.beta_log2FC
S = 1.1379e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.3056628 


r$> cor.test(both_res$log2FoldChange, both_res$BD.beta_log2FC, method='spearman')             

        Spearman's rank correlation rho

data:  both_res$log2FoldChange and both_res$BD.beta_log2FC
S = 1.1827e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
      rho 
0.2776809 


r$> cor.test(both_res$log2FoldChange, both_res$MDD.beta_log2FC, method='spearman')            

        Spearman's rank correlation rho

data:  both_res$log2FoldChange and both_res$MDD.beta_log2FC
S = 4.4619e+11, p-value < 2.2e-16
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
0.07248427 


r$> cor.test(both_res$log2FoldChange, both_res$AAD.beta_log2FC, method='spearman')            

        Spearman's rank correlation rho

data:  both_res$log2FoldChange and both_res$AAD.beta_log2FC
S = 4.8772e+11, p-value = 0.01925
alternative hypothesis: true rho is not equal to 0
sample estimates:
        rho 
-0.01963491 


r$> cor.test(both_res$log2FoldChange, both_res$IBD.beta_log2FC, method='spearman')            

        Spearman's rank correlation rho

data:  both_res$log2FoldChange and both_res$IBD.beta_log2FC
S = 3.8286e+11, p-value = 9.395e-09
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
0.04953456 
```

I get some significant results, but a few things:

 * the paper uses bootstraping to assess significance. Not sure if that's
   built-in to the Spearman correlation function, or if I need to run it myself.
   Maybe using bootstrapping the p-value won't be so inflated?
* I need to see if I get similar high results with Caudate
* Same thing for lncRNA and pseudogenes
* Should also test the RNAseq results in the table

So, I still have 2404 genes in the intersection for lncRNA, and 778 for
pseudogenes. So, still definitely possible to do the analysis. Let's see if the
results change much when using boostrap. It's also not clear whether they used
Spearman as indicated in Figure captions, or Pearson as indicated in the
statistical sections of the method. We can try both.

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

st = 'protein_coding'
met = 'spearman'
dge = as.data.frame(dge_acc[[st]])
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

## Spearman ACC
![](images/2021-02-01-09-47-25.png)
```
ASD pval =  0 
SCZ pval =  0 
BD pval =  0 
MDD pval =  0 
AAD pval =  0.0059 
IBD pval =  0 
```
![](images/2021-01-29-14-57-50.png)
```
ASD pval =  0.3245 
SCZ pval =  0.2743 
BD pval =  0.3955 
MDD pval =  0.0042 
AAD pval =  0.1283 
IBD pval =  0.1122 
```
![](images/2021-01-29-15-02-10.png)
```
ASD pval =  0.3466 
SCZ pval =  0.3359 
BD pval =  0.3401 
MDD pval =  0.3837 
AAD pval =  0.2696 
IBD pval =  0.1304 
```

## Pearson ACC
![](images/2021-01-29-15-19-16.png)
```
ASD pval =  0
SCZ pval =  0
BD pval =  0
MDD pval =  0
AAD pval =  0.1422
IBD pval =  0
```
![](images/2021-01-29-15-18-40.png)
```
ASD pval =  0.2317
SCZ pval =  0.1988
BD pval =  0.4054
MDD pval =  0.0062
AAD pval =  0.013
IBD pval =  0.1475
```
![](images/2021-01-29-15-17-34.png)
```
ASD pval =  0.3365
SCZ pval =  0.4301
BD pval =  0.4927
MDD pval =  0.3279
AAD pval =  0.1995
IBD pval =  0.0418
```

```r
st = 'protein_coding'
met = 'spearman'
dge = as.data.frame(dge_cau[[st]])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)
corrs = list()
disorders = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD', 'IBD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.beta_log2FC', d), met)
}
mylim = max(abs(unlist(corrs)))
boxplot(corrs, ylim=c(-mylim, mylim), ylab=sprintf('%s correlation', met),
        main=sprintf('Caudate %s (n=%d)', st, nrow(both_res)))
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

## Caudate spearman
![](images/2021-02-01-10-07-17.png)
```
ASD pval =  0
SCZ pval =  0.0043
BD pval =  0
MDD pval =  0
AAD pval =  0
IBD pval =  0
```

![](images/2021-01-29-15-09-12.png)
```
ASD pval =  0.1058
SCZ pval =  0.1713
BD pval =  0.3574
MDD pval =  0.2159
AAD pval =  0.0142
IBD pval =  0.4637
```

![](images/2021-01-29-15-08-50.png)
```
ASD pval =  0.1731 
SCZ pval =  0.4971 
BD pval =  0.3343 
MDD pval =  0.335 
AAD pval =  0.2037 
IBD pval =  0.3251 
```

## Caudate pearson
![](images/2021-01-29-15-18-14.png)
```
ASD pval =  0
SCZ pval =  5e-04
BD pval =  7e-04
MDD pval =  0
AAD pval =  0
IBD pval =  0
```
![](images/2021-01-29-15-21-16.png)
```
ASD pval =  0.3414 
SCZ pval =  0.1562 
BD pval =  0.2078 
MDD pval =  0.0505 
AAD pval =  0.0018 
IBD pval =  0.3987 
```
![](images/2021-01-29-15-21-47.png)
```
ASD pval =  0.3764
SCZ pval =  0.268
BD pval =  0.3836
MDD pval =  0.3929
AAD pval =  0.3731
IBD pval =  0.1621
```

What happens if I switch it to a permutation p-value?

```r
do_perm_corrs = function(both_res, log2FC_col, method) {
    corrs = c()
    nperms = 10000
    set.seed(42)
    options(warn=-1)  # remove annoying spearman warnings
    for (p in 1:nperms) {
        idx = sample(nrow(both_res), replace = F)
        corrs = c(corrs, cor.test(both_res[idx, 'log2FoldChange'],
                                  both_res[, log2FC_col],
                                  method=method)$estimate)
    }
    return(corrs)
}

st = 'protein_coding'
met = 'spearman'
dge = as.data.frame(dge_acc[[st]])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)
corrs = list()
disorders = c('ASD', 'SCZ')#, 'BD', 'MDD', 'AAD', 'IBD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_perm_corrs(both_res, sprintf('%s.beta_log2FC', d), met)
}
a = cor.test(both_res[, 'log2FoldChange'], both_res[, log2FC_col],
             method=met)$estimate                                                             
print(sum(corrs[[d]] >= a))
```

Still very significant. OK, what if we remove the AAD and MDD subjects from the
analysis? Does it change the Spearman correlations?

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]

library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021 (1).xlsx',
              'clinical_summary_sheet') 
clean_subjs = df[df$AAD == '', 'original_brain_number']
# clean_subjs = df[df$MDD != 'YES', 'original_brain_number']
data = data[data$original_brain_number %in% clean_subjs, ]

library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]
# data only contains sample metadata, and count_matrix has actual counts

# cleaning up some variables
data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))
data$Diagnosis = factor(data$Diagnosis, levels=c('Control', 'Case'))
data$substance_group = factor(data$substance_group)
data$comorbid_group = factor(data$comorbid_group_update)
data$evidence_level = factor(data$evidence_level)

# removing everything but autosomes
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
# store gene names in geneCounts without version in end of name
tx_meta = data.frame(GENEID = substr(rownames(count_matrix), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
count_matrix = count_matrix[imautosome, ]
tx_meta = tx_meta[imautosome, ]
```

Using the function from 183:

```r
st = 'protein_coding'

res = run_DGE(count_matrix, tx_meta, myregion, st, .05)
met = 'spearman'
dge = as.data.frame(res)
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
meta = readRDS('~/data/post_mortem/aad6469_Gandal_SM_Data-Table-S1_micro.rds')
both_res = merge(dge, meta, by='ensembl_gene_id', all.x=F, all.y=F)
corrs = list()
disorders = c('ASD', 'SCZ', 'BD', 'MDD', 'AAD', 'IBD')
for (d in disorders) {
    cat(d, '\n')
    corrs[[d]] = do_boot_corrs(both_res, sprintf('%s.beta_log2FC', d), met)
}
mylim = max(abs(unlist(corrs)))
boxplot(corrs, ylim=c(-mylim, mylim), ylab=sprintf('%s correlation', met),
        main=sprintf('%s %s (n=%d)', myregion, st, nrow(both_res)))
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

Let's redo Spearman for the 3 datasets. This is data without AAD:

![](images/2021-02-01-10-13-11.png)
```
ASD pval =  0
SCZ pval =  0
BD pval =  0
MDD pval =  0
AAD pval =  0.3216
IBD pval =  0
```

And without MDD:

![](images/2021-02-01-10-22-28.png)
```
ASD pval =  0 
SCZ pval =  0 
BD pval =  0 
MDD pval =  0 
AAD pval =  0 
IBD pval =  0.0015 
```

This is the Caudate without AAD:

![](images/2021-02-01-12-19-20.png)

```
ASD pval =  0 
SCZ pval =  0 
BD pval =  0 
MDD pval =  0.0097 
AAD pval =  0 
IBD pval =  0 
```

And the Caudate without MDD:

![](images/2021-02-01-12-09-14.png)
```
ASD pval =  0 
SCZ pval =  0 
BD pval =  0 
MDD pval =  0.2695 
AAD pval =  0 
IBD pval =  8e-04 
```

# 2021-02-02 11:19:56

Let's see how the TWAS significant genes look in our PM results:

```r
# bw
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_ACC.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)

wei = read.delim('~/data/expression_impute/fusion_twas-master/WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.P01.pos')
cnt = 75
wei$geneid = substring(wei$WGT, cnt, cnt+14)
mf = merge(fusion, wei[, c('ID', 'CHR', 'geneid')], by='ID', all.x=T,
            all.y=F, sort=F)
mf$adjPval = p.adjust(mf$TWAS.P, method='fdr')
twas_genes = mf[mf$adjPval < .05, 'geneid']

load('~/data/post_mortem/DGE_01272021.RData')
dge = as.data.frame(dge_acc[['protein_coding']])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
dge_twas = dge[dge$ensembl_gene_id %in% twas_genes, ]
```

```
r$> dge_twas                                                                                
                    baseMean log2FoldChange      lfcSE        stat    pvalue      padj
ENSG00000088930.8  949.43792   -0.002480057 0.03212861 -0.07719155 0.9384712 1.0000000
ENSG00000140830.9  251.99311    0.017357434 0.07711609  0.22508187 0.8219156 1.0000000
ENSG00000186792.17  97.23679   -0.044933370 0.11356677 -0.39565596 0.6923588 0.9500927
                     weight ensembl_gene_id
ENSG00000088930.8  1.000000 ENSG00000088930
ENSG00000140830.9  1.000000 ENSG00000140830
ENSG00000186792.17 3.108477 ENSG00000186792
```

That didn't work. What if we use all imputed genes?

```r
# bw
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_ACC.ALL.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)

wei = read.delim('~/data/expression_impute/fusion_twas-master/WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.pos')
cnt = 75
wei$geneid = substring(wei$WGT, cnt, cnt+14)
mf = merge(fusion, wei[, c('ID', 'CHR', 'geneid')], by='ID', all.x=T,
            all.y=F, sort=F)
mf$adjPval = p.adjust(mf$TWAS.P, method='fdr')
twas_genes = mf[mf$adjPval < .05, 'geneid']

load('~/data/post_mortem/DGE_01272021.RData')
dge = as.data.frame(dge_acc[['protein_coding']])
dge$ensembl_gene_id = substr(rownames(dge), 1, 15)
dge_twas = dge[dge$ensembl_gene_id %in% twas_genes, ]
```

I have 19 hits, and 1 is p < .05. Is that significant? Nevermind... the gene is
not that interesting: "IFITM2 (Interferon Induced Transmembrane Protein 2) is a
Protein Coding gene. Diseases associated with IFITM2 include Cycloplegia and
Influenza. Among its related pathways are Interferon gamma signaling and Innate
Immune System."

# 2021-02-04 13:41:22

Are the methylation results at all correlated with the DGE results? Nevermind...
will need to re-run the methylation analysis a bit, and adding ENSIDs to the
probes instead of hugo so I don't discartd them. Might revisit data cleaning as
well to be more in line with current practices.