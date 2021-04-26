# 2021-04-23 12:05:30

Philip was going through a list of models for testing covariates that people
have tried, and he found a paper from a good NIMH group that did this (mcMahon's
group):

```
Covariate selection: Since many factors can affect expression measurements in postmortem tissue, we employed a two-tiered strategy aimed at controlling for 1) variables for which diagnostic groups were substantially imbalanced, and 2) variables that were associated with major sources of variance in the gene expression data. We tested 32 known variables for association with diagnosis. These included RNA quality metrics (mapping percentage, GC percentage, 5-3’bias, RIN, extraction batch, library batch), available demographic data (age, gender, reported race, BMI, marital status, and manner of death), post-mortem interval (PMI), source of brain donor, tissue pH, and brain weight, along with post-mortem toxicology (the binary toxicology variables included were alcohol, opioids, major stimulants, hallucinogens, cocaine, sedative hypnotic anxiolitics, non psychiatric drugs, nicotine-cotinine, benzos, cannabinoids, antidepressants, anticholinergics, other psychotropic drugs, antipsychotics, mood stabilizers, anti-epileptics, and THC) results. In addition, 10 ancestry vectors were derived from genotyping on Illumina SNP arrays using Eigensoft [4, 5]. 

Stepwise logistic regression analysis was used to test association with the multinomial dependent variable of diagnostic group (bipolar, major depression, schizophrenia, control). Nominally significant associations were detected between diagnostic groups and race, source of brain donor, and antipsychotic exposure (Supplementary Table 3a). Of these, only antipsychotic exposure remained significant after Bonferroni-correction for the 43 covariates tested. Diagnostic groups were well-balanced in terms of the remaining variables.

Principal component (PC) analysis of the RNA counts data revealed that the top 10 PCs explained ~45% of the variance. Each of these PCs were tested for association with the known covariates by stepwise linear regression. Only race, RIN, and GC percent were significantly associated with one or more PCs at Bonferroni p<0.05 (Supplementary Table 3b). Note that the derived PCs were more strongly associated with reported race than with ancestry vectors.

On the basis of these results, downstream analyses of differential expression between diagnostic groups included race, RIN, and GC percent as covariates, while potential effects of antipsychotic exposure were tested in post-hoc analyses using approximate posterior estimation for GLM coefficients in DESeq2. In QTL analyses, which did not break down the sample by diagnosis and where ancestry-associated differences in linkage disequilibrium are important, 10 ancestry vectors were used along with RIN, race and GC percent.
```

Let's try the same thing on our data, which is somewhat similar to the approach
we were using before anyways. However, before we get into that, let's try the
minimalistic approach that's common across papers: 

```
age, gender, PMI, RINe, race (W vs other- not too fine grained),  and the BBB variable
```

```r
minimal_DGE = function(myregion, add_cov=NA) {
    data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
    rownames(data) = data[,1]
    data[,1] = NULL
    data = round(data)
    sub_name = gsub(x=colnames(data), pattern='X', replacement='')
    colnames(data) = sub_name
    # this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
    # its other replicate
    data = data[, ! colnames(data) %in% c('66552')]
    # outliers based on PCA plots
    outliers = c('68080','68096', '68108', '68084', '68082')
    data = data[, ! colnames(data) %in% outliers]

    library(gdata)
    df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
    data = data[, colnames(data) %in% df$submitted_name]
    df = df[df$submitted_name %in% colnames(data), ]
    df = df[order(df$submitted_name), ]
    data = data[, order(df$submitted_name)]

    keep_me = df$Region == myregion
    data = data[, keep_me]
    df = df[keep_me, ]
    
    # cleaning up some variables
    df$Individual = factor(df$hbcc_brain_id)
    df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
    df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
    df$MoD = factor(df$Manner.of.Death)
    df$Sex = factor(df$Sex)
    df$batch = factor(df$batch)
    df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
    df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
    df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
    df$substance_group = factor(df$substance_group)
    df$comorbid_group = factor(df$comorbid_group_update)
    df$evidence_level = factor(df$evidence_level)
    df$brainbank = factor(df$bainbank)
    # replace the one subject missing population PCs by the median of their
    # self-declared race and ethnicity
    idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
    pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
    med_pop = apply(df[idx, pop_pcs], 2, median)
    df[which(is.na(df$C1)), pop_pcs] = med_pop
    df$BBB = factor(sapply(1:nrow(df),
                            function(x) sprintf('%s_%s',
                                        as.character(df[x,'brainbank']),
                                        as.character(df[x, 'batch']))))
    df$BBB2 = NA                                                                        
    df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
    df[df$batch==3, 'BBB2'] = 2                                                         
    df[df$batch==4, 'BBB2'] = 3      
    df$BBB2 = factor(df$BBB2)
    imWNH = which(df$C1 > 0 & df$C2 < -.075)
    df$POP_BIN = 'other'
    df[imWNH, 'POP_BIN'] = 'WNH'
    df$POP_BIN = factor(df$POP_BIN)        
    # df$RINc = cut(df$RINe, breaks = 4)    
    breaks = quantile(df$RINe, probs = seq(0, 1, by = 0.25))
    df$RINc = cut(df$RINe, breaks=breaks, labels=c('q1', 'q2', 'q3', 'q4'),
                  include.lowest=T)
    print(summary(df$RINc))

    library(GenomicFeatures)
    txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
    txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
                "GENEID")
    bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
    bt_slim = bt[, c('gene_id', 'gene_biotype')]
    bt_slim = bt_slim[!duplicated(bt_slim),]
    txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
    tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
    tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
    imautosome = which(tx_meta$TXCHROM != 'X' &
                    tx_meta$TXCHROM != 'Y' &
                    tx_meta$TXCHROM != 'MT')
    data = data[imautosome, ]
    tx_meta = tx_meta[imautosome, ]

    library("DESeq2")
    if (is.na(add_cov)) {
        fm_str = '~ RINc + BBB2 + Diagnosis'
    } else {
        cov_str = paste0(add_cov, collapse = ' + ')
        fm_str = sprintf('~ %s + RINc + BBB2 + Diagnosis',
                         cov_str)
    }
    # making sure any numeric covariates are scaled
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
            'C1', 'C2', 'C3', 'C4', 'C5')
    for (var in num_vars) {
        df[, var] = scale(df[, var])
    }

    cat('Running', fm_str, '\n')
    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = df,
                                  design = as.formula(fm_str))
    dds = DESeq(dds)

    library(edgeR)
    design = model.matrix(as.formula(fm_str), data=colData(dds))
    isexpr <- filterByExpr(counts(dds), design=design)
    ddsExpr = dds[isexpr, ]
    ddsExpr = DESeq(ddsExpr)

    return(ddsExpr)
}
```

```r
dds.ACC = minimal_DGE('ACC')
dds.Caudate = minimal_DGE('Caudate')
save(dds.ACC, dds.Caudate, file='~/data/post_mortem/minimal_DGE_04232021.RData')
```

OK, let's see what we got:

```
r$> res = results(dds, name = "Diagnosis_Case_vs_Control", alpha=.05)                   

r$> summary(res)                                                                        

out of 25091 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 2, 0.008%
LFC < 0 (down)     : 1, 0.004%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


r$> res = results(dds, name = "Diagnosis_Case_vs_Control", alpha=.1)                    

r$> summary(res)                                                                        

out of 25091 with nonzero total read count
adjusted p-value < 0.1
LFC > 0 (up)       : 5, 0.02%
LFC < 0 (down)     : 1, 0.004%
outliers [1]       : 0, 0%
low counts [2]     : 0, 0%
(mean count < 2)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results

r$> library(IHW)                                                                        

r$> padj.IHW = adj_pvalues(ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha=0.05))
                                                                                        
r$> sum(padj.IHW < .05)                                                                 
[1] 3

r$> padj.IHW = adj_pvalues(ihw(pvalue ~ baseMean,  data=as.data.frame(res), alpha=0.1)) 

r$> sum(padj.IHW < .1)                                                                  
[1] 8
```

Could these be outliers? Yes, the top one was much larger l2FC than everyone
else. Let's try the iterative filtering...

## PCA approach

For the PCA approach, it'll be somewhat similar to what we were doing in the
past:

```r
myregion = 'ACC'

data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name
# this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# its other replicate
data = data[, ! colnames(data) %in% c('66552')]
# outliers based on PCA plots
outliers = c('68080','68096', '68108', '68084', '68082')
data = data[, ! colnames(data) %in% outliers]

library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]

keep_me = df$Region == myregion
data = data[, keep_me]
df = df[keep_me, ]

# cleaning up some variables
df$Individual = factor(df$hbcc_brain_id)
df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
df$MoD = factor(df$Manner.of.Death)
df$Sex = factor(df$Sex)
df$batch = factor(df$batch)
df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
df$substance_group = factor(df$substance_group)
df$comorbid_group = factor(df$comorbid_group_update)
df$evidence_level = factor(df$evidence_level)
df$brainbank = factor(df$bainbank)
# replace the one subject missing population PCs by the median of their
# self-declared race and ethnicity
idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
med_pop = apply(df[idx, pop_pcs], 2, median)
df[which(is.na(df$C1)), pop_pcs] = med_pop
df$BBB = factor(sapply(1:nrow(df),
                        function(x) sprintf('%s_%s',
                                    as.character(df[x,'brainbank']),
                                    as.character(df[x, 'batch']))))
df$BBB2 = NA                                                                        
df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
df[df$batch==3, 'BBB2'] = 2                                                         
df[df$batch==4, 'BBB2'] = 3      
df$BBB2 = factor(df$BBB2)
imWNH = which(df$C1 > 0 & df$C2 < -.075)
df$POP_BIN = 'other'
df[imWNH, 'POP_BIN'] = 'WNH'
df$POP_BIN = factor(df$POP_BIN)        
# df$RINc = cut(df$RINe, breaks = 4)  
# bining so DESeq2 can do its own filyering automatically
breaks = quantile(df$RINe, probs = seq(0, 1, by = 0.25))
df$RINc = cut(df$RINe, breaks=breaks, labels=c('q1', 'q2', 'q3', 'q4'),
              include.lowest=T)

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
            "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                tx_meta$TXCHROM != 'Y' &
                tx_meta$TXCHROM != 'MT')
data = data[imautosome, ]
tx_meta = tx_meta[imautosome, ]

# remove constant genes (including zeros) as it breaks PCA
const_genes = apply(data, 1, sd) == 0
data = data[!const_genes, ]

# checking which PCs are associated with our potential nuiscance variables
set.seed(42)
mypca <- prcomp(t(data), scale=TRUE)
# how many PCs to keep... using Kaiser thredhold, close to eigenvalues < 1
library(nFactors)
eigs <- mypca$sdev^2
nS = nScree(x=eigs)
keep_me = seq(1, nS$Components$nkaiser)

mydata = data.frame(mypca$x[, keep_me])
# create main metadata data frame including metadata and PCs
data.pm = cbind(df, mydata)
rownames(data.pm) = df$hbcc_brain_id
cat('Using', nS$Components$nkaiser, 'PCs from possible', ncol(data), '\n')

# check which PCs are associated at nominal p<.01
num_vars = c('Age', 'PMI',
            'C1', 'C2', 'C3', 'C4', 'C5')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                    dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data.pm[, x], data.pm[, y], method='spearman')
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

categ_vars = c('BBB2', 'MoD', 'substance_group', 'RINc',
            'comorbid_group', 'Sex', 'evidence_level')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                        dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(data.pm[, y], data.pm[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}
mypvals = rbind(categ_pvals, num_pvals)
print(which(mypvals < .01/(ncol(mypvals)*nrow(mypvals)), arr.ind = T) )
print(which(mypvals < .05/(ncol(mypvals)*nrow(mypvals)), arr.ind = T) )
```

For ACC we have 7 PCs:

```
     row col
BBB2   1   1
RINc   4   2
BBB2   1   3
     row col
BBB2   1   1
RINc   4   2
BBB2   1   3
```

and for Caudate we have 8 PCs:

```
     row col
BBB2   1   1
BBB2   1   2
     row col
BBB2   1   1
BBB2   1   2
```

So, let's run a function:

```r
pca_DGE = function(myregion, fm_str) {
    data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
    rownames(data) = data[,1]
    data[,1] = NULL
    data = round(data)
    sub_name = gsub(x=colnames(data), pattern='X', replacement='')
    colnames(data) = sub_name
    # this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
    # its other replicate
    data = data[, ! colnames(data) %in% c('66552')]
    # outliers based on PCA plots
    outliers = c('68080','68096', '68108', '68084', '68082')
    data = data[, ! colnames(data) %in% outliers]

    library(gdata)
    df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
    data = data[, colnames(data) %in% df$submitted_name]
    df = df[df$submitted_name %in% colnames(data), ]
    df = df[order(df$submitted_name), ]
    data = data[, order(df$submitted_name)]

    keep_me = df$Region == myregion
    data = data[, keep_me]
    df = df[keep_me, ]

    # cleaning up some variables
    df$Individual = factor(df$hbcc_brain_id)
    df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
    df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
    df$MoD = factor(df$Manner.of.Death)
    df$Sex = factor(df$Sex)
    df$batch = factor(df$batch)
    df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
    df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
    df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
    df$substance_group = factor(df$substance_group)
    df$comorbid_group = factor(df$comorbid_group_update)
    df$evidence_level = factor(df$evidence_level)
    df$brainbank = factor(df$bainbank)
    # replace the one subject missing population PCs by the median of their
    # self-declared race and ethnicity
    idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
    pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
    med_pop = apply(df[idx, pop_pcs], 2, median)
    df[which(is.na(df$C1)), pop_pcs] = med_pop
    df$BBB = factor(sapply(1:nrow(df),
                            function(x) sprintf('%s_%s',
                                        as.character(df[x,'brainbank']),
                                        as.character(df[x, 'batch']))))
    df$BBB2 = NA                                                                        
    df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
    df[df$batch==3, 'BBB2'] = 2                                                         
    df[df$batch==4, 'BBB2'] = 3      
    df$BBB2 = factor(df$BBB2)
    imWNH = which(df$C1 > 0 & df$C2 < -.075)
    df$POP_BIN = 'other'
    df[imWNH, 'POP_BIN'] = 'WNH'
    df$POP_BIN = factor(df$POP_BIN)        
    # df$RINc = cut(df$RINe, breaks = 4)  
    # bining so DESeq2 can do its own filyering automatically
    breaks = quantile(df$RINe, probs = seq(0, 1, by = 0.25))
    df$RINc = cut(df$RINe, breaks=breaks, labels=c('q1', 'q2', 'q3', 'q4'),
                include.lowest=T)

    library(GenomicFeatures)
    txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
    txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
                "GENEID")
    bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
    bt_slim = bt[, c('gene_id', 'gene_biotype')]
    bt_slim = bt_slim[!duplicated(bt_slim),]
    txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
    tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
    tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
    imautosome = which(tx_meta$TXCHROM != 'X' &
                    tx_meta$TXCHROM != 'Y' &
                    tx_meta$TXCHROM != 'MT')
    data = data[imautosome, ]
    tx_meta = tx_meta[imautosome, ]

    # remove constant genes (including zeros) as it breaks PCA
    const_genes = apply(data, 1, sd) == 0
    data = data[!const_genes, ]

    library("DESeq2")
    # making sure any numeric covariates are scaled
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
            'C1', 'C2', 'C3', 'C4', 'C5')
    for (var in num_vars) {
        df[, var] = scale(df[, var])
    }

    cat('Running', fm_str, '\n')
    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = df,
                                  design = as.formula(fm_str))
    dds <- dds[keep,]
    dds = DESeq(dds)

    library(edgeR)
    design = model.matrix(as.formula(fm_str), data=colData(dds))
    isexpr <- filterByExpr(counts(dds), design=design)
    ddsExpr = dds[isexpr, ]
    ddsExpr = DESeq(ddsExpr)

    return(ddsExpr)
}
```

```r
dds.ACC = pca_DGE('ACC', '~ RINc + BBB2 + Diagnosis')
dds.Caudate = pca_DGE('Caudate', '~ BBB2 + Diagnosis')
save(dds.ACC, dds.Caudate, file='~/data/post_mortem/pca_DGE_04262021.RData')
```

Let's check out our results then:

```
r$> library(IHW)                                                                        

r$> summary(results(dds.ACC, name = "Diagnosis_Case_vs_Control", alpha=.05))            

out of 24387 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 60, 0.25%
LFC < 0 (down)     : 5, 0.021%
outliers [1]       : 13, 0.053%
low counts [2]     : 946, 3.9%
(mean count < 7)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


r$> summary(results(dds.ACC, name = "Diagnosis_Case_vs_Control", alpha=.05, filterFun=ih
    w))                                                                                 

out of 24387 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 71, 0.29%
LFC < 0 (down)     : 6, 0.025%
outliers [1]       : 13, 0.053%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting

r$> summary(results(dds.Caudate, name = "Diagnosis_Case_vs_Control", alpha=.05))        

out of 24059 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 0, 0%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 23, 0.096%
low counts [2]     : 0, 0%
(mean count < 3)
[1] see 'cooksCutoff' argument of ?results
[2] see 'independentFiltering' argument of ?results


r$> summary(results(dds.Caudate, name = "Diagnosis_Case_vs_Control", alpha=.05, filterFu
    n=ihw))                                                                             

out of 24059 with nonzero total read count
adjusted p-value < 0.05
LFC > 0 (up)       : 0, 0%
LFC < 0 (down)     : 0, 0%
outliers [1]       : 23, 0.096%
[1] see 'cooksCutoff' argument of ?results
see metadata(res)$ihwResult on hypothesis weighting
```

A bit too much? Outliers?






vsd <- vst(dds, blind=FALSE)
norm.cts <- assay(vsd)
mat <- limma::removeBatchEffect(norm.cts, vsd$BBB2)

gnames = data.frame(full=rownames(counts(dds)),
                    nov=substr(rownames(counts(dds)), 1, 15))
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
gnames = merge(gnames, mart, by.x='nov', by.y='ensembl_gene_id')
keep_me = gnames$nov %in% gene_list
gene_ids = gnames[keep_me, ]

resid_expr = reshape2::melt(mat[gene_ids$full,])
colnames(resid_expr) = c('gene', 'submitted_name', 'normCount')
junk = colData(vsd)[, c('Diagnosis', 'submitted_name')]
resid_expr = merge(resid_expr, junk, by='submitted_name')
resid_expr = merge(resid_expr, gene_ids, by.x='gene', by.y='full')

library(ggpubr)
library(ggbeeswarm)
quartz()
myplots = list()
clrs = c("green3", "red")
for (g in 1:10) {
    cat(gene_ids[g, 'nov'], '\n')
    d = as.data.frame(resid_expr[resid_expr$gene == gene_list[g],])
    p = (ggplot(d, aes(x=Diagnosis, y=normCount, color = Diagnosis,
                    fill = Diagnosis)) + 
        scale_y_log10() +
        geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8,
                    lwd = 0.5) +
        stat_summary(fun = mean, geom = "point", color = "black",
                    shape = 5, size = 3,
                    position=position_dodge(width = 0.8)) +
        scale_color_manual(values = clrs) +
        scale_fill_manual(values = clrs) +
        geom_quasirandom(color = "black", size = 1, dodge.width = 0.8) +
        theme_bw() + #theme(legend.position = "none") + 
        ggtitle(gene_ids[g, 'hgnc_symbol']))
    myplots[[g]] = p
}
p = ggarrange(plotlist=myplots)
print(p)

```

# TODO
 * does DESeq2 have recommendations for binarizing?
 * figure out best way to bin RINe, and also when to do it. I think we need to
   split by region before creating first dds model to estimate everything
   properly...
 * re-run all results for ACC (caudate seems fine)
