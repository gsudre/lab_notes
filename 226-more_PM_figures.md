# 2021-05-06 05:53:44

Let's make a few more figures.

## PCA details

The goal is to make within-region scree plots, and then a matrix of
-log10P*sign(stat) which will include DX.

```r
# not going to worry about people we don't have data for, or are outliers
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
df$SUB2 = 'no'
df[df$substance_group > 0, 'SUB2'] = 'yes'
df$SUB2 = factor(df$SUB2)
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

df2 = df[!duplicated(df$hbcc_brain_id), ]

# run nonparametric t-tests for numeric variables
num_vars = c('Age', 'PMI', 'C1', 'C2', 'C3', 'C4', 'C5', 'RINe')
mypvals = c()
mystats = c()
for (x in num_vars) {
    res = wilcox.test(as.formula(sprintf('%s ~ Diagnosis', x)), data=df2)
    mypvals = c(mypvals, res$p.value)
    mystats = c(mystats, res$statistic)
}

categ_vars = c('MoD', 'SUB2', 'comorbid_group', 'Sex', 'evidence_level', 'BBB2')
for (x in categ_vars) {
    res = chisq.test(table(df2$Diagnosis, df2[, x]))
    mypvals = c(mypvals, res$p.value)
    mystats = c(mystats, res$statistic)
}
print(c(num_vars, categ_vars)[which(mypvals < .05/length(mypvals))])
```

```
[1] "SUB2"
[1] "C1"             "SUB2"           "comorbid_group"
```

```r
myvars = c(num_vars, categ_vars)
DX_pvals = mypvals
DX_plot = -log10(mypvals) * sign(mystats)
names(DX_plot) = myvars
```

Now we run it for the PCA of ACC. Make sure we save the data for the scree plot
as well:

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
df$SUB2 = 'no'
df[df$substance_group > 0, 'SUB2'] = 'yes'
df$SUB2 = factor(df$SUB2)
df$substance_group = factor(df$substance_group)
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
for (var in num_vars) {
    df[, var] = scale(df[, var])
}

min_subjs = min(table(df$Diagnosis))
keep <- rowSums(data == 0) <= min_subjs
data <- data[keep,]

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
mycorrs = rbind(categ_corrs, num_corrs)
print(which(mypvals < .01/(ncol(mypvals)*nrow(mypvals)), arr.ind = T) )
print(which(mypvals < .05/(ncol(mypvals)*nrow(mypvals)), arr.ind = T) )
```

```
# ACC
Using 7 PCs from possible 53 
There were 21 warnings (use warnings() to see them)
     row col
BBB2   6   1
RINe  14   2
BBB2   6   3
     row col
BBB2   6   1
RINe  14   2
BBB2   6   3

# Caudate
Using 8 PCs from possible 56 
There were 50 or more warnings (use warnings() to see the first 50)
     row col
BBB2   6   1
RINe  14   2
     row col
BBB2   6   1
BBB2   6   2
RINe  14   2
BBB2   6   6
```

```r
ACC_pvals = mypvals
ACC_plot = -log10(mypvals) * sign(mycorrs)
names(ACC_plot) = myvars
ACC_pca = mypca
```

Now, let's work on plotting them. Scree first:

```r
library(ggplot2)
quartz()
pca_obj = ACC_pca
npcs = 12
var_explained_df <- data.frame(PC=factor(colnames(pca_obj$x)[1:npcs],
                                         levels=colnames(pca_obj$x)[1:npcs]),
                               var_explained=100*(pca_obj$sdev[1:npcs])^2/sum((pca_obj$sdev)^2))
ggplot(aes(x=PC, y=var_explained), data=var_explained_df)+
  geom_col()+
  geom_hline(yintercept=var_explained_df[7, 2], linetype="dashed", color = "red") +
  labs(title="Scree plot: ACC RNAseq counts",
        y='Percent variance explained (%)',
        x='Principal components') +
   ylim(0, 50)

pca_obj = mypca
var_explained_df <- data.frame(PC=factor(colnames(pca_obj$x)[1:npcs],
                                         levels=colnames(pca_obj$x)[1:npcs]),
                               var_explained=100*(pca_obj$sdev[1:npcs])^2/sum((pca_obj$sdev)^2))
ggplot(aes(x=PC, y=var_explained), data=var_explained_df)+
  geom_col()+
  geom_hline(yintercept=var_explained_df[8, 2], linetype="dashed", color = "red") +
  labs(title="Scree plot: Caudate RNAseq counts",
        y='Percent variance explained (%)',
        x='Principal components') +
   ylim(0, 50)

```

![](images/2021-05-06-06-43-22.png)

![](images/2021-05-06-06-43-02.png)

Now the correlation plot:

```r
cnames = sapply(colnames(ACC_plot), function(x) sprintf('ACC: %s', x))
plot_df = t(ACC_plot)
rownames(plot_df) = cnames
Caudate_plot = -log10(mypvals) * sign(mycorrs)
cnames = sapply(colnames(Caudate_plot), function(x) sprintf('Caudate: %s', x))
colnames(Caudate_plot) = cnames
plot_df = rbind(plot_df, t(Caudate_plot))
junk = t(data.frame(DX_plot))
plot_df2 = rbind(plot_df, junk[, match(colnames(plot_df), colnames(junk))])
rownames(plot_df2)[nrow(plot_df2)] = 'Diagnosis'
mylim = max(abs(plot_df2))

colnames(plot_df2) = c("Mode of death", "Substance abuse", "Comorbidities",
                       "Sex", "Evidence level", "Brain bank / batch", "Age",
                       "Post-mortem interval", "C1", "C2", "C3", "C4", "C5",
                       "RINe")
library(corrplot)
corrplot(plot_df2, is.corr=F, cl.lim=c(-mylim, mylim), tl.col='black')
```

![](images/2021-05-06-07-11-56.png)

## Expression plots

Let's focus on the expression plots now, and I'll see if I can find something to
move focus away from anything people might think are outliers (even though they
are not):

```r
plot_expression_noDots = function(gene_ids, dds, t_str) {
    library(ggpubr)
    library(ggbeeswarm)
    quartz()
    mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
    myplots = list()
    clrs = c("green3", "red")
    for (g in 1:length(gene_ids)) {
        GENEID = substr(gene_ids[g], 1, 15)
        hgnc = mart[mart$ensembl_gene_id == GENEID, 'hgnc_symbol']
        if (hgnc == '') {
            hgnc = GENEID
        }
        cat(gene_ids[g], hgnc, '\n')
        d <- plotCounts(dds, gene=gene_ids[g], intgroup="Diagnosis",
                        returnData=TRUE)
        p = (ggplot(d, aes(x=Diagnosis, y=count, color = Diagnosis,
                        fill = Diagnosis)) + 
            scale_y_log10() +
            geom_boxplot(alpha = 0.4, outlier.shape = '+', width = 0.8,
                        lwd = 0.5, notch=T) +
            scale_color_manual(values = clrs, labels=c('Unaffected', 'ADHD')) +
            scale_fill_manual(values = clrs, labels=c('Unaffected', 'ADHD')) +
            theme_bw() + theme(axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_blank()) +
            ggtitle(hgnc))
        # trimming out the outliers so they don't influence the y axis
        ylim1 = boxplot.stats(d[d$Diagnosis=='Case','count'])$stats[c(1, 5)]
        ylims = quantile(d$count, c(0, 0.99))
        myplots[[g]] = p + 
                       theme(legend.position = "none") # +
                    #    coord_cartesian(ylim = ylims)
    }
    myplots[[g + 1]] = get_legend(p)
    p = ggarrange(plotlist=myplots)
    print(annotate_figure(p, t_str))
}
load('~/data/post_mortem/pca_DGE_bigger_04292021.RData')
res = results(dds.ACC, name = "Diagnosis_Case_vs_Control", alpha=.05)
plot_expression_noDots(rownames(res)[which(res$padj < .05)], dds.ACC, 'ACC')
res = results(dds.Caudate, name = "Diagnosis_Case_vs_Control", alpha=.05)
plot_expression_noDots(rownames(res)[which(res$padj < .05)], dds.Caudate, 'Caudate')
```

![](images/2021-05-06-12-00-55.png)

![](images/2021-05-06-12-03-43.png)

# 2021-05-07 12:07:07

Let's see if I normalize them differently I can reduce the weight of some of
those outliers:

```r
gene_list = c('HILPDA', 'MYO1G')
r = 'Caudate'

load('~/data/post_mortem/pca_DGE_bigger_04292021.RData')
res_str = sprintf('dds = dds.%s', r)
eval(parse(text=res_str))

res = results(dds, name = "Diagnosis_Case_vs_Control", alpha=.05)
vsd <- rlog(dds, blind=TRUE)
norm.cts <- assay(vsd)

covars = model.matrix(~RINe + C1 + BBB2 + comorbid_group + SUB2,
                      data=colData(dds))
dsn = model.matrix(~Diagnosis, data=colData(dds))
mat <- limma::removeBatchEffect(norm.cts, covariates=covars, design=dsn)
# mat <- limma::removeBatchEffect(norm.cts, covariates=NULL, design=NULL)

gnames = data.frame(full=rownames(counts(dds)),
                    nov=substr(rownames(counts(dds)), 1, 15))
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
gnames = merge(gnames, mart, by.x='nov', by.y='ensembl_gene_id')
keep_me = gnames$hgnc_symbol %in% gene_list
gene_ids = gnames[keep_me, ]

resid_expr = reshape2::melt(mat[gene_ids$full,])
colnames(resid_expr) = c('gene', 'submitted_name', 'normCount')
junk = colData(vsd)[, c('Diagnosis', 'submitted_name')]
resid_expr = merge(resid_expr, junk, by='submitted_name')
resid_expr = merge(resid_expr, gene_ids, by.x='gene', by.y='full')

# plotting each of the significant genes
library(ggpubr)
library(ggbeeswarm)
myplots = list()
clrs = c("green3", "red")
for (g in 1:nrow(gene_ids)) {
    cat(gene_ids[g, 'nov'], '\n')
    d = as.data.frame(resid_expr[resid_expr$hgnc_symbol == gene_list[g],])
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

![](images/2021-05-07-12-32-25.png)

Those two use vsd with blind=F on the left, and TRUE on right. Not much
difference, but looks like the boxplots are a bit more well-behaved. Not much
difference between that and rlog (below), excpet that vst is much faster.

![](images/2021-05-07-12-37-30.png)

Let's then redo the figures using the normalized data:

```r
r = 'ACC'
r='Caudate'
load('~/data/post_mortem/pca_DGE_bigger_04292021.RData')

res_str = sprintf('dds = dds.%s', r)
eval(parse(text=res_str))
res = results(dds, name = "Diagnosis_Case_vs_Control", alpha=.05)
gene_list = rownames(res)[which(res$padj <= .06)]

vsd <- vst(dds, blind=FALSE)
norm.cts <- assay(vsd)

covars = model.matrix(~ RINe + C1 + BBB2 + comorbid_group + SUB2,
                      data=colData(dds))
dsn = model.matrix(~ Diagnosis, data=colData(dds))
mat <- limma::removeBatchEffect(norm.cts, covariates=covars, design=dsn)

gnames = data.frame(full=rownames(counts(dds)),
                    nov=substr(rownames(counts(dds)), 1, 15))
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
gnames = merge(gnames, mart, by.x='nov', by.y='ensembl_gene_id')
keep_me = gnames$full %in% gene_list
gene_ids = gnames[keep_me, ]

resid_expr = reshape2::melt(mat[gene_ids$full,])
colnames(resid_expr) = c('gene', 'submitted_name', 'normCount')
junk = colData(vsd)[, c('Diagnosis', 'submitted_name')]
resid_expr = merge(resid_expr, junk, by='submitted_name')
resid_expr = merge(resid_expr, gene_ids, by.x='gene', by.y='full')

# plotting each of the significant genes
library(ggpubr)
library(ggbeeswarm)
myplots = list()
clrs = c("green3", "red")
for (g in 1:nrow(gene_ids)) {
    hgnc = gene_ids[g, 'hgnc_symbol']
    if (hgnc == '') {
        hgnc = gene_ids[g, 'nov']
    }
    cat(gene_ids[g, 'nov'], hgnc, '\n')
    d = as.data.frame(resid_expr[resid_expr$gene == gene_list[g],])
    p = (ggplot(d, aes(x=Diagnosis, y=normCount, color = Diagnosis,
                        fill = Diagnosis)) + 
            scale_y_log10() +
            geom_boxplot(alpha = 0.4, outlier.shape = '+', width = 0.8,
                        lwd = 0.5, notch=T) +
            scale_color_manual(values = clrs, labels=c('Unaffected', 'ADHD')) +
            scale_fill_manual(values = clrs, labels=c('Unaffected', 'ADHD')) +
            theme_bw() + theme(axis.text.x = element_blank(),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.y = element_blank(),) +
            ggtitle(hgnc))
    myplots[[g]] = p + theme(legend.position = "none")
}
myplots[[g + 1]] = get_legend(p)
p = ggarrange(plotlist=myplots)
print(annotate_figure(p, r))
```

![](images/2021-05-07-12-51-27.png)

![](images/2021-05-07-12-56-17.png)

I plotted 3 because one was screwing up the dimensions and the functions weren't
working. Just crop it.

# 2021-05-10 10:44:57

Let's redo the correlation picture, and make the non-significant results very
dim. I'm using this list from note 222:

```
ASD_Gandal_micro_ACC 0 
SCZ_Gandal_micro_ACC 0 
BD_Gandal_micro_ACC 0 
MDD_Gandal_micro_ACC 0.2606 
AAD_Gandal_micro_ACC 0.2823 
IBD_Gandal_micro_ACC 0 
ASD_Gandal_micro_Caudate 0 
SCZ_Gandal_micro_Caudate 0 
BD_Gandal_micro_Caudate 0 
MDD_Gandal_micro_Caudate 0 
AAD_Gandal_micro_Caudate 0 
IBD_Gandal_micro_Caudate 0.1398 
SCZ_Gandal_RNAseq_ACC 0 
BD_Gandal_RNAseq_ACC 0.0154 
SCZ_Gandal_RNAseq_Caudate 0 
BD_Gandal_RNAseq_Caudate 0 
ASD_Gandal_RNAseq_ACC 0 
ASD_Gandal_RNAseq_Caudate 0 
BD_Akula_ACC 2e-04 
SCZ_Akula_ACC 0 
MDD_Akula_ACC 7e-04 
SCZ_Benjamin_Caudate 0 
BD_Pacifico_Caudate 0.294 
OCD_Piantadosi_Caudate 0 
ASD_Wright_DLPFC_ACC 0 
ASD_Neelroop_FrontalTemporal_ACC 0 
```

```r
library(ggplot2)
quartz()

fname = 'disorders_corrs_bigger_04292021'
corrs = readRDS(sprintf('~/data/post_mortem/%s.rds', fname))
# quick hack to make IBD last, which affects the X labels but not the color
corrs[corrs$disorder == 'IBD', 'disorder'] = 'zIBD'

# this leveling only affects the color ordering
corrs$disorder = factor(corrs$disorder,
                        levels=c('AAD', 'ASD', 'BD', 'MDD', 'OCD', 'SCZ', 'zIBD'))
col_labels = c('Alcohol Abuse or Dependene', 'Autism Spectrum Disorder',
               'Bipolar Disorder', 'Major Depression Disorder',
               'Obsessive Compulsive Disorder', 'Schizophernia',
               'Irritable Bowel Disorder')

# just to share axis
ymax = .4 #max(corrs$corr)
ymin = -.4 #min(corrs$corr)
my_colors = RColorBrewer::brewer.pal(7, "Accent")
imbad = c('MDD_Gandal_micro', 'AAD_Gandal_micro', 'BD_Gandal_RNAseq')

r = 'ACC' 
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s_%s',
                                      mycorrs[i, 'disorder'],
                                      mycorrs[i, 'source']))
mycorrs$alpha = 1
mycorrs$colour = 'black'
mycorrs[mycorrs$id %in% imbad, 'alpha'] = .9
mycorrs[mycorrs$id %in% imbad, 'colour'] = 'grey'
mycorrs$colour = factor(mycorrs$colour)

my_labels = c('AAD [1]',
              'ASD [1]',
              'ASD [2]',
              'ASD [3]',
              'ASD [4]',
              'BD [5]',
              'BD [1]',
              'BD [2]',
              'MDD [5]',
              'MDD [1]',
              'SCZ [5]',
              'SCZ [1]',
              'SCZ [2]',
              'IBS [1]'
              )
myrefs = c('[1] Gandal et al. 2018 (microarray)',
           '[2] Gandal et al. 2018 (RNAseq)',
           '[3] Parikshak et al. 2016',
           '[4] Wright et al. 2017',
           '[5] Akula et al. 2020',
           '[6] Pacifico and Davis, 2017',
           '[7] Piantadosi et al. 2021',
           '[8] Benjamin et al. 2020')
p <- ggplot(mycorrs, aes(x = factor(id), y = corr, fill=disorder,
                         alpha=alpha)) +
    geom_violin(trim=FALSE, aes(colour=colour)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          axis.title.x = element_blank())
p1 = p + ggtitle(r) + geom_hline(yintercept=0, linetype="dotted",
                                color = "red", size=1) +
   ylab('Transcriptome correlation (rho)') + ylim(ymin, ymax) + 
   scale_fill_manual(breaks = levels(corrs$disorder),
                     values = my_colors,
                     labels = col_labels) +
    scale_colour_manual(values = c('black'='#000000', 'grey'='#808080')) +
   scale_x_discrete(labels=my_labels) + 
   guides(colour=F, alpha=F)


r = 'Caudate'
imbad = c('BD_Pacifico', 'zIBD_Gandal_micro')
mycorrs = corrs[corrs$region == r, ]
mycorrs$id = sapply(1:nrow(mycorrs),
                  function(i) sprintf('%s_%s',
                                      mycorrs[i, 'disorder'],
                                      mycorrs[i, 'source']))
mycorrs$alpha = 1
mycorrs$colour = 'black'
mycorrs[mycorrs$id %in% imbad, 'alpha'] = .9
mycorrs[mycorrs$id %in% imbad, 'colour'] = 'grey'
mycorrs$colour = factor(mycorrs$colour)
my_labels = c('AAD [1]',
              'ASD [1]',
              'ASD [2]',
              'BD [1]',
              'BD [2]',
              'BD [6]',
              'MDD [1]',
              'OCD [7]',
              'SCZ [8]',
              'SCZ [1]',
              'SCZ [2]',
              'IBD [1]'
              )
p <- ggplot(mycorrs, aes(x = factor(id), y = corr, fill=disorder,
                         alpha=alpha)) +
    geom_violin(trim=FALSE, aes(colour=colour)) +
    theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
          axis.title.x = element_blank())
p2 = p + ggtitle(r) + geom_hline(yintercept=0, linetype="dotted",
                                color = "red", size=1) +
   ylim(ymin, ymax) + 
   scale_fill_manual(breaks = levels(corrs$disorder),
                     values = my_colors,
                     labels = col_labels,
                     name='Disorders') +
    theme(axis.title.y = element_blank()) + 
   scale_colour_manual(values = c('black'='#000000', 'grey'='#808080')) +
   scale_x_discrete(labels=my_labels) +
   guides(colour=F, alpha=F)

ggarrange(p1, p2, common.legend = T, legend='right',
          legend.grob=get_legend(p2)) 
```

![](images/2021-05-10-15-03-47.png)


## Redoing PCA plots

```r
data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name
# this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# its other replicate
data = data[, ! colnames(data) %in% c('66552')]
# # outliers based on PCA plots
# outliers = c('68080','68096', '68108', '68084', '68082')
# data = data[, ! colnames(data) %in% outliers]

library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]

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
df$SUB2 = 'no'
df[df$substance_group > 0, 'SUB2'] = 'yes'
df$SUB2 = factor(df$SUB2)
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

num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
            'C1', 'C2', 'C3', 'C4', 'C5')
for (var in num_vars) {
    df[, var] = scale(df[, var])
}

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

min_subjs = nrow(df) - 1 
keep <- rowSums(data == 0) <= min_subjs
data <- data[keep,]

library("DESeq2")
fm_str = '~ RINe + BBB2 + comorbid_group + SUB2 + Diagnosis'
dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = df,
                                design = as.formula(fm_str))

vsd <- vst(dds, blind=FALSE)
norm.cts <- assay(vsd)

mysds = apply(norm.cts, 1, sd)
smysds = sort(mysds, decreasing=T, index.return=T)

set.seed(42)
mypca <- prcomp(t(norm.cts[smysds$ix[1:300], ]), scale=TRUE)
pca_obj = mypca

library(ggplot2)
library(ggpubr)
quartz()
pca_obj = mypca
npcs = 12
var_explained_df <- data.frame(PC=factor(colnames(pca_obj$x)[1:npcs],
                                         levels=colnames(pca_obj$x)[1:npcs]),
                               var_explained=100*(pca_obj$sdev[1:npcs])^2/sum((pca_obj$sdev)^2))
p1 = ggplot(aes(x=PC, y=var_explained), data=var_explained_df)+
  geom_col()+
  labs(title="Scree plot: RNAseq counts",
        y='Percent variance explained (%)',
        x='Principal components') +
   ylim(0, 75)
p1

plot_df = data.frame(pca_obj$x)
p2 = ggplot(aes(x=PC1, y=PC2, color=df$Region), data=plot_df) + geom_point() + 
     scale_color_discrete(name='Region', labels=levels(df$Region))
p3 = ggplot(aes(x=PC2, y=PC3, color=df$BBB2), data=plot_df) + geom_point() + 
     scale_color_discrete(name='Batches (Brain bank)',
                          labels=c('1 and 2 (NIMH)', '3 (Pitt and UMBN)',
                                   '4 (Pitt and UMBN)'))
ggarrange(p2, p3, nrow=2, ncol=1, legend='bottom') 
```

![](images/2021-05-11-14-03-58.png)

![](images/2021-05-11-14-04-05.png)


## Getting the correct numbers

```r
data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name

df = read.csv('~/data/rnaseq_derek/UPDATED_file_for_derek_add_cause_of_death.csv')
df = df[!duplicated(df$submitted_name),]
```

```
r$> dim(df)                                                                             
[1] 120  28

r$> table(df$Region)                                                                    

    ACC Caudate 
     56      64 

r$> table(df$Region, df$hbcc_brain_id)                                                  
         
          887 890 993 1227 1292 1331 1435 1524 1530 1683 1709 1908 2080 2371 2485 2510
  ACC       1   1   1    1    1    1    1    1    1    0    1    1    1    1    0    1
  Caudate   1   1   1    1    1    1    1    1    1    1    1    1    1    1    1    1
         
          2546 2623 2636 2646 2724 2753 2781 2842 2873 2877 3002 3003 3005 3006 3007
  ACC        1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
  Caudate    1    1    1    1    1    1    1    1    1    6    1    1    1    0    1
         
          3008 3009 3010 3011 3012 3013 3014 3016 3017 3018 3019 3020 3021 3022 3023
  ACC        1    1    1    1    1    1    1    1    1    1    1    1    1    1    0
  Caudate    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
         
          3024 3025 3027 3028 3049 3050 3051 3052 3053 3055 3056 3057 3058 3059
  ACC        1    1    1    1    1    1    1    1    1    0    1    1    1    1
  Caudate    1    1    1    1    1    1    1    1    1    1    1    1    1    1

r$> df[df$hbcc_brain_id==2877 & df$Region=='Caudate' & 'submitted_name']  
[1] 66552 66553 66554 66556 66557 68068
```

```r
library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]

# this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# its other replicate
data = data[, ! colnames(data) %in% c('66552')]
df = df[df$submitted_name %in% colnames(data), ]

# outliers based on PCA plots
outliers = c('68080','68096', '68108', '68084', '68082')
data = data[, ! colnames(data) %in% outliers]
df = df[df$submitted_name %in% colnames(data), ]
```

```
r$> table(df$Region)                                                                    

    ACC Caudate 
     53      56 

r$> length(unique(df$hbcc_brain_id))                                                    
[1] 60

r$> table(df$Region, df$hbcc_brain_id)                                                  
         
          887 890 993 1227 1292 1331 1435 1524 1530 1683 1709 1908 2080 2371 2485 2510
  ACC       1   1   1    1    1    1    1    1    1    0    1    1    1    1    0    1
  Caudate   1   1   1    1    1    1    1    1    0    1    1    1    1    1    1    1
         
          2546 2623 2636 2646 2724 2753 2781 2842 2873 2877 3002 3003 3005 3006 3007
  ACC        1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
  Caudate    1    1    1    1    1    1    1    1    1    1    1    1    1    0    1
         
          3008 3009 3010 3011 3012 3013 3014 3016 3017 3018 3019 3020 3021 3022 3023
  ACC        1    1    1    1    1    1    1    1    1    1    1    1    1    1    0
  Caudate    1    1    1    1    1    1    1    1    1    1    1    1    1    1    1
         
          3024 3025 3027 3028 3049 3050 3051 3052 3053 3055 3056 3057 3058 3059
  ACC        1    1    1    1    1    1    1    1    1    0    1    1    1    1
  Caudate    1    1    1    1    1    1    1    1    1    1    1    1    1    1
```

And get numbers for the intersection of DGE results:

```r
res = read.csv('~/data/post_mortem/DGE_ACC_bigger_annot_04292021.csv')
resw = read.csv('~/data/post_mortem/DGE_ACC_bigger_WNH_annot_04292021.csv')
m = res[which(res$padj.FDR < .05), c('GENEID', 'hgnc_symbol')]
w = resw[which(resw$padj.FDR < .05), c('GENEID', 'hgnc_symbol')]
m2 = merge(m, w, by='GENEID', all.x=F, all.y=F)
```