# 2021-01-14 17:08:56

I think I figured out how to add the subtype annotations to the transcriptome.
They were in the GTF file I was using all along, but the R function wasn't
importing them. So, in Python:

```bash
conda create --name rnaseq
conda activate rnaseq
conda install ipython
conda install pandas
```

```python
from gtfparse import read_gtf
fn = '/Users/sudregp/data/post_mortem/Homo_sapiens.GRCh38.97.gtf'
df = read_gtf(fn)
out_fn = '/Users/sudregp/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv'
my_cols = ['transcript_id', 'transcript_biotype', 'gene_id', 'gene_biotype']
df[my_cols].to_csv(out_fn)
```

Now we can just add them to our usual analysis and see where it goes. Starting
with DGE:

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
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

At this point, no filtering has been done, except for keeping only autosomes.
And I added all annotations I wanted. Now, it's just a matter of filtering in
any way we want.

```r
subtype = 'pseudogene'

cat('Starting with', nrow(tx_meta), 'variables\n')
keep_me = grepl(tx_meta$gene_biotype, pattern=sprintf('%s$', subtype))
cat('Keeping', sum(keep_me), subtype, 'variables\n')
my_count_matrix = count_matrix[keep_me, ]
my_tx_meta = tx_meta[keep_me, ]

# removing variables with zero or near-zero variance
library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(my_count_matrix), method = pp_order)
X = t(predict(pp, t(my_count_matrix)))
cat('Keeping', nrow(X), 'after NZ and NZV filtering\n')

# checking which PCs are associated with our potential nuiscance variables
set.seed(42)
mypca <- prcomp(t(X), scale=TRUE)
# how many PCs to keep... using Kaiser thredhold, close to eigenvalues < 1
library(nFactors)
eigs <- mypca$sdev^2
nS = nScree(x=eigs)
keep_me = seq(1, nS$Components$nkaiser)

mydata = data.frame(mypca$x[, keep_me])
# create main metadata data frame including metadata and PCs
data.pm = cbind(data, mydata)
rownames(data.pm) = data$hbcc_brain_id
cat('Using', nS$Components$nkaiser, 'PCs from possible', ncol(X), '\n')

# check which PCs are associated at nominal p<.01
num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data.pm[, x], data.pm[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level')
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
use_pcs = unique(c(which(num_pvals < .01, arr.ind = T)[, 'col'],
                   which(categ_pvals < .01, arr.ind = T)[, 'col']))
fm_str = sprintf('~ Diagnosis + %s', paste0(pc_vars[use_pcs], collapse = ' + '))
cat('Found', length(use_pcs), 'PCs p < .01\n')
cat('Using formula:', fm_str, '\n')

# removing variables with low expression
library(edgeR)
design=model.matrix(as.formula(fm_str), data=data.pm)
isexpr <- filterByExpr(X, design=design)
countsExpr = X[isexpr,]
metaExpr = data.frame(GENEID = substr(rownames(countsExpr), 1, 15))
metaExpr = merge(metaExpr, my_tx_meta, by='GENEID', sort=F)
cat('Keeping', nrow(countsExpr), 'after expression filtering\n')

# preparing DESeqData and running main analysis
countdata = round(countsExpr)
colnames(countdata) = rownames(data.pm)
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = data.pm,
                              design = as.formula(fm_str))
dds <- DESeq(dds)
res <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05)
cat('FDR q < .05\n')
print(summary(res))
gene_ids = rownames(res)[res$padj < .05]
print(gene_ids)

# plotting each of the significant genes
library(ggpubr)
library(ggbeeswarm)
quartz()
myplots = list()
clrs = c("green3", "red")
for (g in 1:length(gene_ids)) {
    cat(gene_ids[g], '\n')
    d <- plotCounts(dds, gene=gene_ids[g], intgroup="Diagnosis",
                    returnData=TRUE)
    p = (ggplot(d, aes(x=Diagnosis, y=count, color = Diagnosis,
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
         theme_bw() +
         ggtitle(gene_ids[g]))
    myplots[[g]] = p
}
p = ggarrange(plotlist=myplots)
annotate_figure(p, sprintf('DGE %s %s FDR q<.05', subtype, myregion))

library(IHW)
resIHW <- results(dds, name = "Diagnosis_Case_vs_Control", alpha = 0.05,
                  filterFun=ihw)
cat('IHW q < .05\n')
print(summary(resIHW))
gene_ids = rownames(resIHW)[resIHW$padj < .05]
print(gene_ids)
quartz()
myplots = list()
clrs = c("green3", "red")
for (g in 1:length(gene_ids)) {
    cat(gene_ids[g], '\n')
    d <- plotCounts(dds, gene=gene_ids[g], intgroup="Diagnosis",
                    returnData=TRUE)
    p = (ggplot(d, aes(x=Diagnosis, y=count, color = Diagnosis,
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
         theme_bw() +
         ggtitle(gene_ids[g]))
    myplots[[g]] = p
}
p = ggarrange(plotlist=myplots)
annotate_figure(p, sprintf('DGE %s %s IHW q<.05', subtype, myregion))
```

