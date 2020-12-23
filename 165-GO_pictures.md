# 2020-12-22 14:50:07

Let's make a few plots of specific genes similar to what Derek did in his
presentation. I'll get the genes from the enrichment results, but first we
reproduce the RNAseq analysis:

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]
# we're down from 60K to 38K samples by only looking at the ones with hgnc symbol. We might be losing too much here, so it's a step to reconsider in the future

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

library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])

data2 = cbind(data, mydata)
form = ~ Diagnosis + PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

Now let's make one plot for each gene:

```r
out_dir = '~/data/rnaseq_derek/go_pics/glutamate_raw'
df = read.delim('~/data/post_mortem/summary_results/PM_enrichment/enrichment_results_acc_geneontology_Biological_Process_noRedundant.txt')
junk = df[df$description=='glutamate receptor signaling pathway', 'userId']
go_genes = strsplit(x = junk, split=';')[[1]]
for (g in go_genes) {
    pdf(sprintf('%s/%s.pdf', out_dir, g))
    idx = which(genes$genes$hgnc_symbol == g)
    plot_df = data.frame(lcpm=lcpm[idx, ], DX=data$Diagnosis)
    boxplot(lcpm ~ DX, data=plot_df, notch=TRUE, col=(c("green","red")),
            main=g, xlab="DX", ylab='lCPM')
    dev.off()
}
```

And we can run the residualized version too:

```r
out_dir = '~/data/rnaseq_derek/go_pics/serotonin_resids'
df = read.delim('~/data/post_mortem/summary_results/PM_enrichment/enrichment_results_acc_geneontology_Biological_Process_noRedundant.txt')
junk = df[df$description=='serotonin receptor signaling pathway', 'userId']
go_genes = strsplit(x = junk, split=';')[[1]]
for (g in go_genes) {
    pdf(sprintf('%s/%s.pdf', out_dir, g))
    idx = which(genes$genes$hgnc_symbol == g)
    ens = genes$genes$ensembl_gene_id[idx]
    data2 = data.frame(t(rbind(lcpm, t(mydata))))
    fm_str = '%s ~ PC1 + PC2 + PC7 + PC8 + PC9'
    y = resid(lm(as.formula(sprintf(fm_str, ens)), data=data2))
    plot_df = data.frame(res=y, DX=data$Diagnosis)
    boxplot(res ~ DX, data=plot_df, notch=TRUE, col=(c("green","red")),
            main=g, xlab="DX", ylab='Adjusted lCPM')
    dev.off()
}
```
