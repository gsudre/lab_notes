# 2021-02-01 17:19:13

Let's make some more plots for the presentation to the HBCC bunch.

## ACC developmental sets

```r
df = read.csv('~/data/post_mortem/WG7_dge_acc[[\"protein_coding\"]]_my_acc_sets_10K.csv')
df = df[c(1, 2, 3, 5, 6, 8), c('link', 'normalizedEnrichmentScore', 'pValue')]
df[1, 'link'] = 'pan-developmental'
df$link = as.character(df$link)
df$link = factor(df$link, levels=c("pan-developmental", "prenatal",
                                   "infant (0-2 yrs)", "child (3-11 yrs)",
                                   "adolescent (12-19 yrs)", "adult (>19 yrs)"))
df$Direction = ifelse(df$normalizedEnrichmentScore > 0, 'up', 'down')
df$str = sapply(df$pValue, function(x) sprintf('p = %.1e', x))
df$star_pos = abs(df$normalizedEnrichmentScore) + .1
df[df$pValue ==0, 'str'] = 'p < 1e-5'
stars.df <- df[df$link %in% c("pan-developmental", 'adult (>19 yrs)',
                              "prenatal"), c('link', 'star_pos')]
library(ggplot2)
ggplot(data=df, aes(x=link, y=abs(normalizedEnrichmentScore), fill=Direction)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=str), vjust=1.6, size=5, color='white')+
  theme_minimal() + ylab('Absolute Normalized Enrichment Score') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=14)) +
  geom_text(data = stars.df, aes(y = star_pos, fill=NA), label = "***")
```

## ACC biological processes

```r
df = read.csv('~/data/post_mortem/WG7_dge_acc[[\"protein_coding\"]]_geneontology_Biological_Process_noRedundant_10K.csv')
df = df[df$FDR < .05, c('description', 'normalizedEnrichmentScore')]
df$Direction = ifelse(df$normalizedEnrichmentScore > 0, 'up', 'down')
df$star_pos = abs(df$normalizedEnrichmentScore)
my_order = order(df$star_pos)
df$description = as.character(df$description)
df$description = factor(df$description, levels=df$description[my_order])

library(ggplot2)
ggplot(data=df, aes(x=description, y=abs(normalizedEnrichmentScore), fill=Direction)) +
  geom_bar(stat="identity")+ coord_flip() +
  theme_minimal() + ylab('Absolute Normalized Enrichment Score') + xlab('') +
  theme(axis.text.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=14))
```

## Leading genes for serotonin receptor signaling pathway

```r
df = read.csv('~/data/post_mortem/WG7_dge_acc[[\"protein_coding\"]]_geneontology_Biological_Process_noRedundant_10K.csv')
genes = df[df$description=='serotonin receptor signaling pathway', 'userId']
gene_list = strsplit(genes, ';')[[1]]
subtype = 'protein_coding'

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

cat('Starting with', nrow(tx_meta), 'variables\n')
keep_me = grepl(tx_meta$gene_biotype, pattern=sprintf('%s$', subtype))
cat('Keeping', sum(keep_me), subtype, 'variables\n')
my_count_matrix = count_matrix[keep_me, ]
my_tx_meta = tx_meta[keep_me, ]

# removing variables where more than half of the subjects have zero counts
keep_me = rowSums(my_count_matrix==0) < .25*ncol(my_count_matrix)
my_count_matrix = my_count_matrix[keep_me, ]
cat('Keeping', nrow(my_count_matrix), 'after zero removal\n')

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
        res = cor.test(data.pm[, x], data.pm[, y], method='spearman')
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
# only use the ones not related to Diagnosis
keep_me = c()
for (pc in use_pcs) {
    keep_me = c(keep_me, categ_pvals['Diagnosis', pc] > .05)
}
use_pcs = use_pcs[keep_me]

fm_str = sprintf('~ Diagnosis + %s', paste0(pc_vars[use_pcs],
                                            collapse = ' + '))
cat('Found', length(use_pcs), 'PCs p < .01\n')
cat('Using formula:', fm_str, '\n')

# scaling PCs to assure convergence
for (var in pc_vars[use_pcs]) {
    data.pm[, var] = scale(data.pm[, var])
}

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
# because DESeq doesn't remove outliers if there are continuous variables
# in the formula, we need to do this iteratively
nOutliers = Inf
myCounts = round(countsExpr)
while (nOutliers > 0) {
    dds <- DESeqDataSetFromMatrix(countData = myCounts,
                                colData = data.pm,
                                design = as.formula(fm_str))
    cat('Processing', nrow(dds), 'variables.\n')
    dds <- DESeq(dds)
    maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
    # outlier cut-off uses the 99% quantile of the F(p,m-p) distribution (with 
    # p the number of parameters including the intercept and m number of
    # samples).
    m <- ncol(dds)
    # number or parameters (PCs + Diagnosis + intercept)
    p <- length(use_pcs) + 2
    co = qf(.99, p, m - p)
    keep_me = which(maxCooks < co)
    nOutliers = nrow(myCounts) - length(keep_me)
    cat('Found', nOutliers, 'outliers.\n')
    myCounts = round(myCounts)[keep_me, ]
}

gnames = data.frame(full=rownames(counts(dds)),
                    nov=substr(rownames(counts(dds)), 1, 15))
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
gnames = merge(gnames, mart, by.x='nov', by.y='ensembl_gene_id')
keep_me = gnames$nov %in% gene_list
gene_ids = gnames[keep_me, ]

# plotting each of the significant genes
library(ggpubr)
library(ggbeeswarm)
quartz()
myplots = list()
clrs = c("green3", "red")
for (g in 1:nrow(gene_ids)) {
    cat(gene_ids[g, 'nov'], '\n')
    d <- plotCounts(dds, gene=gene_ids[g, 'full'], intgroup="Diagnosis",
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
        theme_bw() + #theme(legend.position = "none") + 
        ggtitle(gene_ids[g, 'hgnc_symbol']))
    myplots[[g]] = p
}
p = ggarrange(plotlist=myplots)
print(p)
```

## Caudate developmental sets

```r
df = read.csv('~/data/post_mortem/WG7_dge_cau[[\"protein_coding\"]]_my_caudate_sets_10K.csv')
df = df[c(1, 2, 4, 5, 7, 8), c('link', 'normalizedEnrichmentScore', 'pValue')]
df[3, 'link'] = 'pan-developmental'
df$link = as.character(df$link)
df$link = factor(df$link, levels=c("pan-developmental", "prenatal",
                                   "infant (0-2 yrs)", "child (3-11 yrs)",
                                   "adolescent (12-19 yrs)", "adult (>19 yrs)"))
df$Direction = ifelse(df$normalizedEnrichmentScore > 0, 'up', 'down')
df$str = sapply(df$pValue, function(x) sprintf('p = %.2f', x))
df[df$pValue ==0, 'str'] = 'p < 1e-5'
library(ggplot2)
ggplot(data=df, aes(x=link, y=abs(normalizedEnrichmentScore), fill=Direction)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=str), vjust=1.6, size=5, color='white')+
  theme_minimal() + ylab('Absolute Normalized Enrichment Score') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=14))
```

## Caudate biological processes

```r
df = read.csv('~/data/post_mortem/WG7_dge_cau[[\"protein_coding\"]]_geneontology_Biological_Process_noRedundant_10K.csv')
df$score = abs(df$normalizedEnrichmentScore)
df = df[order(df$score, decreasing=T)[1:10], c('description', 'score', 'normalizedEnrichmentScore')]
df$Direction = ifelse(df$normalizedEnrichmentScore > 0, 'up', 'down')
df$star_pos = df$score + .1
my_order = order(df$star_pos)
df$description = as.character(df$description)
df$description = factor(df$description, levels=df$description[my_order])
stars.df <- df[df$description %in% c("sperm motility",
                                     "microtubule bundle formation"),]
library(ggplot2)
ggplot(data=df, aes(x=description, y=abs(normalizedEnrichmentScore), fill=Direction)) +
  geom_bar(stat="identity")+ coord_flip() +
  theme_minimal() + ylab('Absolute Normalized Enrichment Score') + xlab('') +
  theme(axis.text.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=14)) +
  geom_text(data = stars.df, aes(y = star_pos, fill=NA), label = "***")
```

## Caudate and ACC overlaps

```r
library(GeneOverlap)
library(VennDiagram)
load('~/data/post_mortem/DGE_01272021.RData')

all_res = c()
subtypes = list(pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('pc', 'lnc', 'pg')) {
    res.acc = dge_acc[[subtypes[[st]]]]
    res.cau = dge_cau[[subtypes[[st]]]]
    
    both_res = merge(as.data.frame(res.acc), as.data.frame(res.cau), by=0,
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    t = .005
    prs_genes = both_res[both_res$pvalue.prs < t, 'Row.names']
    dx_genes = both_res[both_res$pvalue.dx < t, 'Row.names']
    go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                genome.size=nrow(both_res))
    go.obj <- testGeneOverlap(go.obj)
    inter = intersect(prs_genes, dx_genes)
    pval1 = getPval(go.obj)
    pval2 = NA
    this_res = c(subtypes[[st]], t, 'abs', length(prs_genes),
                    length(dx_genes), length(inter), pval1, pval2)
    all_res = rbind(all_res, this_res)

    out_fname = sprintf('~/tmp/%s.png', st)
    venn.plot = venn.diagram(list(ACC = dx_genes, Caudate = prs_genes),
                             euler.d=TRUE, fill=c('red','green'),
                             main=subtypes[[st]],
                             filename=out_fname, imagetype='png',
                             main.cex=5, main.fontfamily='sans', cex=4,
                             fontfamily='sans', cat.fontfamily='sans',
                             cat.cex=4, category.names=c('', ''))
 }
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
```

## ACC PRS DX overlap

```r
library(GeneOverlap)
load('~/data/post_mortem/DGE_PRS_01272021.RData')
load('~/data/post_mortem/DGE_01272021.RData')

p = "PRS0.005000"
all_res = c()
st = 'pc'
subtypes = list(pc='protein_coding', lnc='lncRNA', pg='pseudogene')
res.dx = dge_acc[[subtypes[[st]]]]
res_str = sprintf('res.prs = dgePRS_acc_%s[["%s"]]', st, p)
eval(parse(text=res_str))

both_res = merge(as.data.frame(res.dx), as.data.frame(res.prs), by=0,
                    all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
t = .01
prs_genes = both_res[both_res$pvalue.prs < t, 'Row.names']
dx_genes = both_res[both_res$pvalue.dx < t, 'Row.names']
out_fname = sprintf('~/tmp/acc_%s.png', st)
venn.plot = venn.diagram(list(DX = dx_genes, PRS = prs_genes),
                            euler.d=TRUE, fill=c('red','green'),
                            main=subtypes[[st]],
                            filename=out_fname, imagetype='png',
                            main.cex=5, main.fontfamily='sans', cex=3,
                            fontfamily='sans', cat.fontfamily='sans',
                            cat.cex=4, category.names=c('', ''))
go.obj <- newGeneOverlap(prs_genes, dx_genes,
                            genome.size=nrow(both_res))
go.obj <- testGeneOverlap(go.obj)
inter = intersect(prs_genes, dx_genes)
pval1 = getPval(go.obj)
pval2 = NA
this_res = c(subtypes[[st]], p, t, 'abs', length(prs_genes),
                length(dx_genes), length(inter), pval1, pval2)
all_res = rbind(all_res, this_res)
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
```

## Caudate PRS DX overlap

```r
library(GeneOverlap)
load('~/data/post_mortem/DGE_PRS_01272021.RData')
load('~/data/post_mortem/DGE_01272021.RData')

p = "PRS0.005000"
all_res = c()
st = 'pc'
subtypes = list(pc='protein_coding', lnc='lncRNA', pg='pseudogene')
res.dx = dge_cau[[subtypes[[st]]]]
res_str = sprintf('res.prs = dgePRS_cau_%s[["%s"]]', st, p)
eval(parse(text=res_str))

both_res = merge(as.data.frame(res.dx), as.data.frame(res.prs), by=0,
                    all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
t = .01
prs_genes = both_res[both_res$pvalue.prs < t, 'Row.names']
dx_genes = both_res[both_res$pvalue.dx < t, 'Row.names']
out_fname = sprintf('~/tmp/cau_%s.png', st)
venn.plot = venn.diagram(list(DX = dx_genes, PRS = prs_genes),
                            euler.d=TRUE, fill=c('red','green'),
                            main=subtypes[[st]],
                            filename=out_fname, imagetype='png',
                            main.cex=5, main.fontfamily='sans', cex=3,
                            fontfamily='sans', cat.fontfamily='sans',
                            cat.cex=4, category.names=c('', ''))
go.obj <- newGeneOverlap(prs_genes, dx_genes,
                            genome.size=nrow(both_res))
go.obj <- testGeneOverlap(go.obj)
inter = intersect(prs_genes, dx_genes)
pval1 = getPval(go.obj)
pval2 = NA
this_res = c(subtypes[[st]], p, t, 'abs', length(prs_genes),
                length(dx_genes), length(inter), pval1, pval2)
all_res = rbind(all_res, this_res)
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
```
