# 2020-11-24 20:13:56

Let's run the interaction of Diagnosis and Region.

```r
library(ggplot2)
library(caret)

data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)

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

num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data[, x], mydata[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex', 'Region')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(mydata[, y], data[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}
```

Now let's evaluate which PCs to keep:

```
r$> which(categ_pvals < .01, arr.ind = T)                                     
                row col
batch             1   1
Region            8   1
batch             1   2
substance_group   4   2
batch             1   3
batch             1   4
POP_CODE          6   4
batch             1   6
Sex               7   7
batch             1   8
MoD               3   8
batch             1   9
MoD               3   9
POP_CODE          6   9
batch             1  10
MoD               3  11
batch             1  12
POP_CODE          6  13
MoD               3  14
substance_group   4  14

r$> which(categ_pvals < .05, arr.ind = T)                                    
                row col
batch             1   1
Region            8   1
batch             1   2
substance_group   4   2
batch             1   3
batch             1   4
POP_CODE          6   4
MoD               3   5
substance_group   4   5
batch             1   6
Diagnosis         2   6
substance_group   4   6
batch             1   7
POP_CODE          6   7
Sex               7   7
batch             1   8
MoD               3   8
batch             1   9
MoD               3   9
POP_CODE          6   9
batch             1  10
MoD               3  11
batch             1  12
POP_CODE          6  13
MoD               3  14
substance_group   4  14
POP_CODE          6  14
Sex               7  14
```

If we want to keep the Diagnosis and Region information, we can remove 2, 3, 4,
6, 7, 8, 9, 10, 11, 12, 13, 14. But 6 is related to Diagnosis at p< .05, so I'll
keep it.

It's a shame that the first PC is conflated with batch and region, because I
know batch is a huge effect.

And just for the record, here's the same result using the numeric variables:

```
r$> which(num_pvals < .01, arr.ind = T)                      
                        row col
RINe                      4   1
clusters                  2   2
RINe                      4   2
PMI                       5   2
RINe                      4   3
C6                       11   3
clusters                  2   4
C1                        6   4
C2                        7   4
C3                        8   4
C4                        9   4
C7                       12   4
pcnt_optical_duplicates   1   8
clusters                  2   8
PMI                       5   8
C7                       12   8
C1                        6   9
C2                        7   9
C3                        8   9
C4                        9   9
C7                       12   9
Age                       3  11
clusters                  2  12
RINe                      4  13
C1                        6  13
C2                        7  13
C3                        8  13
C4                        9  13
C7                       12  13
Age                       3  14
C2                        7  14
C3                        8  14
C4                        9  14
```

It doesn't add anything to the previous PCs. It would be hard to actually.

```r
data2 = cbind(data, mydata)
form = ~ Diagnosis*Region + PC2 + PC3 + PC4 + PC7 + PC8 + PC9 + PC10 + PC11 + PC12 + PC13 + PC14
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
res = topTable(fit2, coef='DiagnosisCase:RegionCaudate', number=Inf, sort.by='P')
write.csv(res, file='~/data/rnaseq_derek/results_interaction.csv', row.names=F)
```

Now let's check if there's anything interesting about the interaction term:

```r
library(WebGestaltR)
library(gage)
compute_gage = function(ranks, out_root, useFDR=T) {
    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'pathway_KEGG')
    for (db in DBs) {
        cat(out_root, db, '\n')
        suffix = ifelse(useFDR, '_q1', '_p06')
        myc = ifelse(useFDR, .1, .06)
        myqp = ifelse(useFDR, 'q.val', 'p.val')

        gs = loadGeneSet(enrichDatabase=db)
        gmt = gs$geneSet
        a = idMapping(inputGene=gmt$gene, sourceIdType='entrezgene',
                    targetIdType='genesymbol')
        gmt2 = merge(gmt, a$mapped[, c('userId', 'geneSymbol')], by.x = 'gene',
                    by.y='userId', all.x=F, all.y=F)
        # and convert it to lists
        mylist = list()
        for (s in unique(gmt2$geneSet)) {
            mylist[[s]] = unique(gmt2$geneSymbol[gmt2$geneSet==s])
        }
        resSameDir = gage(ranks, gsets = mylist, compare='unpaired',
                          set.size=c(5, 800), same.dir=T)
        sigSameDir = sigGeneSet(resSameDir, cutoff=myc, qpval=myqp)
        resOneDir = as.data.frame(rbind(sigSameDir$less, sigSameDir$greater))
        resOneDir$Enrichment = ifelse(resOneDir$stat.mean > 0,
                                      "Up-regulated", "Down-regulated")
        resOneDir = merge(resOneDir, gs$geneSetDes, by.x=0, by.y='geneSet',
                          sort=F)
        out_fname = sprintf('~/data/post_mortem/gage_%s_OneDir_%s%s.csv',
                            out_root, db, suffix)
        write.csv(resOneDir, file=out_fname, row.names=F)
        resBothDir = gage(ranks, gsets = mylist, compare='unpaired',
                          set.size=c(5, 800), same.dir=F)
        sigBothDir = sigGeneSet(resBothDir, cutoff=myc, qpval=myqp)
        sigBothDir = merge(sigBothDir$greater, gs$geneSetDes, by.x=0,
                           by.y='geneSet', sort=F)
        out_fname = sprintf('~/data/post_mortem/gage_%s_BothDir_%s%s.csv',
                            out_root, db, suffix)
        write.csv(sigBothDir, file=out_fname, row.names=F)
    }
    # my own lists
    DBs = c('adhd_genes')
    for (db in DBs) {
        cat(out_root, db, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        gmt = readGmt(db_file) # already in gene symbols
        # and convert it to lists
        mylist = list()
        for (s in unique(gmt$geneSet)) {
            mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
        }
        gmt_desc = gmt[, c('geneSet', 'description')]
        gmt_desc = gmt_desc[!duplicated(gmt_desc),]

        resSameDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=T)
        sigSameDir = sigGeneSet(resSameDir, cutoff=.06, qpval='p.val')
        resOneDir = as.data.frame(rbind(sigSameDir$less, sigSameDir$greater))
        resOneDir$Enrichment = ifelse(resOneDir$stat.mean > 0, "Up-regulated", "Down-regulated")
        resOneDir = merge(resOneDir, gmt_desc, by.x=0, by.y='geneSet', sort=F)
        out_fname = sprintf('~/data/post_mortem/gage_%s_OneDir_%s_p06.csv',
                            out_root, db)
        write.csv(resOneDir, file=out_fname, row.names=F)
        resBothDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=F)
        sigBothDir = sigGeneSet(resBothDir, cutoff=.06, qpval='p.val')
        sigBothDir = merge(sigBothDir$greater, gmt_desc, by.x=0, by.y='geneSet', sort=F)
        out_fname = sprintf('~/data/post_mortem/gage_%s_BothDir_%s_p06.csv',
                            out_root, db)
        write.csv(sigBothDir, file=out_fname, row.names=F)
    }
}

tmp = res
dup_genes = tmp$hgnc_symbol[duplicated(tmp$hgnc_symbol)]
res = tmp[!tmp$hgnc_symbol %in% dup_genes, ]
for (g in dup_genes) {
  gene_data = tmp[tmp$hgnc_symbol==g, ]
  best_res = which.min(gene_data$P.Value)
  res = rbind(res, gene_data[best_res, ])
}
ranks = -log(res$P.Value) * sign(res$logFC)
names(ranks) = res$hgnc_symbol
ranks = sort(ranks, decreasing=T)
compute_gage(ranks, 'rnaseq_DXbyRegion', useFDR=F)
```

We get a GABAergic synapse KEGG pathway, but hard to focus on it among the many
results. It's nominal only. We also get a nominal hit for neurotransmitter
receptor activity, but that's it.