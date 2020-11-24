# 2020-11-23 20:44:58

LEt's run all our results through GAGE, as it seems to preserve most of our
results and its documentation is quite up-to-date.

```r
library(WebGestaltR)
library(gage)

load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
tmp = rnaseq_acc

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

db = 'geneontology_Biological_Process_noRedundant'
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

resSameDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=T)
sigSameDir = sigGeneSet(resSameDir, cutoff=.05, qpval='q.val')
resOneDir = as.data.frame(rbind(sigSameDir$less, sigSameDir$greater))
resOneDir$Enrichment = ifelse(resOneDir$stat.mean > 0, "Up-regulated", "Down-regulated")
resOneDir = merge(resOneDir, gs$geneSetDes, by.x=0, by.y='geneSet', sort=F)
out_fname = sprintf('~/data/rnaseq_derek/gage_rnaseq_acc_OneDir_%s_q05.csv', db)
write.csv(resOneDir, file=out_fname, row.names=F)
resBothDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=F)
sigBothDir = sigGeneSet(resBothDir, cutoff=.05, qpval='q.val')
sigBothDir = merge(sigBothDir$greater, gs$geneSetDes, by.x=0, by.y='geneSet', sort=F)
out_fname = sprintf('~/data/rnaseq_derek/gage_rnaseq_acc_BothDir_%s_q05.csv', db)
write.csv(sigBothDir, file=out_fname, row.names=F)
```

Then, if using one of our sets:

```r
db_file = '~/data/post_mortem/acc_developmental.gmt'
gmt = readGmt(db_file) # already in gene symbols
# and convert it to lists
mylist = list()
for (s in unique(gmt$geneSet)) {
    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
}
gmt_desc = gmt[, c('geneSet', 'description')]
gmt_desc = gmt_desc[!duplicated(gmt_desc),]

resSameDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=T)
sigSameDir = sigGeneSet(resSameDir, cutoff=.05, qpval='p.val')
resOneDir = as.data.frame(rbind(sigSameDir$less, sigSameDir$greater))
resOneDir$Enrichment = ifelse(resOneDir$stat.mean > 0, "Up-regulated", "Down-regulated")
resOneDir = merge(resOneDir, gmt_desc, by.x=0, by.y='geneSet', sort=F)
out_fname = '~/data/rnaseq_derek/gage_rnaseq_acc_OneDir_dev_p05.csv'
write.csv(resOneDir, file=out_fname, row.names=F)
resBothDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=F)
sigBothDir = sigGeneSet(resBothDir, cutoff=.05, qpval='p.val')
sigBothDir = merge(sigBothDir$greater, gmt_desc, by.x=0, by.y='geneSet', sort=F)
out_fname = '~/data/rnaseq_derek/gage_rnaseq_acc_BothDir_dev_p05.csv'
write.csv(sigBothDir, file=out_fname, row.names=F)
```

And finally, with the ADHD gene sets:

```r
db_file = '~/data/post_mortem/adhd_genes.gmt'
gmt = readGmt(db_file) # already in gene symbols
# and convert it to lists
mylist = list()
for (s in c('GWAS1', 'GWAS', 'TWAS1', 'TWAS2', 'TWAS', 'CNV1', 'CNV2')) {
    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
}
gmt_desc = gmt[, c('geneSet', 'description')]
gmt_desc = gmt_desc[!duplicated(gmt_desc),]

resSameDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=T)
sigSameDir = sigGeneSet(resSameDir, cutoff=.06, qpval='p.val')
resOneDir = as.data.frame(rbind(sigSameDir$less, sigSameDir$greater))
resOneDir$Enrichment = ifelse(resOneDir$stat.mean > 0, "Up-regulated", "Down-regulated")
resOneDir = merge(resOneDir, gmt_desc, by.x=0, by.y='geneSet', sort=F)
out_fname = '~/data/rnaseq_derek/gage_rnaseq_acc_OneDir_adhdGenes_p06.csv'
write.csv(resOneDir, file=out_fname, row.names=F)
resBothDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=F)
sigBothDir = sigGeneSet(resBothDir, cutoff=.06, qpval='p.val')
sigBothDir = merge(sigBothDir$greater, gmt_desc, by.x=0, by.y='geneSet', sort=F)
out_fname = '~/data/rnaseq_derek/gage_rnaseq_acc_BothDir_adhdGenes_p06.csv'
write.csv(sigBothDir, file=out_fname, row.names=F)
```

# 2020-11-24 06:23:44

Let's script this out:

```r
compute_gage = function(ranks, out_root, isACC=T, useFDR=T) {
    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'pathway_KEGG')
    for (db in DBs) {
        cat(out_root, db, '\n')
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
        resSameDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=T)
        sigSameDir = ifelse(useFDR,
                            sigGeneSet(resSameDir, cutoff=.05, qpval='q.val'),
                            sigGeneSet(resSameDir, cutoff=.06, qpval='p.val'))
        suffix = ifelse(useFDR, '_q05', 'p06')
        resOneDir = as.data.frame(rbind(sigSameDir$less, sigSameDir$greater))
        resOneDir$Enrichment = ifelse(resOneDir$stat.mean > 0, "Up-regulated", "Down-regulated")
        resOneDir = merge(resOneDir, gs$geneSetDes, by.x=0, by.y='geneSet', sort=F)
        out_fname = sprintf('~/data/post_mortem/gage_%s_OneDir_%s%s.csv',
                            out_root, db, suffix)
        write.csv(resOneDir, file=out_fname, row.names=F)
        resBothDir = gage(ranks, gsets = mylist, compare='unpaired', set.size=c(5, 800), same.dir=F)
        sigBothDir = ifelse(useFDR,
                            sigGeneSet(resBothDir, cutoff=.05, qpval='q.val'),
                            sigGeneSet(resBothDir, cutoff=.06, qpval='p.val'))
        sigBothDir = merge(sigBothDir$greater, gs$geneSetDes, by.x=0, by.y='geneSet', sort=F)
        out_fname = sprintf('~/data/post_mortem/gage_%s_BothDir_%s_%s.csv',
                            out_root, db, suffix)
        write.csv(sigBothDir, file=out_fname, row.names=F)
    }
    # my own lists
    DBs = ifelse(isACC, 'acc_developmental', 'caudate_developmental')
    DBs = c(DBs, 'adhd_genes')
    for (db in DBs) {
        cat(out_root, db, isACC, '\n')
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
```

```r
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

tmp = rnaseq_acc
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
compute_gage(ranks, 'rnaseq_acc', T)

tmp = rnaseq_caudate
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
compute_gage(ranks, 'rnaseq_caudate', F)
```

And now we run the same stuff for the imputations:

```r
data_dir = '~/data/expression_impute/'
phenotypes = read.table(sprintf('%s/phenos.txt', data_dir))[,1]
G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]

for (phen in phenotypes) {
    isACC = !(grepl(x=phen, pattern='Caudate') || grepl(x=phen, pattern='ATR'))
    for (md in c('EN', 'MASHR')) {
        res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen),
                        header=1)
        id_num = sapply(res$gene, function(x) strsplit(x=x, split='\\.')[[1]][1])
        dups = duplicated(id_num)
        id_num = id_num[!dups]
        res$id_num = id_num
        imnamed = res$id_num %in% G_list$ensembl_gene_id
        res = res[imnamed, ]
        G_list2 = merge(G_list, res, by.x='ensembl_gene_id', by.y='id_num')
        imautosome = which(G_list2$chromosome_name != 'X' &
                        G_list2$chromosome_name != 'Y' &
                        G_list2$chromosome_name != 'MT')
        G_list2 = G_list2[imautosome, ]

        # signs for zscore match the ones for effect
        tmp = G_list2[, c('hgnc_symbol', 'pvalue', 'zscore')]
        dup_genes = tmp$hgnc_symbol[duplicated(tmp$hgnc_symbol)]
        res = tmp[!tmp$hgnc_symbol %in% dup_genes, ]
        for (g in dup_genes) {
            gene_data = tmp[tmp$hgnc_symbol==g, ]
            best_res = which.min(gene_data$pvalue)
            res = rbind(res, gene_data[best_res, ])
        }
        ranks = -log(res$pvalue) * sign(res$zscore)
        names(ranks) = res$hgnc_symbol
        ranks = sort(ranks, decreasing=T)
        out_root = sprintf('imp_%s_%s', md, phen)
        compute_gage(ranks, out_root, isACC, useFDR=F)
    }
}
```