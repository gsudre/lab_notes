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

