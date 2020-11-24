# 2020-11-24 14:44:39

Let's look at the results David sent through gage first.

```r
compute_gage = function(ranks, out_root, isACC=T, useFDR=T) {
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
for (r in c('acc', 'caudate')) {
    df = read.csv(sprintf('~/data/post_mortem/david_pca_%s.csv', r))
    for (p in colnames(df)[2:4]) {
        tmp = df[!is.na(df[, p]), c('geneName', p)]
        colnames(tmp) = c('hgnc_symbol', 'P.Value')
        dup_genes = tmp$hgnc_symbol[duplicated(tmp$hgnc_symbol)]
        res = tmp[!tmp$hgnc_symbol %in% dup_genes, ]
        for (g in dup_genes) {
            gene_data = tmp[tmp$hgnc_symbol==g, ]
            best_res = which.min(gene_data$P.Value)
            res = rbind(res, gene_data[best_res, ])
        }
        ranks = -log(res$P.Value)
        names(ranks) = res$hgnc_symbol
        ranks = sort(ranks, decreasing=T)
        out_name = sprintf('iso_%s_%s', r, p)
        compute_gage(ranks, out_name, isACC=(r=='acc'), useFDR=F)
    }
}
```