# 2020-11-25 08:57:32

Let me see how the WG results change for RNAseq if using -logP*sign(logFC) for
rank.

```r
# bw
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
load(sprintf('%s/rnaseq_results_11122020.rData', data_dir))
ncpu=2

region='acc'
# region='caudate'
eval(parse(text=sprintf('res = rnaseq_%s', region)))

ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

# my own GMTs
for (db in c('adhd_genes', sprintf('%s_developmental', region))) {
    cat(region, db, '\n')
    project_name = sprintf('%s_%s', region, db)
    db_file = sprintf('~/data/post_mortem/%s.gmt', db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabaseFile=db_file,
                                enrichDatabaseType="genesymbol",
                                interestGene=tmp2,
                                outputDirectory = data_dir,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                minNum=3, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800)
    out_fname = sprintf('%s/WG3_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    project_name = sprintf('%s_%s', region, db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800)
    out_fname = sprintf('%s/WG3_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```
