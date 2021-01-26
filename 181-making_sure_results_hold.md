# 2021-01-26 10:25:17

Let's make sure our WG and PRs results still hold when using the results from
note 180.

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=2

load('~/data/post_mortem/DGE_01262021.RData')

region='acc'

ranks = -log(dge_acc_pc$pvalue) * sign(dge_acc_pc$log2FoldChange)
tmp2 = data.frame(geneid=substring(rownames(dge_acc_pc), 1, 15), rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

DBs = c(sprintf('my_%s_sets', region), # just to get GWAS and TWAS sets
        sprintf('%s_manySets_co0.990', region),
        sprintf('%s_manySets_co0.950', region),
        sprintf('%s_manySets', region))
for (db in DBs) {
    cat(region, db, '\n')
    db_file = sprintf('~/data/post_mortem/%s.gmt', db)
    project_name = sprintf('WG6_%s_%s_10K', region, db)
    enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                        organism="hsapiens",
                        enrichDatabaseFile=db_file,
                        enrichDatabaseType="genesymbol",
                        interestGene=tmp2,
                        outputDirectory = data_dir,
                        interestGeneType="ensembl_gene_id",
                        sigMethod="top", topThr=50,
                        minNum=3, projectName=project_name,
                        isOutput=T, isParallel=T,
                        nThreads=ncpu, perNum=10000, maxNum=800))
    out_fname = sprintf('%s/WG6_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}

DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    project_name = sprintf('WG6_%s_%s_10K', region, db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=50,
                                outputDirectory = data_dir,
                                minNum=5, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000)
    out_fname = sprintf('%s/WG6_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

They do! I'll script it out and run it for Caudate as well:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=2

load('~/data/post_mortem/DGE_01262021.RData')

subtypes = c('pc', 'pg', 'lnc')

for (region in c('caudate', 'acc')) {
    for (st in subtypes) {
        res_str = ifelse(region == 'acc', sprintf('dge_acc_%s', st),
                         sprintf('dge_cau_%s', st))
        ranks_str = sprintf('ranks = -log(%s$pvalue) * sign(%s$log2FoldChange)',
                            res_str, res_str)
        gid_str = sprintf('geneid=substring(rownames(%s), 1, 15)', res_str)
        
        eval(parse(text=ranks_str))
        eval(parse(text=gid_str))

        tmp2 = data.frame(geneid=geneid, rank=ranks)
        tmp2 = tmp2[order(ranks, decreasing=T),]

        DBs = c(sprintf('my_%s_sets', region), # just to get GWAS and TWAS sets
                sprintf('%s_manySets_co0.990', region),
                sprintf('%s_manySets_co0.950', region),
                sprintf('%s_manySets', region))
        for (db in DBs) {
            cat(res_str, db, '\n')
            db_file = sprintf('~/data/post_mortem/%s.gmt', db)
            project_name = sprintf('WG6_%s_%s_10K', res_str, db)
            enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabaseFile=db_file,
                                enrichDatabaseType="genesymbol",
                                interestGene=tmp2,
                                outputDirectory = data_dir,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=50,
                                minNum=3, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800))
            if (class(enrichResult) != "try-error") {
                out_fname = sprintf('%s/WG6_%s_%s_10K.csv', data_dir, res_str, db)
                write.csv(enrichResult, file=out_fname, row.names=F)
            }
        }

        DBs = c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant')
        for (db in DBs) {
            cat(res_str, db, '\n')
            project_name = sprintf('WG6_%s_%s_10K', res_str, db)

            enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase=db,
                                        interestGene=tmp2,
                                        interestGeneType="ensembl_gene_id",
                                        sigMethod="top", topThr=50,
                                        outputDirectory = data_dir,
                                        minNum=5, projectName=project_name,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, perNum=10000))
            if (class(enrichResult) != "try-error") {
                out_fname = sprintf('%s/WG6_%s_%s_10K.csv', data_dir, res_str, db)
                write.csv(enrichResult, file=out_fname, row.names=F)
            }
        }
    }
}
```

So, it turns out that both our sets and the GO sets only work for the
protein_coding genes! There are no intersections for the pseudogenes or lncRNA!

