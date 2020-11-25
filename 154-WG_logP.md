# 2020-11-25 08:57:32

Let me see how the WG results change if using -logP*sign(logFC) for
rank.

```r
# bw
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
load(sprintf('%s/rnaseq_results_11122020.rData', data_dir))
ncpu=8

# region='acc'
region='caudate'
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

And do the same thing for the imputation results:

```r
# bw
library(WebGestaltR)

data_dir = '~/data/expression_impute/'
phenotypes = read.table(sprintf('%s/phenos.txt', data_dir))[,1]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
ncpu=8

for (phen in phenotypes) {
    if (grepl(x=phen, pattern='Caudate') || grepl(x=phen, pattern='ATR')) {
        region = 'caudate'
    } else {
        region = 'ACC'
    }
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
        tmp = G_list2[, c('hgnc_symbol', 'pvalue', 'zscore')]
        ranks = -log(tmp$pvalue) * sign(tmp$zscore)
        tmp2 = data.frame(hgnc_symbol=tmp$hgnc_symbol, rank=ranks)
        tmp2 = tmp2[order(ranks, decreasing=T),]
        
        for (db in c('adhd_genes', sprintf('%s_developmental', region))) {
            cat(md, phen, db, '\n')
            project_name = sprintf('logP_%s_%s_%s', md, phen, db)
            # make sure no dots in the name
            project_name = gsub(x=project_name, pattern='\\.',
                                replacement='_')
            db_file = sprintf('~/data/post_mortem/%s.gmt', db)
            enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabaseFile=db_file,
                                        enrichDatabaseType="genesymbol",
                                        interestGene=tmp2,
                                        interestGeneType="genesymbol",
                                        sigMethod="top", topThr=150000,
                                        minNum=3,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, perNum=10000,
                                        outputDirectory = data_dir,
                                        projectName=project_name,
                                        maxNum=800)
            out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                                md, phen, db)
            out_fname = gsub(x=out_fname, pattern='\\.',
                                replacement='_')
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
        DBs = c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant')
        for (db in DBs) {
            cat(region, db, '\n')
            project_name = sprintf('logP_%s_%s_%s', md, phen, db)
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
            out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                                md, phen, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }
}
```

And then for splices:

```r
# bw
library(WebGestaltR)
ncpu=8
data_dir = '~/data/post_mortem/'

r = 'acc'
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
    tmp2 = data.frame(hgnc_symbol=names(ranks), rank=ranks)
    for (db in c('adhd_genes', sprintf('%s_developmental', r))) {
        cat(r, p, db, '\n')
        project_name = sprintf('logPiso_%s_%s_%s', r, p, db)
        # make sure no dots in the name
        project_name = gsub(x=project_name, pattern='\\.',
                            replacement='_')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabaseFile=db_file,
                                    enrichDatabaseType="genesymbol",
                                    interestGene=tmp2,
                                    interestGeneType="genesymbol",
                                    sigMethod="top", topThr=150000,
                                    minNum=3,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000,
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    maxNum=800)
        out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                            r, p, db)
        out_fname = gsub(x=out_fname, pattern='\\.',
                            replacement='_')
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(r, db, '\n')
        project_name = sprintf('logPiso_%s_%s_%s', r, p, db)
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
        out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                            r, p, db)
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
}
```
