# 2021-03-08 10:00:31

Let's rerun some our our main results, this time using the all subtype results:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=30

load('~/data/post_mortem/DGE_03022021.RData')

for (region in c('acc', 'cau')) {
    res_str = sprintf('res = dge_%s[["all"]][["res"]]', region)
    eval(parse(text=res_str))
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)

    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    for (db in c('KEGG', 'Panther', 'Reactome', 'Wikipathway')) {
        cat(region, db, '\n')
        project_name = sprintf('WGP1_all_%s_%s_10K', region, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=sprintf('pathway_%s', db),
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", minNum=3,
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, topThr=20, perNum=10000))
    }
}
```

Let's also run the other DBs, not just the pathways:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=30

load('~/data/post_mortem/DGE_03022021.RData')

for (region in c('caudate', 'acc')) {
    res_str = ifelse(region == 'acc', 'dge_acc[["all"]]$res',
                     'dge_cau[["all"]]$res')
    ranks_str = sprintf('ranks = -log(%s$pvalue) * sign(%s$log2FoldChange)',
                        res_str, res_str)
    gid_str = sprintf('geneid=substring(rownames(%s), 1, 15)', res_str)
    
    eval(parse(text=ranks_str))
    eval(parse(text=gid_str))

    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    DBs = c(sprintf('my_%s_sets', region))
    for (db in DBs) {
        cat(res_str, db, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        project_name = sprintf('WG7_%s_%s_10K', res_str, db)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                            organism="hsapiens",
                            enrichDatabaseFile=db_file,
                            enrichDatabaseType="genesymbol",
                            interestGene=tmp2,
                            outputDirectory = data_dir,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=20,
                            minNum=3, projectName=project_name,
                            isOutput=T, isParallel=T,
                            nThreads=ncpu, perNum=10000, maxNum=800))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG7_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(res_str, db, '\n')
        project_name = sprintf('WG7_%s_%s_10K', res_str, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=db,
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", topThr=20,
                                    outputDirectory = data_dir,
                                    minNum=5, projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG7_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }
}
```

I also created the 03082021 version of DTE and DTU by loading the most recent
one and adding the 'all' version. That what I'll load to run the GSEA:

```r
library(WebGestaltR)
library(dplyr)

data_dir = '~/data/post_mortem/'
ncpu=30

load('~/data/post_mortem/DTE_03082021.RData')

for (myregion in c('caudate', 'acc')) {
    res_str = ifelse(myregion == 'acc', 'res = dte_acc[["all"]]$res',
                     'res = dte_cau[["all"]]$res')
    eval(parse(text=res_str))

    res_str = ifelse(myregion == 'acc', 'tmp = dte_acc[["all"]]$stageRObj',
                     'tmp = dte_cau[["all"]]$stageRObj')
    eval(parse(text=res_str))
    tx_meta = getTx2gene(tmp)

    res$TXNAME = substr(rownames(res), 1, 15)
    res$rank = -log(res$pvalue) * sign(res$log2FoldChange)
    m = merge(as.data.frame(res), tx_meta, by='TXNAME')

    ranks = m %>% group_by(GENEID) %>% slice_min(n=1, pvalue, with_ties=F)
    tmp2 = data.frame(geneid=ranks$GENEID, rank=ranks$rank)
    tmp2 = tmp2[order(tmp2$rank, decreasing=T),]

    res_str = sprintf('DTE_all_%s', myregion)
    DBs = c(sprintf('my_%s_sets', myregion))
    for (db in DBs) {
        cat(res_str, db, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        project_name = sprintf('WG8_%s_%s_10K', res_str, db)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                            organism="hsapiens",
                            enrichDatabaseFile=db_file,
                            enrichDatabaseType="genesymbol",
                            interestGene=tmp2,
                            outputDirectory = data_dir,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=20,
                            minNum=3, projectName=project_name,
                            isOutput=T, isParallel=T,
                            nThreads=ncpu, perNum=10000, maxNum=800))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG8_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(res_str, db, '\n')
        project_name = sprintf('WG8_%s_%s_10K', res_str, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=db,
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", topThr=20,
                                    outputDirectory = data_dir,
                                    minNum=5, projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG8_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('KEGG', 'Panther', 'Reactome', 'Wikipathway')
    for (db in DBs) {
        cat(myregion, db, '\n')
        project_name = sprintf('WGP2_all_%s_%s_10K', myregion, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=sprintf('pathway_%s', db),
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", minNum=3,
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, topThr=20, perNum=10000))
    }
}
```

And finally we run GSEA for DTU:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=2

load('~/data/post_mortem/DTU_03082021.RData')

for (myregion in c('caudate', 'acc')) {
    res_str = ifelse(myregion == 'acc', 'res = dtu_acc[["all"]]$res.g',
                     'res = dtu_cau[["all"]]$res.g')
    eval(parse(text=res_str))

    tmp2 = data.frame(geneid=res$gene_id, rank=-log(res$pvalue))
    tmp2 = tmp2[order(tmp2$rank, decreasing=T),]

    res_str = sprintf('DTU_all_%s', myregion)
    DBs = c(sprintf('my_%s_sets', myregion))
    for (db in DBs) {
        cat(res_str, db, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        project_name = sprintf('WG9_%s_%s_10K', res_str, db)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                            organism="hsapiens",
                            enrichDatabaseFile=db_file,
                            enrichDatabaseType="genesymbol",
                            interestGene=tmp2,
                            outputDirectory = data_dir,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=20,
                            minNum=3, projectName=project_name,
                            isOutput=T, isParallel=T,
                            nThreads=ncpu, perNum=10000, maxNum=800))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG9_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(res_str, db, '\n')
        project_name = sprintf('WG9_%s_%s_10K', res_str, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=db,
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", topThr=20,
                                    outputDirectory = data_dir,
                                    minNum=5, projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG9_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('KEGG', 'Panther', 'Reactome', 'Wikipathway')
    for (db in DBs) {
        cat(myregion, db, '\n')
        project_name = sprintf('WGP3_all_%s_%s_10K', myregion, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=sprintf('pathway_%s', db),
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", minNum=3,
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, topThr=20, perNum=10000))
    }
}
```

The thing is that for DTU we could also use the stageR screening pvalues:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=2

load('~/data/post_mortem/DTU_03082021.RData')

for (myregion in c('caudate', 'acc')) {
    res_str = ifelse(myregion == 'acc', 'res = dtu_acc[["all"]]',
                     'res = dtu_cau[["all"]]')
    eval(parse(text=res_str))

    df = stageR::getAdjustedPValues(res$stageRObj, 
                                    onlySignificantGenes=FALSE, order=TRUE)
    df2 = df[, c('geneID', 'gene')]
    df3 = df2[!duplicated(df2$geneID), ]
 
    tmp2 = data.frame(geneid=df3$geneID, rank=-log(df3$gene))
    tmp2 = tmp2[order(tmp2$rank, decreasing=T),]

    res_str = sprintf('DTUstageR_all_%s', myregion)
    DBs = c(sprintf('my_%s_sets', myregion))
    for (db in DBs) {
        cat(res_str, db, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        project_name = sprintf('WG9_%s_%s_10K', res_str, db)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                            organism="hsapiens",
                            enrichDatabaseFile=db_file,
                            enrichDatabaseType="genesymbol",
                            interestGene=tmp2,
                            outputDirectory = data_dir,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=20,
                            minNum=3, projectName=project_name,
                            isOutput=T, isParallel=T,
                            nThreads=ncpu, perNum=10000, maxNum=800))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG9_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(res_str, db, '\n')
        project_name = sprintf('WG9_%s_%s_10K', res_str, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=db,
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", topThr=20,
                                    outputDirectory = data_dir,
                                    minNum=5, projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG9_%s_%s_10K.csv', data_dir, res_str, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }

    DBs = c('KEGG', 'Panther', 'Reactome', 'Wikipathway')
    for (db in DBs) {
        cat(myregion, db, '\n')
        project_name = sprintf('WGP3stageR_all_%s_%s_10K', myregion, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=sprintf('pathway_%s', db),
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", minNum=3,
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, topThr=20, perNum=10000))
    }
}
```

## PRS

Time to re-run all the PRS analysis with the all data (instead of subtypes). I'm
using the code form note 181, which I updated to use NA as subtype, use all
possible covariates in the nuisancePC check, and to receive the data matrix
instead of myregion, which was never used:

```r
load('~/data/post_mortem/DGE_PRS_01272021.RData')
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
dgePRS_acc_all = list()
for (prs in prs_names) {
    dgePRS_acc_all[[prs]] = run_DGE_PRS(count_matrix, tx_meta, data, NA, prs,
                                        .05)
}

# reload data

dgePRS_cau_all = list()
for (prs in prs_names) {
    dgePRS_cau_all[[prs]] = run_DGE_PRS(count_matrix, tx_meta, data, NA, prs,
                                        .05)
}

save(dgePRS_acc_lnc, dgePRS_acc_pc, dgePRS_acc_pg,
     dgePRS_cau_lnc, dgePRS_cau_pc, dgePRS_cau_pg,
     dgePRS_acc_all, dgePRS_cau_all,
     file='~/data/post_mortem/DGE_PRS_03082021.RData')
```

And we do the same for DTE. The code is in note 86, and I had to update sending
the data matrix (instead of myregion, which was never used), added pH, RIN, and
brainbank, and all subtype selection. Then, we do:

```r
load('~/data/post_mortem/DTE_PRS_02032021.RData')
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
dtePRS_acc_all = list()
for (prs in prs_names) {
    dtePRS_acc_all[[prs]] = run_DTE_PRS(count_matrix, tx_meta, data, NA, prs,
                                        .05)
}

# reload data

dtePRS_cau_all = list()
for (prs in prs_names) {
    dtePRS_cau_all[[prs]] = run_DTE_PRS(count_matrix, tx_meta, data, NA, prs,
                                        .05)
}

save(dtePRS_acc_lnc, dtePRS_acc_pc, dtePRS_acc_pg,
     dtePRS_cau_lnc, dtePRS_cau_pc, dtePRS_cau_pg,
     dtePRS_acc_all, dtePRS_cau_all,
     file='~/data/post_mortem/DTE_PRS_03082021.RData')
```

Let's then compute the PRS and DX overlaps:

```r
library(GeneOverlap)
load('~/data/post_mortem/DGE_PRS_03082021.RData')
load('~/data/post_mortem/DGE_03022021.RData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc', 'pg')) {
    # res.dx = dge_cau[[subtypes[[st]]]]$res
    res.dx = dge_acc[[subtypes[[st]]]]$res
    for (p in prs_names) {
        cat(st, p, '\n')
        # res_str = sprintf('res.prs = dgePRS_cau_%s[["%s"]]', st, p)
        res_str = sprintf('res.prs = dgePRS_acc_%s[["%s"]]', st, p)
        eval(parse(text=res_str))

        both_res = merge(as.data.frame(res.dx), as.data.frame(res.prs), by=0,
                         all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs > 0,
                                 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx > 0,
                                'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            allUp = union(both_res[both_res$stat.prs > 0, 'Row.names'],
                          both_res[both_res$stat.dx > 0, 'Row.names'])
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
            go.obj <- testGeneOverlap(go.obj)
            pval2 = getPval(go.obj)
            this_res = c(subtypes[[st]], p, t, 'up', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs < 0,
                                 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx < 0,
                                'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            allDown = union(both_res[both_res$stat.prs < 0, 'Row.names'],
                            both_res[both_res$stat.dx < 0, 'Row.names'])
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
            go.obj <- testGeneOverlap(go.obj)
            pval2 = getPval(go.obj)
            this_res = c(subtypes[[st]], p, t, 'down', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t, 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t, 'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            pval2 = NA
            this_res = c(subtypes[[st]], p, t, 'abs', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
    }
}
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
# out_fname = '~/data/post_mortem/allDGE_caudateUpDown_prs_overlap_results_03082021.csv'
out_fname = '~/data/post_mortem/allDGE_accUpDown_prs_overlap_results_03082021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

And also ACC and Caudate overlaps:

```r
library(GeneOverlap)
load('~/data/post_mortem/DGE_03022021.RData')

all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc', 'pg')) {
    res.acc = dge_acc[[subtypes[[st]]]]$res
    res.cau = dge_cau[[subtypes[[st]]]]$res
    
    both_res = merge(as.data.frame(res.acc), as.data.frame(res.cau), by=0,
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs > 0,
                                'Row.names']
        dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx > 0,
                            'Row.names']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allUp = union(both_res[both_res$stat.prs > 0, 'Row.names'],
                        both_res[both_res$stat.dx > 0, 'Row.names'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(subtypes[[st]], t, 'up', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs < 0,
                                'Row.names']
        dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx < 0,
                            'Row.names']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allDown = union(both_res[both_res$stat.prs < 0, 'Row.names'],
                        both_res[both_res$stat.dx < 0, 'Row.names'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(subtypes[[st]], t, 'down', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
    for (t in c(.05, .01, .005, .001)) {
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
    }
}
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDGE_upDown_overlap_results_03082021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

And we do the same thing, but using DTE instead. I'll start with ACC and Caudate
overlaps because PRS hasn't finished calculating yet.

```r
library(GeneOverlap)
load('~/data/post_mortem/DTE_03082021.RData')

all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc', 'pg')) {
    res.acc = dte_acc[[subtypes[[st]]]]$res
    res.cau = dte_cau[[subtypes[[st]]]]$res
    
    both_res = merge(as.data.frame(res.acc), as.data.frame(res.cau), by=0,
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs > 0,
                                'Row.names']
        dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx > 0,
                            'Row.names']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allUp = union(both_res[both_res$stat.prs > 0, 'Row.names'],
                        both_res[both_res$stat.dx > 0, 'Row.names'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(subtypes[[st]], t, 'up', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs < 0,
                                'Row.names']
        dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx < 0,
                            'Row.names']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allDown = union(both_res[both_res$stat.prs < 0, 'Row.names'],
                        both_res[both_res$stat.dx < 0, 'Row.names'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(subtypes[[st]], t, 'down', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
    for (t in c(.05, .01, .005, .001)) {
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
    }
}
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDTE_upDown_overlap_results_03082021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

And now we go for PRS results:

```r
library(GeneOverlap)
load('~/data/post_mortem/DTE_PRS_03082021.RData')
load('~/data/post_mortem/DTE_03082021.RData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc', 'pg')) {
    res.dx = dte_cau[[subtypes[[st]]]]$res
    # res.dx = dte_acc[[subtypes[[st]]]]$res
    for (p in prs_names) {
        cat(st, p, '\n')
        res_str = sprintf('res.prs = dtePRS_cau_%s[["%s"]]', st, p)
        # res_str = sprintf('res.prs = dtePRS_acc_%s[["%s"]]', st, p)
        eval(parse(text=res_str))

        both_res = merge(as.data.frame(res.dx), as.data.frame(res.prs), by=0,
                         all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs > 0,
                                 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx > 0,
                                'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            allUp = union(both_res[both_res$stat.prs > 0, 'Row.names'],
                          both_res[both_res$stat.dx > 0, 'Row.names'])
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
            go.obj <- testGeneOverlap(go.obj)
            pval2 = getPval(go.obj)
            this_res = c(subtypes[[st]], p, t, 'up', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs < 0,
                                 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx < 0,
                                'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            allDown = union(both_res[both_res$stat.prs < 0, 'Row.names'],
                            both_res[both_res$stat.dx < 0, 'Row.names'])
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
            go.obj <- testGeneOverlap(go.obj)
            pval2 = getPval(go.obj)
            this_res = c(subtypes[[st]], p, t, 'down', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t, 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t, 'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            pval2 = NA
            this_res = c(subtypes[[st]], p, t, 'abs', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
    }
}
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDTE_caudateUpDown_prs_overlap_results_03082021.csv'
# out_fname = '~/data/post_mortem/allDTE_accUpDown_prs_overlap_results_03082021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

# 2021-03-09 07:17:41

Let's do the overlap and PRS comparisons for DTU as well. Just like for GSEA,
I'll do it with the stageR adjusted pvalues, and then just the DRIMseq ones.
Also, note that we don't have up and down here (not part of the DTU analysis).

```r
library(GeneOverlap)
load('~/data/post_mortem/DTU_03082021.RData')

all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc')) {
    res.acc = dtu_acc[[subtypes[[st]]]]$res.g
    res.cau = dtu_cau[[subtypes[[st]]]]$res.g
    
    both_res = merge(as.data.frame(res.acc), as.data.frame(res.cau), by='gene_id',
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t, 'gene_id']
        dx_genes = both_res[both_res$pvalue.dx < t, 'gene_id']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        pval2 = NA
        this_res = c(subtypes[[st]], t, 'abs', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDTU_upDown_overlap_results_03092021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

And repeat it for the stageR p-values:

```r
library(GeneOverlap)
load('~/data/post_mortem/DTU_03082021.RData')

all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc')) {
    res.acc = dtu_acc[[subtypes[[st]]]]
    res.cau = dtu_cau[[subtypes[[st]]]]

    df = stageR::getAdjustedPValues(res.acc$stageRObj, 
                                    onlySignificantGenes=FALSE, order=TRUE)
    df2 = df[, c('geneID', 'gene')]
    df.acc = df2[!duplicated(df2$geneID), ]

    df = stageR::getAdjustedPValues(res.cau$stageRObj, 
                                    onlySignificantGenes=FALSE, order=TRUE)
    df2 = df[, c('geneID', 'gene')]
    df.cau = df2[!duplicated(df2$geneID), ]
    
    both_res = merge(as.data.frame(df.acc), as.data.frame(df.cau), by='geneID',
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$gene.prs < t, 'geneID']
        dx_genes = both_res[both_res$gene.dx < t, 'geneID']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        pval2 = NA
        this_res = c(subtypes[[st]], t, 'abs', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDTUstageR_upDown_overlap_results_03092021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

And finally we do the same analysis for PRS and DTU. I only ran this analysis
for "all", because it took forever.

```r
library(GeneOverlap)
load('~/data/post_mortem/DTU_PRS_03082021.RData')
load('~/data/post_mortem/DTU_03082021.RData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
st = 'all'
res.dx = dtu_cau[[subtypes[[st]]]]$res.g
# res.dx = dtu_acc[[subtypes[[st]]]]$res.g
for (p in prs_names) {
    cat(st, p, '\n')
    res_str = sprintf('res.prs = dtuPRS_cau[["%s"]]', p)
    # res_str = sprintf('res.prs = dtuPRS_acc[["%s"]]', p)
    eval(parse(text=res_str))

    both_res = merge(res.dx, res.prs, by='gene_id',
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t, 'gene_id']
        dx_genes = both_res[both_res$pvalue.dx < t, 'gene_id']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        pval2 = NA
        this_res = c(subtypes[[st]], p, t, 'abs', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDTU_caudateUpDown_prs_overlap_results_03092021.csv'
# out_fname = '~/data/post_mortem/allDTU_accUpDown_prs_overlap_results_03092021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

## Making supplemental tables

Let's extract our all results into a table, adding gene names and subtype:

```r
load('~/data/post_mortem/DGE_03022021.RData')
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
mydir = '~/data/post_mortem/'

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]

for (r in c('acc', 'cau')) {
    for (st in c('all', 'protein_coding', 'lncRNA', 'pseudogene')) {
        res_str = sprintf('res = dge_%s[["%s"]]', r, st)
        eval(parse(text=res_str))
        fname = sprintf('%s/DGE_%s_%s_annot_03112021.csv', mydir, r, st)

        df = as.data.frame(res$res)
        colnames(df)[ncol(df)] = 'padj.FDR'
        df$GENEID = substr(rownames(df), 1, 15)
        df2 = merge(df, mart, sort=F,
                    by.x='GENEID', by.y='ensembl_gene_id', all.x=T, all.y=F)
        df2 = merge(df2, bt_slim, sort=F,
                    by.x='GENEID', by.y='gene_id', all.x=T, all.y=F)
        df2 = df2[order(df2$pvalue), ]
        df3 = as.data.frame(res$resIHW)
        df3$GENEID = substr(rownames(df3), 1, 15)
        df2 = merge(df2, df3[, c('GENEID', 'padj')], sort=F, by='GENEID',
                    all.x=T, all.y=F)
        colnames(df2)[ncol(df2)] = 'padj.IHW'
        write.csv(df2, row.names=F, file=fname)
    }
}

load('~/data/post_mortem/DTE_03082021.RData')
for (r in c('acc', 'cau')) {
    for (st in c('all', 'protein_coding', 'lncRNA', 'pseudogene')) {
        res_str = sprintf('res = dte_%s[["%s"]]', r, st)
        eval(parse(text=res_str))
        fname = sprintf('%s/DTE_%s_%s_annot_03112021.csv', mydir, r, st)

        df = res$res
        colnames(df)[ncol(df)] = 'padj.FDR'
        df$TXNAME = substr(rownames(df), 1, 15)
        tx2gene = stageR::getTx2gene(res$stageRObj)
        df2 = merge(as.data.frame(df), tx2gene, sort=F,
                    by.x='TXNAME', all.x=T, all.y=F)
        df2 = merge(df2, mart, sort=F,
                    by.x='GENEID', by.y='ensembl_gene_id', all.x=T, all.y=F)
        df2 = merge(df2, bt_slim, sort=F,
                    by.x='GENEID', by.y='gene_id', all.x=T, all.y=F)
        df2 = df2[order(df2$pvalue), ]
        df3 = as.data.frame(res$resIHW)
        df3$TXNAME = substr(rownames(df3), 1, 15)
        df2 = merge(df2, df3[, c('TXNAME', 'padj')], sort=F, by='TXNAME',
                    all.x=T, all.y=F)
        colnames(df2)[ncol(df2)] = 'padj.IHW'
        write.csv(df2, row.names=F, file=fname)
    }
}

load('~/data/post_mortem/DTU_03082021.RData')
for (r in c('acc', 'cau')) {
    for (st in c('all', 'protein_coding', 'lncRNA')) {
        res_str = sprintf('res = dtu_%s[["%s"]]', r, st)
        eval(parse(text=res_str))
        fname = sprintf('%s/DTU_%s_%s_annot_03112021.csv', mydir, r, st)

        df = stageR::getAdjustedPValues(res$stageRObj,
                                        onlySignificantGenes=FALSE, order=TRUE)
        colnames(df)[3:4] = c('padj_gene', 'padj_transcript')
        df2 = merge(df, mart, sort=F,
                    by.x='geneID', by.y='ensembl_gene_id', all.x=T, all.y=F)
        df2 = merge(df2, bt_slim, sort=F,
                    by.x='geneID', by.y='gene_id', all.x=T, all.y=F)
        write.csv(df2, row.names=F, file=fname)
    }
}
```

Finally, let's write out our methylation results.

```r
for (r in c('acc', 'cau')) {
    region = ifelse(r=='acc', 'ACC', 'Caudate')
    load(sprintf('~/data/methylation_post_mortem/res_%s_02262021.RData',
                 region))
    for (st in c("all", "Island", "Shelf", "Shore", "Sea", "enhancer", "body",
                 "promoter1", "promoter2")) {
        res_str = sprintf('res = res_%s[["%s"]]', r, st)
        eval(parse(text=res_str))
        fname = sprintf('%s/methylDMP_%s_%s_annot_03082021.csv', mydir, r, st)

        df = as.data.frame(res$DMPs)
        df2 = df[order(df$P.Value), ]
        write.csv(df2, row.names=F, file=fname)
    }
}
```

I already exported the results for methyl and PRS in the past.


# TODO:
 * spit out methylation overlaps