# 2021-02-08 19:17:01

One of the points people made in the HBCC presentation was exploring pathways
and protein interactions. It's possible to do that within GSEA and WG itself, so
let's see what we get there.

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=3

load('~/data/post_mortem/DGE_02082021.RData')

for (region in c('acc', 'cau')) {
    res_str = sprintf('res = dge_%s[["protein_coding"]][["res"]]',
                    region)
    eval(parse(text=res_str))
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)

    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    for (db in c('KEGG', 'Panther', 'Reactome', 'Wikipathway')) {
        cat(region, db, '\n')
        project_name = sprintf('WGP1_pc_%s_%s_10K', region, db)

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

OK, those are all running, so it's just a matter of scripting it for all DGE,
DTE, and DTU. Maybe they show differences?

For PPI BIOGRID we need a gene list only. I'll threshold at different nominal pvalues:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=2

load('~/data/post_mortem/DGE_02082021.RData')

region = 'acc'  # 'cau'
res_str = sprintf('res = dge_%s[["protein_coding"]][["res"]]',
                  region)
eval(parse(text=res_str))

geneid = substring(rownames(res), 1, 15)
db = 'network_PPI_BIOGRID'

for (p in c(.05, .01, .005) {
    tmp2 = geneid[which(res$pvalue < p)]
    cat(region, db, p, '\n')
    project_name = sprintf('WGP1_%s_%s_%.3f_1K', region, db, p)
    enrichResult <- try(WebGestaltR(enrichMethod="NTA",
                            organism="hsapiens",
                            enrichDatabase=db,
                            interestGene=tmp2,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", minNum=3,
                            outputDirectory = data_dir,
                            projectName=project_name,
                            isOutput=T, isParallel=T,
                            nThreads=ncpu,
                            highlightSeedNum=10,
                            networkConstructionMethod="Network_Retrieval_Prioritization"))
}
```

This also runs, but the interpretations might not be as straight forward.

Let's focus on re-creating the DTE and DTU results with the new variables (191),
and the methylation, before we go back to this network analysis.

# 2021-02-09 09:10:08

Let's run the same thing for DTE:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=3

load('~/data/post_mortem/DTE_02082021.RData')

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "TXNAME"), columns=c('GENEID','TXCHROM'),
               "TXNAME")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('transcript_id', 'transcript_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]

library(dplyr)
for (region in c('acc', 'cau')) {
    res_str = sprintf('res = dte_%s[["protein_coding"]][["res"]]',
                    region)
    eval(parse(text=res_str))
    res$TXNAME = substr(rownames(res), 1, 15)
    res$rank = -log(res$pvalue) * sign(res$log2FoldChange)

    # keep only the remaining transcripts and their corresponding genes
    txi_str = sprintf('counts = counts(dte_%s[["protein_coding"]][["dds"]])',
                    region)
    eval(parse(text=txi_str))
    txdf.sub = txdf[match(substr(rownames(counts), 1, 15), txdf$TXNAME),]
    
    tx_meta = merge(txdf.sub, bt_slim, by.x='TXNAME', by.y='transcript_id')
    m = merge(as.data.frame(res), tx_meta, by='TXNAME', all.x=T, all.y=F)

    ranks = m %>% group_by(GENEID) %>% slice_min(n=1, pvalue, with_ties=F)
    tmp2 = data.frame(geneid=ranks$GENEID, rank=ranks$rank)
    tmp2 = tmp2[order(tmp2$rank, decreasing=T),]

    for (db in c('KEGG', 'Panther', 'Reactome', 'Wikipathway')) {
        cat(region, db, '\n')
        project_name = sprintf('WGP2_pc_%s_%s_10K', region, db)

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

Not sure if lncRNA or pseudogenes would work here... let's see:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=2

load('~/data/post_mortem/DGE_02082021.RData')

for (region in c('acc', 'cau')) {
    res_str = sprintf('res = dge_%s[["pseudogene"]][["res"]]',
                    region)
    eval(parse(text=res_str))
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)

    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    for (db in c('KEGG', 'Panther', 'Reactome', 'Wikipathway')) {
        cat(region, db, '\n')
        project_name = sprintf('WGP1_pg_%s_%s_10K', region, db)

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

Everything failed for lncRNA and pseudogene, so now go there.


# TODO
 * DTU
 * network analysis?


# Useful links:
 * https://www.nature.com/articles/s41596-018-0103-9