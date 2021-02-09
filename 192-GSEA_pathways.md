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

# TODO
 * anything on lncRNA and pseudogenes?
 * DTE
 * DTU
 * network analysis?


# Useful links:
 * https://www.nature.com/articles/s41596-018-0103-9