# 2020-11-05 07:19:41

Let me see if the results change much with using cameraPR instead of WebGestalt.

```r
library(limma)
library(WebGestaltR)  # for readGmt()

get_enrich_order2 = function( res, gene_sets ){
  if( !is.null(res$z.std) ){
    stat = res$z.std
  }else if( !is.null(res$F.std) ){
    stat = res$F.std
  }else if( !is.null(res$t) ){
    stat = res$t
  }else{
    stat = res$F
  }
  names(stat) = res$hgnc_symbol
  stat = stat[!is.na(names(stat))]
  index = ids2indices(gene_sets, names(stat))
  cameraPR( stat, index )
}

data_dir = '~/data/rnaseq_derek/'
load(sprintf('%s/xmodal_results_10152020.RData', data_dir))

region='acc'
# region='caudate'
eval(parse(text=sprintf('res = rnaseq_%s', region)))

tmp2 = res[, c('hgnc_symbol', 't')]
for (db in c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant',
                'pathway_KEGG', 'disease_Disgenet',
                'phenotype_Human_Phenotype_Ontology',
                'network_PPI_BIOGRID')) {
    cat(region, db, '\n')
    # assigning hgnc to our gene set lists
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
    res_camera = get_enrich_order2( tmp2, mylist )
    if (is.null(gs$geneSetDes)) {
        # PPI doesn't have descriptions
        m = cbind(rownames(res_camera), res_camera, NA)
        colnames(m)[1] = 'Row.names'
        colnames(m)[ncol(m)] = 'description'
    } else {
        # attach gene set description
        m = merge(res_camera, gs$geneSetDes, by.x=0, by.y=1)
        m = m[order(m$PValue), ]
    }
    out_fname = sprintf('%s/camera_%s_%s.csv', data_dir, region, db)
    write.csv(m, file=out_fname, quote=F, row.names=F)
}
# my own GMTs
for (db in c('disorders', sprintf('%s_developmental', region))) {
    cat(region, db, '\n')
    db_file = sprintf('~/data/post_mortem/%s.gmt', db)
    gmt = readGmt(db_file) # already in gene symbols
    # and convert it to lists
    mylist = list()
    for (s in unique(gmt$geneSet)) {
        mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
    }
    res_camera = get_enrich_order2( tmp2, mylist )
    # some massaging to look like the other results
    desc = gmt[, c(1,2)]
    desc = desc[!duplicated(desc$geneSet), ]
    m = merge(res_camera, desc, by.x=0, by.y=1)
    m = m[order(m$PValue), ]
    out_fname = sprintf('%s/camera_%s_%s.csv', data_dir, region, db)
    write.csv(m, file=out_fname, quote=F, row.names=F)
}
```

Now we do the same thing for methylation:

```r
data_dir = '~/data/methylation_post_mortem/'

# region='acc'
region='caudate'
res = readRDS(sprintf('%s/%s_methyl_results_11032020.rds', data_dir, region))
idx = res$gene != ''
genes = res[idx, ]

imautosome = which(genes$CHR != 'X' &
                    genes$CHR != 'Y' &
                    genes$CHR != 'MT')
genes = genes[imautosome, ]

tmp = genes[, c('gene', 't')]
tmp2 = c()
# will do it this way because of the two tails of T distribution
for (g in unique(tmp$gene)) {
    gene_data = tmp[tmp$gene==g, ]
    best_res = which.max(abs(gene_data$t))
    tmp2 = rbind(tmp2, gene_data[best_res, ])
}
colnames(tmp2)[1] = 'hgnc_symbol'
for (db in c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant',
                'pathway_KEGG', 'disease_Disgenet',
                'phenotype_Human_Phenotype_Ontology',
                'network_PPI_BIOGRID')) {
    cat(region, db, '\n')
    # assigning hgnc to our gene set lists
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
    res_camera = get_enrich_order2( tmp2, mylist )
    if (is.null(gs$geneSetDes)) {
        # PPI doesn't have descriptions
        m = cbind(rownames(res_camera), res_camera, NA)
        colnames(m)[1] = 'Row.names'
        colnames(m)[ncol(m)] = 'description'
    } else {
        # attach gene set description
        m = merge(res_camera, gs$geneSetDes, by.x=0, by.y=1)
        m = m[order(m$PValue), ]
    }
    out_fname = sprintf('%s/camera_%s_%s.csv', data_dir, region, db)
    write.csv(m, file=out_fname, quote=F, row.names=F)
}
# my own GMTs
for (db in c('disorders', sprintf('%s_developmental', region))) {
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
                                sigMethod="top", topThr=10,
                                minNum=3, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000)
    out_fname = sprintf('%s/WG_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, quote=F,
                row.names=F)
}
```



# TODO
