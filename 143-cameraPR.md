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

# region='acc'
region='caudate'
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


```r

load('~/data/rnaseq_derek/adhd_genesets_philip.RDATA')
load('~/data/rnaseq_derek/c5_gene_sets.RData')
load('~/data/rnaseq_derek/brain_disorders_gene_sets.RData')
load('~/data/rnaseq_derek/data_for_alex.RData')
co = .9 
idx = anno$age_category==1 & anno$cutoff==co
genes_overlap = unique(anno[idx, 'anno_gene'])
for (s in 2:5) {
  idx = anno$age_category==s & anno$cutoff==co
  g2 = unique(anno[idx, 'anno_gene'])
  genes_overlap = intersect(genes_overlap, g2)
}
genes_unique = list()
for (s in 1:5) {
  others = setdiff(1:5, s)
  idx = anno$age_category==s & anno$cutoff==co
  g = unique(anno[idx, 'anno_gene'])
  for (s2 in others) {
    idx = anno$age_category==s2 & anno$cutoff==co
    g2 = unique(anno[idx, 'anno_gene'])
    rm_me = g %in% g2
    g = g[!rm_me]
  }
  genes_unique[[sprintf('dev%s_c%.1f', s, co)]] = unique(g)
}
genes_unique[['overlap']] = unique(genes_overlap)

data2 = cbind(data, mydata)
form = ~ Diagnosis + PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
res = topTable(fit2, coef='DiagnosisControl', number=Inf)

adhd_camera = get_enrich_order2( res, t2 ) 
c5_camera = get_enrich_order2( res, c5_all)
dis_camera = get_enrich_order2( res, disorders)
dev_camera = get_enrich_order2( res, genes_unique )
```



# TODO
