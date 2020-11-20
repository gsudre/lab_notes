# 2020-11-18 15:02:03

The idea here is to run the same gene set enrichment analysis, but using Kwangmi's data. Let's start with camera, because it runs faster, and then we can set up WebGestalt.

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

data_dir = '~/data/kwangmi/clinic/'
phenotypes = c('maxOverTimeSX_hi', 'everADHD_DX')
stat = c('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')
transforms = c('neg', 'negLogP')
DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant',
        'pathway_KEGG', 'disease_Disgenet',
        'phenotype_Human_Phenotype_Ontology',
        'network_PPI_BIOGRID',
        'geneontology_Biological_Process',
        'geneontology_Cellular_Component',
        'geneontology_Molecular_Function')

for (phen in phenotypes) {
    df = read.table(sprintf('%s/white_%s.SAIGE.gene.0.1.txt',
                            data_dir, phen), header=1)
    for (st in stat) {
        for (tr in transforms) {
            tmp2 = df[, c('Gene', st)]
            tmp2 = tmp2[!is.na(tmp2[, st]), ]
            if (tr == 'neg') {
                tmp2[, st] = 1-tmp2[, st]
            } else {
                tmp2[, st] = 1-log(tmp2[, st])
            }
            colnames(tmp2) = c('hgnc_symbol', 't')
            # my own GMTs
            for (db in c('disorders', 'adhd_genes')) {
                cat(phen, db, st, tr, '\n')
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
                out_fname = sprintf('%s/camera_%s_%s_%s_%s.csv', data_dir, phen, db, st, tr)
                write.csv(m, file=out_fname, quote=F, row.names=F)
            }
            for (db in DBs) {
                cat(phen, db, st, tr, '\n')
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
                    m$description = gsub(x=m$description, pattern=',', replacement=';')
                }
                out_fname = sprintf('%s/camera_%s_%s_%s_%s.csv', data_dir, phen, db, st, tr)
                write.csv(m, file=out_fname, quote=F, row.names=F)
            }
        }
    }
}
```

And we run the same stuff for WebGestalt:

```r
library(WebGestaltR)
ncpu=2
nperms=100
data_dir = '~/data/kwangmi/clinic/'
phenotypes = c('maxOverTimeSX_hi', 'everADHD_DX')
stat = c('Pvalue', 'Pvalue_Burden', 'Pvalue_SKAT')
transforms = c('neg', 'negLogP')
DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant',
        'pathway_KEGG', 'disease_Disgenet',
        'phenotype_Human_Phenotype_Ontology',
        'network_PPI_BIOGRID')

for (phen in phenotypes) {
    df = read.table(sprintf('%s/white_%s.SAIGE.gene.0.1.txt',
                            data_dir, phen), header=1)
    for (st in stat) {
        for (tr in transforms) {
            tmp2 = df[, c('Gene', st)]
            tmp2 = tmp2[!is.na(tmp2[, st]), ]
            if (tr == 'neg') {
                tmp2[, st] = 1-tmp2[, st]
            } else {
                tmp2[, st] = 1-log(tmp2[, st])
            }
            colnames(tmp2) = c('hgnc_symbol', 't')
            # my own GMTs
            for (db in c('disorders', 'adhd_genes')) {
                cat(phen, db, st, tr, '\n')
                db_file = sprintf('~/data/post_mortem/%s.gmt', db)
                project_name = sprintf('%s_%s_%s_%s', phen, db, st, tr)
                enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabaseFile=db_file,
                                            enrichDatabaseType="genesymbol",
                                            interestGene=tmp2,
                                            interestGeneType="genesymbol",
                                            sigMethod="top", topThr=150000,
                                            minNum=3,
                                            isOutput=T, isParallel=T,
                                            nThreads=ncpu, perNum=nperms,
                                            outputDirectory = data_dir,
                                            projectName=project_name)
                out_fname = sprintf('%s/WG_%s_%s_%s_%s.csv', data_dir, phen, db, st, tr)
                write.csv(enrichResult, file=out_fname, quote=F, row.names=F)
            }
            for (db in DBs) {
                cat(phen, db, st, tr, '\n')
                project_name = sprintf('%s_%s_%s_%s', phen, db, st, tr)
                enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabase=db,
                                            interestGene=tmp2,
                                            interestGeneType="genesymbol",
                                            sigMethod="top", topThr=150000,
                                            outputDirectory = data_dir,
                                            minNum=5, projectName=project_name,
                                            isOutput=T, isParallel=T,
                                            nThreads=ncpu, perNum=nperms)
                out_fname = sprintf('%s/WG_%s_%s_%s_%s.csv', data_dir, phen, db, st, tr)
                write.csv(enrichResult, file=out_fname, quote=F, row.names=F)
            }
        }
    }
}
```

# 2020-11-19 11:12:55

For gene overlap, we can do:

```r
library(GeneOverlap)
library(corrplot)
data_dir = '~/data/kwangmi/clinic/'
res_var = c('maxOverTimeSX_hi', 'everADHD_DX', 'maxOverTimeSX_inatt')
res_str = c('hi', 'DX', 'inatt')
thres=.005
pvals = matrix(nrow=length(res_var), ncol=length(res_var),
               dimnames=list(res_str, res_str))
st = 'Pvalue'
for (i in 1:length(res_var)) {​
    for (j in 1:length(res_var)) {​
        res1 = read.table(sprintf('%s/white_%s.SAIGE.gene.0.1.txt',
                                  data_dir, res_var[i]),
                          header=1)[, c('Gene', st)]
        res2 = read.table(sprintf('%s/white_%s.SAIGE.gene.0.1.txt',
                                  data_dir, res_var[j]),
                          header=1)[, c('Gene', st)]
        res1 = res1[!is.na(res1[, st]), ]
        res2 = res2[!is.na(res2[, st]), ]
        colnames(res1) = c('hgnc_symbol', 'pvalue')
        colnames(res2) = c('hgnc_symbol', 'pvalue')
        uni = intersect(res1$hgnc_symbol, res2$hgnc_symbol)
        # only evaluate genes in both sets
        res1 = res1[res1$hgnc_symbol %in% uni, ]
        res2 = res2[res2$hgnc_symbol %in% uni, ]
        go.obj <- newGeneOverlap(res1$hgnc_symbol[res1$pvalue < thres],
                                 res2$hgnc_symbol[res2$pvalue < thres],
                                 genome.size=length(uni))
        go.obj <- testGeneOverlap(go.obj)
        pvals[res_str[i], res_str[j]] = getPval(go.obj)
    }​
}​
quartz()
b = 1-pvals
b[b<.95] = NA
corrplot(b, method='color', tl.cex=.8, cl.cex=1, type='upper',
         is.corr=F, na.label='x')
```