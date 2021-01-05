# 2021-01-05 13:01:48

Let's try specifying the brian development sets without the overlap:

```r
library('ABAEnrichment')
tmp = get_id('striatum')
tmp[tmp$ontology=='developmental',]
```

```
     structure      ontology structure_id
1 STR_striatum developmental  Allen:10333
```

```r
tmp = get_id(' cingulate')
tmp[tmp$ontology=='developmental',]
```

```
                                                    structure      ontology structure_id
1 MFC_anterior (rostral) cingulate (medial prefrontal) cortex developmental  Allen:10278
```

Let's then create the ACC set first:

```r
r = 'acc'
co = .95

struc_id = ifelse(r == 'acc', 'Allen:10278', 'Allen:10333') 
anno = get_annotated_genes(structure_ids=struc_id, dataset='5_stages', 
                            cutoff_quantiles=c(co))
```

Now we create the new GMT:

```r
library(ActivePathways)

dev_names = c('prenatal', 'infant (0-2 yrs)', 'child (3-11 yrs)',
              'adolescent (12-19 yrs)', 'adult (>19 yrs)')
gmt = read.GMT('~/data/post_mortem/hsapiens_disease_Disgenet_entrezgene.gmt')
junk = gmt[1:7]
for (i in 1:5) {
    g = unique(anno[anno$age_category==i, 'anno_gene'])
    cat(i, length(g), '\n')
    tmp = list(id = sprintf('dev%d_c%.2f', i, co), genes = unique(g),
             name = dev_names[i])
    junk[[i]] = tmp
}
# and we add the TWAS and GWAS sets too
gmt2 = read.GMT('~/data/post_mortem/adhd_genes.gmt')
junk[[6]] = gmt2[['GWAS1']]
junk[[7]] = gmt2[['TWAS']]
gmt_name = sprintf('~/data/post_mortem/%s_noOverlap.gmt', r)
write.GMT(junk, gmt_name)
```

And just run the same thing for Caudate. Now, a quick GSEA analysis:

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

ncpu=6

for (region in c('acc', 'caudate')) {
    eval(parse(text=sprintf('res = rnaseq_%s', region)))

    ranks = -log(res$P.Value) * sign(res$logFC)
    tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    # my own GMTs
    db_file = sprintf('~/data/post_mortem/%s_noOverlap.gmt', region)
    enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                        organism="hsapiens",
                        enrichDatabaseFile=db_file,
                        enrichDatabaseType="genesymbol",
                        interestGene=tmp2,
                        outputDirectory = data_dir,
                        interestGeneType="genesymbol",
                        sigMethod="top", topThr=150000,
                        minNum=3, projectName=project_name,
                        isOutput=F, isParallel=T,
                        nThreads=ncpu, perNum=10000, maxNum=1000))
    if (class(enrichResult) != "try-error") {
        out_fname = sprintf('%s/WGnoOverLap_%s_10K.csv', data_dir, region)
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
}
```

The results were very extreme... let's try a different cut-off.

Here's a good way of understanding the GSEA metrics, including the FDR:

https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?Interpreting_GSEA
