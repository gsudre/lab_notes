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
        out_fname = sprintf('%s/WGnoOverlap_%s_10K.csv', data_dir, region)
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
}
```

The results were very extreme... let's try a different cut-off.

Here's a good way of understanding the GSEA metrics, including the FDR:

https://www.gsea-msigdb.org/gsea/doc/GSEAUserGuideFrame.html?Interpreting_GSEA

I'm not sure if it's fair to add the GWAS and TWAS sets here... the FDR
correction shoiuld be made among sets in a similar domain, right?

Before we do something like that, let's see if ACC genes enrich caudate, or
vice-versa:

```r
res = rnaseq_acc
db_file = '~/data/post_mortem/my_caudate_sets.gmt'

ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

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
    out_fname = sprintf('%s/WG_pmACC_caudate_sets_10K.csv', data_dir, region)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

And then I inverted caudate and ACC. But the answer is yes, they were both
enriched. Likely because there is a big overlap. So, let's create a few lists
where there are genes specific to each region. First, let's check the overlaps:

```r
gmta = read.GMT('~/data/post_mortem/my_acc_sets.gmt')
gmtc = read.GMT('~/data/post_mortem/my_caudate_sets.gmt')
for (i in 1:5) {
    mydev = sprintf('dev%d_c0.9', i)
    inter = length(intersect(gmta[[mydev]]$genes, gmtc[[mydev]]$genes))
    nacc = length(gmta[[mydev]]$genes)
    ncau = length(gmtc[[mydev]]$genes)
    cat(mydev, nacc, ncau, inter, '\n')
}
mydev = 'overlap_c0.9'
inter = length(intersect(gmta[[mydev]]$genes, gmtc[[mydev]]$genes))
nacc = length(gmta[[mydev]]$genes)
ncau = length(gmtc[[mydev]]$genes)
cat(mydev, nacc, ncau, inter, '\n')
```

Let's make better sets then. Some including the other brain region and other
don't. So, basically:

 * ACC_dev1:5
 * ACC dev1:5_devSpec + overlap
 * ACC_dev1:5_regSpec
 * ACC_dev1:5_regSpec_devSpec + overlap_regSpec

```r
library(ActivePathways)
dev_names = c('prenatal', 'infant (0-2 yrs)', 'child (3-11 yrs)',
              'adolescent (12-19 yrs)', 'adult (>19 yrs)')
gmt = read.GMT('~/data/post_mortem/hsapiens_disease_Disgenet_entrezgene.gmt')

r = 'acc'
co = .95

struc_id = ifelse(r == 'acc', 'Allen:10278', 'Allen:10333')
other_id = ifelse(r != 'acc', 'Allen:10278', 'Allen:10333')
anno = get_annotated_genes(structure_ids=struc_id, dataset='5_stages', 
                            cutoff_quantiles=c(co))
anno2 = get_annotated_genes(structure_ids=other_id, dataset='5_stages', 
                            cutoff_quantiles=c(co))
other_genes = unique(anno2[, 'anno_gene'])

junk = gmt[1:22]
for (i in 1:5) {
    g = unique(anno[anno$age_category==i, 'anno_gene'])
    tmp = list(id = sprintf('dev%d_c%.3f', i, co), genes = unique(g),
             name = dev_names[i])
    cat(tmp$id, length(g), '\n')
    junk[[i]] = tmp
}
cnt = 6
for (i in 1:5) {
    g = unique(anno[anno$age_category==i, 'anno_gene'])
    g = g[! g %in% other_genes]
    tmp = list(id = sprintf('dev%d_c%.3f_regSpec', i, co), genes = unique(g),
             name = dev_names[i])
    cat(tmp$id, length(g), '\n')
    junk[[cnt]] = tmp
    cnt = cnt + 1
}
# now the unique sets
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
  genes_unique[[s]] = unique(g)
}
# writing them to GMT
for (i in 1:5) {
    g = genes_unique[[i]]
    tmp = list(id = sprintf('dev%d_c%.3f_devSpec', i, co), genes = unique(g),
             name = dev_names[i])
    cat(tmp$id, length(g), '\n')
    junk[[cnt]] = tmp
    cnt = cnt + 1
}
g = genes_overlap
tmp = list(id = sprintf('overlap_c%.3f', co), genes = unique(g),
             name = 'overlap')
cat(tmp$id, length(g), '\n')
junk[[cnt]] = tmp
cnt = cnt + 1
for (i in 1:5) {
    g = genes_unique[[i]]
    g = g[! g %in% other_genes]
    tmp = list(id = sprintf('dev%d_c%.3f_devSpec_regSpec', i, co),
                genes = unique(g), name = dev_names[i])
    cat(tmp$id, length(g), '\n')
    junk[[cnt]] = tmp
    cnt = cnt + 1
}
g = genes_overlap
g = g[! g %in% other_genes]
tmp = list(id = sprintf('overlap_c%.3f_regSpec', co), genes = unique(g),
             name = 'overlap')
cat(tmp$id, length(g), '\n')
junk[[cnt]] = tmp
gmt_name = sprintf('~/data/post_mortem/%s_manySets_co%.3f.gmt', r, co)
write.GMT(junk, gmt_name)
```

These are my ACC set sizes for co=.9:

```
dev1_c0.90 1758 
dev2_c0.90 1625 
dev3_c0.90 1531 
dev4_c0.90 1639 
dev5_c0.90 1639 
dev1_c0.90_regSpec 150 
dev2_c0.90_regSpec 95 
dev3_c0.90_regSpec 54 
dev4_c0.90_regSpec 84 
dev5_c0.90_regSpec 149 
dev1_c0.90_devSpec 471 
dev2_c0.90_devSpec 38 
dev3_c0.90_devSpec 41 
dev4_c0.90_devSpec 39 
dev5_c0.90_devSpec 88 
overlap_c0.90 927 
dev1_c0.90_devSpec_regSpec 125 
dev2_c0.90_devSpec_regSpec 20 
dev3_c0.90_devSpec_regSpec 4 
dev4_c0.90_devSpec_regSpec 6 
dev5_c0.90_devSpec_regSpec 55 
overlap_c0.90_regSpec 7 
```

and this is for the Caudate:

```
dev1_c0.90 1733 
dev2_c0.90 1805 
dev3_c0.90 1724 
dev4_c0.90 1601 
dev5_c0.90 1701 
dev1_c0.90_regSpec 88 
dev2_c0.90_regSpec 167 
dev3_c0.90_regSpec 219 
dev4_c0.90_regSpec 165 
dev5_c0.90_regSpec 170 
dev1_c0.90_devSpec 388 
dev2_c0.90_devSpec 95 
dev3_c0.90_devSpec 114 
dev4_c0.90_devSpec 34 
dev5_c0.90_devSpec 27 
overlap_c0.90 926 
dev1_c0.90_devSpec_regSpec 63 
dev2_c0.90_devSpec_regSpec 53 
dev3_c0.90_devSpec_regSpec 90 
dev4_c0.90_devSpec_regSpec 21 
dev5_c0.90_devSpec_regSpec 19 
overlap_c0.90_regSpec 11 
```

Again, devSpec means they are developmental stage-specific (overlap removed),
and regSpec means the genes for the other brain region were removed.

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
    db_file = sprintf('~/data/post_mortem/%s_manySets_co%.3f.gmt', region, co)
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
                        nThreads=ncpu, perNum=10000, maxNum=2000))
    if (class(enrichResult) != "try-error") {
        out_fname = sprintf('%s/WGmanySets_%s_co%.3f_10K.csv', data_dir,
                            region, co)
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
}
```

# 2021-01-06 10:05:56

We got some mixed patterns for caudate... let me see if running co=.95 makes the
pattern more clear. Here are the set sizes, ACC first:

```
dev1_c0.950 856 
dev2_c0.950 810 
dev3_c0.950 756 
dev4_c0.950 809 
dev5_c0.950 803 
dev1_c0.950_regSpec 91 
dev2_c0.950_regSpec 70 
dev3_c0.950_regSpec 34 
dev4_c0.950_regSpec 58 
dev5_c0.950_regSpec 100 
dev1_c0.950_devSpec 284 
dev2_c0.950_devSpec 31 
dev3_c0.950_devSpec 26 
dev4_c0.950_devSpec 26 
dev5_c0.950_devSpec 51 
overlap_c0.950 416 
dev1_c0.950_devSpec_regSpec 76 
dev2_c0.950_devSpec_regSpec 14 
dev3_c0.950_devSpec_regSpec 1 
dev4_c0.950_devSpec_regSpec 6 
dev5_c0.950_devSpec_regSpec 28 
overlap_c0.950_regSpec 3 

dev1_c0.950 873 
dev2_c0.950 898 
dev3_c0.950 876 
dev4_c0.950 783 
dev5_c0.950 877 
dev1_c0.950_regSpec 68 
dev2_c0.950_regSpec 102 
dev3_c0.950_regSpec 141 
dev4_c0.950_regSpec 113 
dev5_c0.950_regSpec 127 
dev1_c0.950_devSpec 239 
dev2_c0.950_devSpec 51 
dev3_c0.950_devSpec 63 
dev4_c0.950_devSpec 15 
dev5_c0.950_devSpec 18 
overlap_c0.950 400 
dev1_c0.950_devSpec_regSpec 50 
dev2_c0.950_devSpec_regSpec 31 
dev3_c0.950_devSpec_regSpec 38 
dev4_c0.950_devSpec_regSpec 6 
dev5_c0.950_devSpec_regSpec 14 
overlap_c0.950_regSpec 4 
```

