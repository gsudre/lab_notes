# 2020-09-27 14:20:52

Let's run some gene set analysis in the isoform results Laura sent. The
spreadsheet is called PCAtest_results.xlsx and I downloaded it to
data/methylation. They ran the analysis independently for ACC and Caudate, and
they have results for RSEM counts and Kallisto, which according to Valer should
indeed have low correlation as they use different methods for counting. 

Philip said we should:

* run some pathway analyses to see what’s enriched by these genes (probably at
  the the differing levels of sig)
* run an overlap test for enrichment with the genes implicated in ADHD through
  GWAS/EWAS/RVAS and TWAS
* look at the really strong hits individually--- see if these genes are doing
  anything that sounds very ‘ADHD-y’

Let me run the first 2 and maybe have one of the IRTAs do the latter. I'll use
cameraPR for now because we don't need to specify a priori p-value cut-off, and
I'll make the pre-computed ranking just 1-pvalue. Also, I'll just use the
geometric mean of the p-value forboth test they reported (welch's t-test and
coin general independence test). 

```r
library(limma)
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
  # print(head(stat))
  index = ids2indices(gene_sets, names(stat))
  cameraPR( stat, index )
}
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

for (r in c('acc', 'caudate')) {
    for (m in c('rsem', 'kallisto')) {
        fname = sprintf('~/data/methylation/%s_%s.csv', r, m)
        df = read.csv(fname)
        colnames(df)[1] = 'hgnc_symbol'
        df$F = 1-df$pval
        adhd_cameraPR = get_enrich_order2( df, t2 ) 
        c5_cameraPR = get_enrich_order2( df, c5_all)
        dis_cameraPR = get_enrich_order2( df, disorders)
        dev_cameraPR = get_enrich_order2( df, genes_unique )
        print(sprintf('==== %s %s ====', r, m))
        print(adhd_cameraPR)
        print(dis_cameraPR)
        print(dev_cameraPR)
        ngood = sum(c5_cameraPR$PValue<.01)
        print(sprintf('Gene ontology nominal at p < .01: %d', ngood))
        out_fname = sprintf('~/data/methylation/res_c5_%s_%s.csv', r, m)
        write.csv(c5_cameraPR, file=out_fname)
    }
}
```

And this is what we get:

```
[1] "==== acc rsem ===="
        NGenes Direction    PValue       FDR
TWAS         8      Down 0.1296079 0.4245604
TWASnom     27      Down 0.1698242 0.4245604
CNV         41      Down 0.5079205 0.7891446
GWAS        16      Down 0.6556487 0.7891446
EWAS       110        Up 0.7891446 0.7891446
          NGenes Direction    PValue       FDR
ASD_GWAS      19      Down 0.1175165 0.6709865
ADHD_TWAS      8      Down 0.1296079 0.6709865
BD_GWAS       27      Down 0.1677466 0.6709865
BD_TWAS       10      Down 0.2720425 0.7273145
ASD_EWAS       4        Up 0.3035411 0.7273145
BD_EWAS       10      Down 0.3636573 0.7273145
SCZ_EWAS      41      Down 0.4434956 0.7602783
ADHD_GWAS     16      Down 0.6556487 0.9156570
SCZ_TWAS      36      Down 0.6893826 0.9156570
ADHD_EWAS    110        Up 0.7891446 0.9156570
ASD_TWAS       2        Up 0.8393523 0.9156570
SCZ_GWAS      13        Up 0.9859010 0.9859010
          NGenes Direction      PValue        FDR
overlap      857      Down 0.006134329 0.03680597
dev1_c0.9    347      Down 0.082706725 0.24812018
dev2_c0.9     60        Up 0.171652756 0.27127669
dev5_c0.9     68      Down 0.180851124 0.27127669
dev3_c0.9     69      Down 0.443558232 0.53226988
dev4_c0.9     21        Up 0.706900446 0.70690045
[1] "Gene ontology nominal at p < .01: 57"

[1] "==== acc kallisto ===="
        NGenes Direction    PValue       FDR
TWASnom     27      Down 0.1518427 0.7592135
CNV         41        Up 0.6038706 0.9753871
EWAS       110      Down 0.6489730 0.9753871
GWAS        16      Down 0.8181655 0.9753871
TWAS         8      Down 0.9753871 0.9753871
          NGenes Direction     PValue       FDR
BD_GWAS       27      Down 0.06195506 0.4782511
ASD_EWAS       5        Up 0.10741069 0.4782511
BD_EWAS       10      Down 0.11956276 0.4782511
ASD_TWAS       2        Up 0.54985346 0.9753871
SCZ_TWAS      36      Down 0.58560950 0.9753871
SCZ_EWAS      41      Down 0.58998193 0.9753871
ADHD_EWAS    110      Down 0.64897295 0.9753871
ASD_GWAS      20      Down 0.73848621 0.9753871
ADHD_GWAS     16      Down 0.81816549 0.9753871
BD_TWAS       10        Up 0.82303827 0.9753871
SCZ_GWAS      13      Down 0.92121583 0.9753871
ADHD_TWAS      8      Down 0.97538710 0.9753871
          NGenes Direction     PValue       FDR
overlap      861      Down 0.03588689 0.2153213
dev5_c0.9     68      Down 0.12198165 0.3659450
dev1_c0.9    349      Down 0.19175098 0.3835020
dev4_c0.9     22      Down 0.63489443 0.7694352
dev2_c0.9     60      Down 0.64119604 0.7694352
dev3_c0.9     69        Up 0.92294910 0.9229491
[1] "Gene ontology nominal at p < .01: 34"

[1] "==== caudate rsem ===="
        NGenes Direction    PValue       FDR
TWAS         8        Up 0.2004646 0.8009421
TWASnom     27        Up 0.4321536 0.8009421
GWAS        16      Down 0.5750725 0.8009421
EWAS       108      Down 0.6612845 0.8009421
CNV         41      Down 0.8009421 0.8009421
          NGenes Direction    PValue       FDR
BD_GWAS       28        Up 0.1870594 0.5460161
ADHD_TWAS      8        Up 0.2004646 0.5460161
BD_TWAS       10        Up 0.2023051 0.5460161
SCZ_GWAS      13      Down 0.2891505 0.5460161
ASD_GWAS      19      Down 0.3242903 0.5460161
SCZ_EWAS      39        Up 0.3366407 0.5460161
BD_EWAS        9      Down 0.3628382 0.5460161
ASD_TWAS       2      Down 0.3640107 0.5460161
SCZ_TWAS      36        Up 0.4811169 0.6273518
ASD_EWAS       4      Down 0.5422754 0.6273518
ADHD_GWAS     16      Down 0.5750725 0.6273518
ADHD_EWAS    108      Down 0.6612845 0.6612845
          NGenes Direction    PValue       FDR
dev4_c0.9     21      Down 0.3307942 0.9371264
overlap      858      Down 0.3966144 0.9371264
dev3_c0.9     69      Down 0.5596260 0.9371264
dev2_c0.9     58      Down 0.7260127 0.9371264
dev5_c0.9     67        Up 0.8828513 0.9371264
dev1_c0.9    348        Up 0.9371264 0.9371264
[1] "Gene ontology nominal at p < .01: 43"

[1] "==== caudate kallisto ===="
        NGenes Direction    PValue       FDR
EWAS       110        Up 0.2369190 0.8288190
TWAS         8        Up 0.3315276 0.8288190
TWASnom     27      Down 0.5426501 0.8910998
CNV         41      Down 0.7128798 0.8910998
GWAS        16      Down 0.9245132 0.9245132
          NGenes Direction     PValue       FDR
BD_TWAS       10        Up 0.03041902 0.3650282
ADHD_EWAS    110        Up 0.23691903 0.8728471
SCZ_GWAS      13      Down 0.24654243 0.8728471
ADHD_TWAS      8        Up 0.33152760 0.8728471
ASD_GWAS      20      Down 0.42550792 0.8728471
ASD_EWAS       5      Down 0.43642355 0.8728471
BD_GWAS       27        Up 0.63348550 0.9380155
SCZ_EWAS      41        Up 0.71784637 0.9380155
SCZ_TWAS      36      Down 0.72562827 0.9380155
BD_EWAS       10        Up 0.79048603 0.9380155
ADHD_GWAS     16      Down 0.92451317 0.9380155
ASD_TWAS       2      Down 0.93801548 0.9380155
          NGenes Direction    PValue       FDR
dev3_c0.9     69        Up 0.1683400 0.7842696
overlap      861      Down 0.4732712 0.7842696
dev2_c0.9     60        Up 0.4950908 0.7842696
dev5_c0.9     68      Down 0.5228464 0.7842696
dev4_c0.9     22      Down 0.8989833 0.9784994
dev1_c0.9    349      Down 0.9784994 0.9784994
[1] "Gene ontology nominal at p < .01: 43"
```

# 2020-10-02 09:44:17

Note that I changed the folder name from methylation to isoforms!

Let's give it a shot using fgsea, like:

http://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html

Mostly because it can use a pre-ranked list of genes. The question then is
mostly whether the rank can be just the p-values. I'll have to read the paper
later. But let's see if it works for now.

```r
library(fgsea)
load('~/data/rnaseq_derek/c5_gene_sets.RData')
for (r in c('acc', 'caudate')) {
    for (m in c('rsem', 'kallisto')) {
        fname = sprintf('~/data/isoforms/%s_%s.csv', r, m)
        df = read.csv(fname)
        colnames(df)[1] = 'hgnc_symbol'
        ranks = 1-df$pval
        names(ranks) = df$hgnc_symbol
        ranks = ranks[ranks > 0]
        res = fgsea(c5_all, ranks, scoreType='pos')
        ngood = sum(res$pval<.01)
        print(sprintf('Gene ontology nominal at p < .01: %d', ngood))
        out_fname = sprintf('~/data/isoforms/fgsea_c5_%s_%s.csv', r, m)
        # not writing out last column which is a list of genes and it's messing up write.csv... we can check it later if needed
        write.csv(res[order(pval), 1:(ncol(res)-1)], file=out_fname)
    }
}
```

To use goana for example, I'll need to Entrez GeneID.

```r
library(limma)
library('org.Hs.eg.db')
for (r in c('acc', 'caudate')) {
    for (m in c('rsem', 'kallisto')) {
        fname = sprintf('~/data/isoforms/%s_%s.csv', r, m)
        df = read.csv(fname)
        colnames(df)[1] = 'hgnc_symbol'
        ranks = 1-df$pval
        names(ranks) = df$hgnc_symbol
        ranks = ranks[ranks > 0]
        symbols = names(ranks)[ranks > .95]
        geneids = mapIds(org.Hs.eg.db, symbols, 'ENTREZID', 'SYMBOL')
        res_go = goana(geneids, species='Hs')
        res_go$FDR = p.adjust(res_go$P.DE, method='fdr')
        res_kegg = kegga(geneids, species='Hs')
        res_kegg$FDR = p.adjust(res_kegg$P.DE, method='fdr')
        out_fname = sprintf('~/data/isoforms/goana_p05_%s_%s.csv', r, m)
        write.csv(res_go[order(res_go$P.DE),], file=out_fname)
        out_fname = sprintf('~/data/isoforms/keggana_p05_%s_%s.csv', r, m)
        write.csv(res_kegg[order(res_kegg$P.DE),], file=out_fname)
    }
}
```

# TODO
*  how about something like this, thresholding in p-values? http://web.mit.edu/~r/current/arch/i386_linux26/lib/R/library/limma/html/goana.html
*  could also try fry/roast, but there might be better test-statistics that
   LAura's group could give us