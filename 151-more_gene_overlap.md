# 2020-11-20 15:04:46

Let's play a bit more with GeneOverlap just using our rnaseq results and our
ready-made gene sets:

```r
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
library(GeneOverlap)
library(WebGestaltR)
db_file = '~/data/post_mortem/acc_developmental.gmt'
gmt = readGmt(db_file) # already in gene symbols
# and convert it to lists
mylist = list()
for (s in unique(gmt$geneSet)) {
    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
}
t = .001
ra = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
rc = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
gom.obj <- newGOM(list(rnaseq_acc=ra, rnaseq_caudate=rc),
                  mylist, spec='hg19.gene')
getMatrix(gom.obj, name='pval')
```

Have to play a bit more with the thresholds, and potentially add lists with
overlapping genes in them?

