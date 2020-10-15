# 2020-10-15 06:52:24

Let's try different types of gene set analysis on the results of the 3 different
modalities: rnaseq, methylation, and isoforms. The latter need some extra
filtering based on Valer's observations, so I'll keep that for last.

I'll use this:
https://www.bioconductor.org/packages/devel/bioc/vignettes/GeneOverlap/inst/doc/GeneOverlap.pdf
mostly because it implements the hypergeometric test in a fast way.

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/xmodal_results_10152020.RData')
rnaseq_acc_p05 = rnaseq_acc[rnaseq_acc$P.Value < .05, 'hgnc_symbol']
rnaseq_caudate_p05 = rnaseq_caudate[rnaseq_caudate$P.Value < .05, 'hgnc_symbol']
gsizeu = length(union(rnaseq_acc$hgnc_symbol, rnaseq_caudate$hgnc_symbol))
gsizei = length(intersect(rnaseq_acc$hgnc_symbol, rnaseq_caudate$hgnc_symbol))
go.obj <- testGeneOverlap(newGeneOverlap(rnaseq_acc_p05, rnaseq_caudate_p05,
                          genome.size=gsizeu))
```

I'm using the union as the gene universe because otherwise we might have genes
in one list that are not in the universe. The universe for the different lists
will change because the genes we're testing in each are different. But we will
need to be careful that they're unique in each result variable!

```r
rnaseq_acc_p01 = rnaseq_acc[rnaseq_acc$P.Value < .01, 'hgnc_symbol']
rnaseq_caudate_p01 = rnaseq_caudate[rnaseq_caudate$P.Value < .01, 'hgnc_symbol']
gom.obj <- newGOM(list(acc_p01=rnaseq_acc_p01, acc_p05=rnaseq_acc_p05),
                  list(caudate_p01=rnaseq_caudate_p01, caudate_p05=rnaseq_caudate_p05),
                  genome.size=gsizeu)
drawHeatmap(gom.obj)
getMatrix(gom.obj, name='pval')
getMatrix(gom.obj, name='odds.ratio')
```

![](images/2020-10-15-07-31-13.png)

# TODO
 * maybe try using spec=hg19 to check on differences?