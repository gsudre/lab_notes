# 2020-11-20 12:23:42

I've been doing some digging, and we need to be careful with the gene sets we're
feeding in. Looking at these:

* https://bioinformaticsbreakdown.com/how-to-gsea/
* https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#gsea-analysis
* http://www.bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html
* http://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/RNA-seqWorkflow.pdf
* http://bioconductor.org/packages/release/bioc/vignettes/gage/inst/doc/gage.pdf

So, it looks like GO sets are usually in the same direction (e.g. all over or
under-expressed), but that's not true for KEGG pathways, and possibly not true
either for the gene sets we cooked up. So we need to be able to specify that in
our analysis. It looks like both fgsea and gage (the latter for sure) let us do
that. 

I also have a concern about camera because I'm not estimating the inter-gene
correlation, which is one of its main things. I'm using the default, and the
results somewhat match WG results, but I'm not sure that's enough. When I tried
to estimate the correlation using the rnaseq_acc data it wiped out the results,
which could just mean that the 50 or so samples we have are not enough to
estimate the numbers properly. So, it might be better to use some technique that
doesn't require that.

Another difference is that the recommendation is to use log2(fold change) or,
another option, -log10({p value}) * sign({Fold Change}). The sign of fold change
is the same as sign Tstat, if that's not available (i.e. imputation).

OK, so let's try running fgsea for the rnaseq results and see what we get:

```r
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
library(fgsea)
library(org.Hs.eg.db)
pathways <- gmtPathways('~/data/post_mortem/hsapiens_geneontology_Molecular_Function_entrezgene.gmt')
a = mapIdsList(org.Hs.eg.db, keys=rnaseq_acc$hgnc_symbol, column="ENTREZID", keytype="SYMBOL")
entrezids = unlist(a, use.names=F)
ranks = rnaseq_acc$logFC
names(ranks) = entrezids
fgseaRes <- fgsea(pathways=pathways, stats=ranks,
                  scoreType='std')
collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.05], 
                                      pathways, ranks)
mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
                         order(-NES), pathway]
plotGseaTable(pathways[mainPathways], ranks, fgseaRes, gseaParam = 0.5)
```

Do the results look similar using the other type of rank?

```r
ranks2 = -log(rnaseq_acc$P.Value) * sign(rnaseq_acc$logFC)
names(ranks2) = entrezids
fgseaRes2 <- fgsea(pathways=pathways, stats=ranks2,
                  scoreType='std')
collapsedPathways2 <- collapsePathways(fgseaRes2[order(pval)][padj < 0.05], 
                                      pathways, ranks)
mainPathways2 <- fgseaRes2[pathway %in% collapsedPathways2$mainPathways][
                         order(-NES), pathway]
plotGseaTable(pathways[mainPathways2], ranks2, fgseaRes2, gseaParam = 0.5)
```

The correlation between the two rank vectors is .77 (p < 2.2e-16), and .95 for
Spearman. 


# TODO
 * What's the effect of setting the score type to pos or neg?
 * add go set description so they are displayed in results
 * deal with duplicate ids
 * gage
 * try overlap with gene lists too (our lists)