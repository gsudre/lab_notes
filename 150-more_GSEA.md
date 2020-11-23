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

# 2020-11-23 06:45:50

Let's first use the gene description so it's easier to tell the results apart:

```r
library(WebGestaltR)
library(fgsea)

load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
ranks = rnaseq_acc$logFC
names(ranks) = rnaseq_acc$hgnc_symbol

db = 'geneontology_Molecular_Function_noRedundant'
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
fgseaRes <- fgsea(pathways=mylist, stats=ranks,
                  scoreType='std')
m = merge(fgseaRes, gs$geneSetDes, by.x='pathway', by.y='geneSet')
m = m[order(m$pval), ]
```

So now we can appreciate the differences between running logFC and the pvalue a
bit better:

```r
library(WebGestaltR)
library(fgsea)

load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
ranks = rnaseq_acc$logFC
names(ranks) = rnaseq_acc$hgnc_symbol

db = 'geneontology_Biological_Process_noRedundant'
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
fgseaRes <- fgsea(pathways=mylist, stats=ranks, scoreType='std')
m = merge(fgseaRes, gs$geneSetDes, by.x='pathway', by.y='geneSet')
m = m[order(m$pval), ]

cp <- collapsePathways(m[m$padj < 0.05, ], mylist, ranks)
mp <- m[pathway %in% cp$mainPathways, ]
```

This is my list even after collapsing (from 70 to 45):

```
 1:                                  neutrophil mediated immunity
 2:                                                  phagocytosis
 3:                      response to molecule of bacterial origin
 4:                                          synapse organization
 5:                                      adaptive immune response
 6:                                   acute inflammatory response
 7:                  immune response-regulating signaling pathway
 8:                          regulation of innate immune response
 9:                          cellular response to biotic stimulus
10:                        regulation of trans-synaptic signaling
11:                                         forebrain development
12:                                                  angiogenesis
13:                            defense response to other organism
14:                         regulation of vasculature development
15:                                      translational initiation
16:                 central nervous system neuron differentiation
17:                                              axon development
18:                                                     cognition
19:                                           leukocyte migration
20:                           regulation of inflammatory response
21:                                        platelet degranulation
22:                            regulation of leukocyte activation
23:                                      response to ammonium ion
24:          positive regulation of response to external stimulus
25:                                      interleukin-2 production
26:                                       NIK/NF-kappaB signaling
27:                                      interleukin-8 production
28: vascular endothelial growth factor receptor signaling pathway
29:                          serotonin receptor signaling pathway
30:                   multicellular organismal response to stress
31:                                          mast cell activation
32:                                         neuromuscular process
33:                                                   coagulation
34:                         activation of protein kinase activity
35:                               regulation of body fluid levels
36:                               response to mechanical stimulus
37:                                natural killer cell activation
38:                                     response to oxygen levels
39:                          multicellular organismal homeostasis
40:                          extracellular structure organization
41:                 protein localization to endoplasmic reticulum
42:                                     response to interleukin-1
43:                         peripheral nervous system development
44:                                         RNA catabolic process
45:                          multi-multicellular organism process
```

Then, if using `ranks = -log(rnaseq_acc$P.Value) * sign(rnaseq_acc$logFC)` we
get:

```
                                                      description
 1:                                          synapse organization
 2:                        regulation of trans-synaptic signaling
 3:                   regulation of neuron projection development
 4:                                              axon development
 5:                                                     cognition
 6:             negative regulation of nervous system development
 7:                                      response to ammonium ion
 8:           negative regulation of cell projection organization
 9:                              regulation of membrane potential
10:                 central nervous system neuron differentiation
11:                                      translational initiation
12:                     regulation of ion transmembrane transport
13:                    microtubule organizing center organization
14:     cell-cell adhesion via plasma-membrane adhesion molecules
15:                          serotonin receptor signaling pathway
16:                                         forebrain development
17:                              membrane lipid metabolic process
18:                  immune response-regulating signaling pathway
19:                   multicellular organismal response to stress
20:                                    double-strand break repair
21:                                                 glycosylation
22:                                         RNA catabolic process
23:                            defense response to other organism
24:                                          spindle organization
25:                                         neuromuscular process
26:                         regulation of neurotransmitter levels
27:                  chemical synaptic transmission, postsynaptic
28:                                                signal release
29:                                          exploration behavior
30:                                              feeding behavior
31:                           mitotic cell cycle phase transition
32:                                                adult behavior
33:                                         cell cycle checkpoint
34:                                      adaptive immune response
35:                 protein localization to endoplasmic reticulum
36:                                        chromosome segregation
37:                                                   gliogenesis
38: vascular endothelial growth factor receptor signaling pathway
```

This seems to be a more curated list?

Also, there seems to be a clear difference between looking at positive only,
negative only, and both sides:

```
r$> head(fgseaResP[order(pval),])                                                                                                
      pathway         pval        padj   log2err        ES      NES size                               leadingEdge
1: GO:0006413 5.910473e-06 0.005023902 0.6105269 0.4252248 2.108822  175   EIF1,RPS29,RPS12,METTL3,POLR2G,RPL3,...
2: GO:0031023 2.499296e-05 0.010622006 0.5756103 0.4559576 2.217141  120 CEP152,XRCC2,CEP135,CENPJ,HAUS3,CNTLN,...
3: GO:0002764 6.260800e-05 0.017738933 0.5384341 0.3284636 1.675616  386   TLR5,FCGR3A,ALPK1,FCGR2A,SOS1,RUNX1,...
4: GO:0006302 2.074548e-04 0.037194479 0.5188481 0.3692989 1.843308  198    XRCC2,RAD52,PARP2,CHEK2,UIMC1,EME1,...
5: GO:0007051 2.187911e-04 0.037194479 0.5188481 0.3893469 1.923952  152  HAUS3,CHEK2,CENPE,CEP192,POC1A,AURKA,...
6: GO:0006310 2.756957e-04 0.037523494 0.4984931 0.3583753 1.799487  226   XRCC2,RAD52,TNFSF4,THOC1,EME1,RBBP8,...

r$> head(fgseaResN[order(pval),])                                                                                                
      pathway         pval         padj   log2err         ES       NES size                                 leadingEdge
1: GO:0050803 1.000000e-10 2.833333e-08        NA -0.4962117 -2.559987  190   AMIGO2,LRRC4,PCDH8,LRRTM1,RAB17,LRRN1,...
2: GO:0050808 1.000000e-10 2.833333e-08        NA -0.4518863 -2.401316  337   AMIGO2,LRRC4,EFNB2,PCDH8,LRRTM1,RAB17,...
3: GO:0099177 1.000000e-10 2.833333e-08        NA -0.4400671 -2.351769  379 LRRC4,CALB1,GUCY1B1,BAIAP3,LRRTM1,PTK2B,...
4: GO:0097485 1.104329e-08 2.346698e-06 0.7477397 -0.4375328 -2.273277  217       EFNB2,GFRA2,CHL1,CDH4,FEZ2,VSTM2L,...
5: GO:0010975 1.518659e-08 2.581720e-06 0.7337620 -0.3654818 -1.960130  437       EFNB2,NDEL1,CDH4,ID1,RAB17,TRIM46,...
6: GO:0061564 2.608957e-08 3.696022e-06 0.7337620 -0.3630050 -1.945744  419      EFNB2,TRIM32,GFRA2,NDEL1,CHL1,CDH4,...

r$> head(fgseaResS[order(pval),])                                                                                                
      pathway         pval         padj   log2err         ES       NES size                                 leadingEdge
1: GO:0050808 1.000000e-10 4.250000e-08        NA -0.4518863 -2.075079  337   AMIGO2,LRRC4,EFNB2,PCDH8,LRRTM1,RAB17,...
2: GO:0099177 1.000000e-10 4.250000e-08        NA -0.4400671 -2.049806  379 LRRC4,CALB1,GUCY1B1,BAIAP3,LRRTM1,PTK2B,...
3: GO:0050803 2.103509e-10 5.959943e-08 0.8266573 -0.4962117 -2.138185  190   AMIGO2,LRRC4,PCDH8,LRRTM1,RAB17,LRRN1,...
4: GO:0010975 2.686641e-08 5.709112e-06 0.7337620 -0.3654818 -1.723070  437       EFNB2,NDEL1,CDH4,ID1,RAB17,TRIM46,...
5: GO:0097485 7.383692e-08 1.255228e-05 0.7049757 -0.4375328 -1.909181  217       EFNB2,GFRA2,CHL1,CDH4,FEZ2,VSTM2L,...
6: GO:0061564 9.967538e-08 1.321325e-05 0.7049757 -0.3630050 -1.701787  419      EFNB2,TRIM32,GFRA2,NDEL1,CHL1,CDH4,...
```

Which leads me to believe we need to do all of them. This will become especially
important when we compare it to GAGE results. Before we do that, let's clean up
the duplicated genes in the results

# TODO
 * What's the effect of setting the score type to pos or neg?
 * what's the effect of just using the positive ranks or the negative ranks?
 * deal with duplicate ids
 * play with EPS for fgsea
 * gage
 * try overlap with gene lists too (our lists)