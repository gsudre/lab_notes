# 2020-12-23 09:56:06

## PM_enrichment
This is our main result, running the model `geneExpression ~ DX + nuiscancePCs`
in the post mortem (PM) data. It gives a list of genes with associated p-value
and statistics, and I ran that for ACC and Caudate separately. In most papers
they do a FDR cut-off in that list and run with the selected genes. Nothing
survives if we do that, so we explore these results in different ways.
 * ORA (over-representation analysis): given two lists of genes A and B,
   selected from a universe U of genes, calculate whether the overlap between A
   and B is significant.
 * GSEA (gene set enrichment analysis): given a **ranked** list of genes A and a
   gene set B, gives the probability that the genes in B are ranked higher than
   genes not in B. That's a simplistic explanation, and the different algorithms
   (WebGestalt (WG), camera, GAGE, fgsea) perform different math to calculate
   that probability. Another variable here is how to calculate the rank for each
   gene in A. For these results, I used `-log(P)*sign(t)` as the rank, and WG to
   calculate GSEA.

So, in this folder you'll see GSEA results for pmACC and pmCaudate, using
GeneOntology sets and also gene sets we devised ourselves.

### go_pics
I then selected 2 gene sets from Gene Ontology and made boxplots of the leading
genes in those sets. I made them for the raw `geneExpression ~ DX` but also for
`resid(geneExpression ~ nuiscancePCs) ~ DX`, but the results are similar. I
think it's easier to visualize them in the raw though.  As usual, in these
boxplots the lower and upper hinges correspond to the first and third quartiles
(the 25th and 75th percentiles). The whiskers extend from the hinge to the
largest/smallest value no further than 1.5*IQR from the hinge (where IQR is the
inter-quartile range, or distance 
between the first and third quartiles). Data beyond the end of the whiskers are
called "outlying" points and are plotted individually. The notches extend
`1.58*IQR / sqrt(n)`. This gives a roughly 95% confidence interval for comparing
medians. If notches don't overlap, it suggests that the medians are
significantly different.

## PM_overlap
I used different nominal p-value thresholds to select individual genes from
pmACC and pmCaudate results, and evaluated whether the ACC and Caudate overlap
was significant using ORA. I asked Gauri to check the 10 genes in the top result
(.png), and she compiled the list in the .xlsx.

## PRS
I used the genotypes in the PM dataset and the ADHD GWAS to calculate PRS. The
first result is PRS_models.txt, which only looks at `DX ~ PRS` to establish
whether PRS is related to ADHD at all in this cohort, regardless of gene
expression. We normally use that to select a PRS threshold based on R^2. In this
case, we'd select PRS0.5.

We can then replace DX by PRS to run `geneExpression ~ PRS + nuiscancePCs` in
the PM cohort, again split between ACC and caudate. That also results in a list
of genes with p-values and statistics. I thresholded the lists using different
nominal p-values, and did an ORA analysis against the genes in the
`geneExpression ~ DX + nuiscancePCs` result, and that's listed in the two .csv
files.

However, that result doesn't take into consideration whether genes are over or
under expressed. I then created *UpDown*.csv results, which calculates the
overlap only within the list of genes up or down regulated, based on their
tstats. The issue then is how to calculate the gene universe for the
hypergeometric test. We select a certain number of PRsgenes and PMgenes, and
check if their overlap is significant. The bigger the universe of possible
selections, the more significant the overlap. When not splitting between
up/down, the universe is straight-forward: it's simply all genes we're working
with. But when we split it, it could be argued that the universe is only the
genes up (or down), or it's still the entire gene universe (as being up or down
could be seen as another form of selection). I provided p-values for both
(pvalWhole is the latter scenario).

Even though the up/down splits weren't veyr significant, we can still visualize
them by splitting the more significant results that used just the pvalues into
up and down-regulated, and displaing them in Venn diagrams (venn_pics.zip).

## ABCD_overlap
Here we use the ABCD genotypes to impute the gene expression in the brain. That
gives two big tables (ACC and Caudate), listing the imputed gene expression for
each subject. We then use the brain scans of each subject in the model
`impGeneExpression ~ brain + age + sex + fsqc_qu_motion + fsqc_qu_pialover + fsqc_qu_wmunder + fsqc_qu_inhomogeneity + fsqc_qu_artifact`. As usual, that
generates a list of genes with associated stats, which we can compare to the PM
gene lists.

Similarly to what I did for the PRS analysis, I ran ORA for different nominal
p-value cut-offs (imp*.pnc). I also listed the 
Another way to check the overlap between PRS and DX result is to calculate
Spearman correlation between the two ranked list of genes. Here, we can rank the
lists based on p-value, or the same rank used for GSEA (`-log(P) * sign(t)`).
# TODO
[8:33 AM] Shaw, Philip (NIH/NHGRI) [E]
    1) GWAS (xPrediXcan) results for gene expression in ACC/caudate--- the correlation (p and logP*signFC); the hyper.   Just to know there ins't much there.  
​[8:34 AM] Shaw, Philip (NIH/NHGRI) [E]
    2) For the ABCD--- the hypergeo, when split by up/down regulation.  (Just to see how far away from sig this is). 


[Yesterday 9:25 PM] Sudre, Gustavo (NIH/NHGRI) [E]
    FYI, here are the items I still have in my list:

	make Venn diagrams for PRS PM overlaps
	run DX*brain for PM brain
	beef up WGCNA results (check ACC clusters robustness, make pictures of the genes there, repeat everything for Caudate)
	check PRS to ABCD imputation overlap (regardless of circularity)
	check how diagnosis for each brain bank was conducted
	run more traditional DTE and DTU pipelines while paper is being written 

