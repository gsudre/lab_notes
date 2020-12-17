# 2020-12-16 20:05:18

IF there's a chance we won't be using real brain data, we might as well use
S-PrediXcan and just get the ffect of predicted expression in different brain
regions, using the GWAS results. Let' give that a try.

The first step is harmonization of the GWAS:

https://github.com/hakyimlab/MetaXcan/wiki/Tutorial:-GTEx-v8-MASH-models-integration-with-a-Coronary-Artery-Disease-GWAS

Maybe the idea could be: I have a list of genes expressed in ACC and then one
for Caudate, and their significance for ADHD as they come from the GWAS
imputation. Then we have a list of genes and their significance in regressing
against our brain phenotypes. Hopefully there is some overlap there. That would
mean that those genes are (predicted to be) differentially expressed in ADHD
ACC/caudate, and they are significantly expressed in our brain phenotype.
Ideally they would also be related to ADHD in our cohort, but that might be
hard. Maybe in ABCD?

I could also use these sPrediXcan results for the PM results, checking for
overlaps with that list. That's somewhat similar to the PRS regression results
though, except that it's region specific and does not use a nominal p-value in
the GWAS (at least not for the imputation part), while PRS has different
thresholds.

I could also potentially just carve out a region of the brain that is associated
with DX (voxelwise), hopefully falling within ACC, and go backwards from that.
That would only work for the ACC, but that's where our main results are anyways.

