# 2020-12-29 07:44:26

Let's try to play a bit with FUSION to see if their TWAS results are a bit
different than S-PredXcan:

http://gusevlab.org/projects/fusion/

I first run their tool to clean up the file, and it spits out Z score byt itself:

```bash
# bw
cd ~/data/expression_impute/fusion_twas-master/
cd ldsc-master
module load python/2.7
./munge_sumstats.py --sumstats ~/pgc2017/adhd_eur_jun2017 --N 55374 \
    --N-cas 20183 --N-con 35191 --signed-sumstats OR,1 \
    --out ../adhd_eur_jun2017
cd ..
module load R
for c in {1..22}; do
    echo $c;
    Rscript FUSION.assoc_test.R \
        --sumstats adhd_eur_jun2017.sumstats.gz \
        --weights ./WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.P01.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr $c \
        --out ADHD_ACC.${c}.dat;
done
```

And do the same thing for the Caudate:

```bash
for c in {1..22}; do
    echo $c;
    Rscript FUSION.assoc_test.R \
        --sumstats adhd_eur_jun2017.sumstats.gz \
        --weights ./WEIGHTS/Brain_Caudate_basal_ganglia.P01.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr $c \
        --out ADHD_Caudate.${c}.dat;
done
```

Now it's just a matter of combining all results and comparing them to what we
see in the PM data in R:

```r
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_ACC.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)

library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
both_res = merge(rnaseq_acc, fusion, by.x='hgnc_symbol', by.y='ID',
                 all.x=F, all.y=F)
thresh = c(.05, .01, .005, .001)
imp_nums = vector(length=length(thresh), mode='numeric')
rna_nums = vector(length=length(thresh), mode='numeric')
pvals = matrix(data=NA, nrow=length(thresh), ncol=length(thresh))
rownames(pvals) = sapply(thresh, function(x) sprintf('PM_%.3f', x))
colnames(pvals) = sapply(thresh, function(x) sprintf('IMP_%.3f', x))
inter_nums = pvals
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$TWAS.P < thresh[ti], 'hgnc_symbol']
    imp_nums[ti] = length(imp_genes)
    for (tr in 1:length(thresh)) {
        rna_genes = both_res[both_res$P.Value < thresh[tr], 'hgnc_symbol']
        rna_nums[tr] = length(rna_genes)
        go.obj <- newGeneOverlap(imp_genes, rna_genes,
                                 genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(imp_genes, rna_genes)
        pval = getPval(go.obj)
        pvals[tr, ti] = pval
        inter_nums[tr, ti] = length(intersect(rna_genes, imp_genes))
    }
}
```

Close...

```
r$> pvals                                                                                               
          IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050 0.07757113 0.4869580 0.8153215 0.2507528
PM_0.010 0.37652053 0.7077565 0.5505896 0.2265755
PM_0.005 0.41763284 0.5076920 0.3691940 0.1375818
PM_0.001 1.00000000 1.0000000 1.0000000 1.0000000

r$> inter_nums                                                                                          
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050        21         5         2         2
PM_0.010         5         1         1         1
PM_0.005         3         1         1         1
PM_0.001         0         0         0         0
```

Let's at least check if it's the same for the caudate:

```r
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_Caudate.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
both_res = merge(rnaseq_caudate, fusion, by.x='hgnc_symbol', by.y='ID',
                 all.x=F, all.y=F)
thresh = c(.05, .01, .005, .001)
imp_nums = vector(length=length(thresh), mode='numeric')
rna_nums = vector(length=length(thresh), mode='numeric')
pvals = matrix(data=NA, nrow=length(thresh), ncol=length(thresh))
rownames(pvals) = sapply(thresh, function(x) sprintf('PM_%.3f', x))
colnames(pvals) = sapply(thresh, function(x) sprintf('IMP_%.3f', x))
inter_nums = pvals
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$TWAS.P < thresh[ti], 'hgnc_symbol']
    imp_nums[ti] = length(imp_genes)
    for (tr in 1:length(thresh)) {
        rna_genes = both_res[both_res$P.Value < thresh[tr], 'hgnc_symbol']
        rna_nums[tr] = length(rna_genes)
        go.obj <- newGeneOverlap(imp_genes, rna_genes,
                                 genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(imp_genes, rna_genes)
        pval = getPval(go.obj)
        pvals[tr, ti] = pval
        inter_nums[tr, ti] = length(intersect(rna_genes, imp_genes))
    }
}
```

As expected, not as strong as ACC:

```
r$> pvals                                                                                               
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050 0.7408160 0.8359501 0.9806482 0.7155434
PM_0.010 0.2731887 1.0000000 1.0000000 1.0000000
PM_0.005 0.6344386 1.0000000 1.0000000 1.0000000
PM_0.001 1.0000000 1.0000000 1.0000000 1.0000000

r$> inter_nums                                                                                          
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050        13         4         1         1
PM_0.010         5         0         0         0
PM_0.005         2         0         0         0
PM_0.001         0         0         0         0
```

Let's take a look at the correlation analysis, first the absolute but also the
ranked version:

```
r$> cor.test(both_res$P.Value, both_res$TWAS.P, method='spearman')                                      

        Spearman's rank correlation rho

data:  both_res$P.Value and both_res$TWAS.P
S = 2047665896, p-value = 0.4026
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
0.01738134 

Warning message:
In cor.test.default(both_res$P.Value, both_res$TWAS.P, method = "spearman") :
  Cannot compute exact p-value with ties

r$> cor.test(sign(both_res$TWAS.Z)*both_res$TWAS.P, sign(both_res$t)*both_res$P.Value, method='spearman'
    )                                                                                                   

        Spearman's rank correlation rho

data:  sign(both_res$TWAS.Z) * both_res$TWAS.P and sign(both_res$t) * both_res$P.Value
S = 1994053265, p-value = 0.03783
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
0.04310857 

Warning message:
In cor.test.default(sign(both_res$TWAS.Z) * both_res$TWAS.P, sign(both_res$t) *  :
  Cannot compute exact p-value with ties
```

Interestingly, the ranked correlation for Caudate worked... how does it look for
ACC?

```
r$> cor.test(both_res$P.Value, both_res$TWAS.P, method='spearman')                                      

        Spearman's rank correlation rho

data:  both_res$P.Value and both_res$TWAS.P
S = 756966488, p-value = 0.7081
alternative hypothesis: true rho is not equal to 0
sample estimates:
         rho 
-0.009222311 

Warning message:
In cor.test.default(both_res$P.Value, both_res$TWAS.P, method = "spearman") :
  Cannot compute exact p-value with ties

r$> cor.test(sign(both_res$TWAS.Z)*both_res$TWAS.P, sign(both_res$t)*both_res$P.Value, method='spearman'
    )                                                                                                   

        Spearman's rank correlation rho

data:  sign(both_res$TWAS.Z) * both_res$TWAS.P and sign(both_res$t) * both_res$P.Value
S = 738158849, p-value = 0.5198
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
0.01585289 
```

No...

## Up/Down regulation

Even though I don't think much is going to come out of this, it's worth trying
the up/down analysis too:

```r
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$TWAS.P < thresh[ti] & both_res$TWAS.Z > 0,
                         'hgnc_symbol']
    imp_nums[ti] = length(imp_genes)
    for (tr in 1:length(thresh)) {
        rna_genes = both_res[both_res$P.Value < thresh[tr] & both_res$t > 0,
                             'hgnc_symbol']
        rna_nums[tr] = length(rna_genes)
        go.obj <- newGeneOverlap(imp_genes, rna_genes,
                                 genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(imp_genes, rna_genes)
        pval = getPval(go.obj)
        pvals[tr, ti] = pval
        inter_nums[tr, ti] = length(intersect(rna_genes, imp_genes))
    }
}
```

```
r$> inter_nums                                                                                          
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050        10         4         2         2
PM_0.010         4         1         1         1
PM_0.005         3         1         1         1
PM_0.001         0         0         0         0

r$> pvals                                                                                               
          IMP_0.050  IMP_0.010 IMP_0.005  IMP_0.001
PM_0.050 0.04412797 0.05261103 0.2439800 0.05370831
PM_0.010 0.07377422 0.35080305 0.2454928 0.10919924
PM_0.005 0.10083557 0.26159542 0.1794159 0.07796568
PM_0.001 1.00000000 1.00000000 1.0000000 1.00000000
```

Maybe there's a hint here for upregulated, but nothing for down (not shown).

# 2021-01-06 20:40:45

Philip asked me to look a bit closer at these results, maybe going through with
the analysis further. Before I do that, there are a few weird things. For
example:

```
r$> dim(rnaseq_acc)                                                                                             
[1] 17677     9

r$> dim(fusion)                                                                                                 
[1] 2603    3

r$> dim(both_res)                                                                                               
[1] 1651   11
```

So, I lose almost half of my fusion genes when I merge them with the PM results.
That's a bit weird. What's going on there?

```
r$> sum(! fusion$ID %in% G_list$hgnc_symbol)                                                                    
[1] 768
```

So, close to 800 genes in fusion are not even in our initial lsit of 38K genes.
In other words, it's not our cleaning that's doing this. 

Would we get a different result if we used ALL the genes, instead of just the
p<.01 most heritable ones?

http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Anterior_cingulate_cortex_BA24.ALL.tar.bz2
http://gusevlab.org/projects/fusion/weights/GTEx.Brain_Caudate_basal_ganglia.ALL.tar.bz2

# 2021-01-07 07:23:02

To be precise, in the weights file for ACC there are 2109 genes, but 804 of
those (30%) are not listed in our universe of PM genes. And that's even before
any sort of cleaning. What's going on here? Are they even in Derek's initial
file?

As a clue, it looks like gtex v7 used hg19:

https://gtexportal.org/home/datasets

They list the gtf files they used, so maybe that's a good thing to start with?

I think I might have better luck if I replace the IDs in the results file by
their Gene ID (ENS).

```r
wei = read.delim('~/data/expression_impute/fusion_twas-master/WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.P01.pos')
cnt = 75
wei$geneid = substring(wei$WGT, cnt, cnt+14)
```

Now we're down to 524, but that's still 20% of the predicted genes. Is it
because they were removed for not having hugo ids? Yes, if I just look at the
headers in the data form Derek, I'd be missing only 111 (5%). 

The issue here is that it looks like Derek aligned his counts to GrCh38, and
GTex v7 is defined on hg19. So, a few variants don't exist anymore, or have
switched gene ids. BTW, Gtex v8 is in hg38. Here's an example:

```
r$> wei[180,]                                                                                   
                                   PANEL
180 Brain_Anterior_cingulate_cortex_BA24
                                                                                                     WGT
180 Brain_Anterior_cingulate_cortex_BA24/Brain_Anterior_cingulate_cortex_BA24.ENSG00000090920.9.wgt.RDat
       ID CHR       P0       P1   N          geneid
180 FCGBP  19 40353963 40440533 109 ENSG00000090920

r$> which(G_list$hgnc_symbol=='FCGBP')                                                          
[1] 37023

r$> G_list[37023,]                                                                              
      ensembl_gene_id hgnc_symbol chromosome_name
53281 ENSG00000275395       FCGBP              19

r$> which(id_num=='ENSG00000275395')                                                            
ENSG00000275395.6 
            53392 
```

In other words, in the weight file (gtex v7), FCGBP is mapped to
ENSG00000090920, which is confirmed by the ensembl website:

http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000090920;r=19:40353963-40440533

However, if we look for that same ID in hg38 it says it's been deprecated. If we
look for FCGBP in hg38, we get ENSG00000275395, which is in Derek's header, and
also in our G_list. This complicates things a bit, and I'll need to make sure
the rest of the analysis was all done with the same reference.

For reference, the GWAS was imputed to 1000 Genomes Project Phase 3 reference
panel, which is in hg19/GRCh37. So, at least the GWAS going to fusion is in the same
reference as the fusion database.

I'm somewhat curious to see what would happen to the DGE analysis if it was all
done in hg19 as well... if anything, the PRS for the PM cohort was. 

Before I re-run this whole thing, let me revisit the FUSION results again... so,
it turns out that 4 hits in fusion acc survive bonferroni and 9 survive FDR. But
of those 9, only 4 are in the results matrix (and in the G_list), but all 9 are
actually in the matrix Derek sent... so, that would be another argument for
re-running the whole analysis without removing those genes a priori. 

I re-ran it and now I have 8 out of the 9 that survive FDR, so I might be able
to do some GSEA with that?


# 2021-01-08 07:12:50

So, it seems to me that we can try to really narrow down the FUSION results, and
then see if those genes have anything interesting in PM, or increase the pool of
FUSION genes. Let's try both. Also, another option is try to go through
SPreiXcan again.

```bash
# bw
cd ~/data/expression_impute/fusion_twas-master/
module load R
for c in {1..22}; do
    echo $c;
    Rscript FUSION.assoc_test.R \
        --sumstats adhd_eur_jun2017.sumstats.gz \
        --weights ./WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr $c \
        --out ADHD_ACC.ALL.${c}.dat;
done
```

And do the same thing for the Caudate:

```bash
for c in {1..22}; do
    echo $c;
    Rscript FUSION.assoc_test.R \
        --sumstats adhd_eur_jun2017.sumstats.gz \
        --weights ./WEIGHTS/Brain_Caudate_basal_ganglia.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr $c \
        --out ADHD_Caudate.ALL.${c}.dat;
done
```

I'm getting some warnings like this:

```
If a large number of genes were skipped, verify that your GWAS Z-scores, expression weights, and LDREF data use the same SNPs (or nearly)
Or consider pre-imputing your summary statistics to the LDREF markers using summary-imputation software such as [https://github.com/bogdanlab/fizi]
```

Here I'm seeing about 200 out of 1100 SNPs skipped, and it was much less for the
P01 set. But would this be something interesting to try regardless?

Let's see what's our gene overlap now, first ACC:

```r
# bw
mydir = '~/data/expression_impute/fusion_twas-master/'
fusion = c()
for (c in 1:22) {
    fname = sprintf('%s/ADHD_ACC.ALL.%d.dat', mydir, c)
    tmp = read.delim(fname)
    fusion = rbind(fusion, tmp[, c('ID', 'TWAS.Z', 'TWAS.P')])
}
fusion$TWAS.P = as.numeric(fusion$TWAS.P)
fusion = fusion[!is.na(fusion$TWAS.P), ]
fusion$TWAS.Z = as.numeric(fusion$TWAS.Z)

wei = read.delim('~/data/expression_impute/fusion_twas-master/WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.pos')
cnt = 75
wei$geneid = substring(wei$WGT, cnt, cnt+14)
mf = merge(fusion, wei, by='ID', all.x=T, all.y=F, sort=F)
mf$adjPval = p.adjust(mf$TWAS.P, method='fdr')
# we now have 11269 genes in the FUSION results for ACC.
# 32 are good for FDR and 5 for Bonferroni

load('~/tmp/ensid2.rdata')

library(GeneOverlap)
both_res = merge(rnaseq_acc, mf, by.x='GENEID', by.y='geneid',
                 all.x=F, all.y=F)
# of those, 9K are also in the rnaseq_acc results
thresh = c(.05, .01, .005, .001)
imp_nums = vector(length=length(thresh), mode='numeric')
rna_nums = vector(length=length(thresh), mode='numeric')
pvals = matrix(data=NA, nrow=length(thresh), ncol=length(thresh))
rownames(pvals) = sapply(thresh, function(x) sprintf('PM_%.3f', x))
colnames(pvals) = sapply(thresh, function(x) sprintf('IMP_%.3f', x))
inter_nums = pvals
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$TWAS.P < thresh[ti], 'GENEID']
    imp_nums[ti] = length(imp_genes)
    for (tr in 1:length(thresh)) {
        rna_genes = both_res[both_res$P.Value < thresh[tr], 'GENEID']
        rna_nums[tr] = length(rna_genes)
        go.obj <- newGeneOverlap(imp_genes, rna_genes,
                                 genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(imp_genes, rna_genes)
        pval = getPval(go.obj)
        pvals[tr, ti] = pval
        inter_nums[tr, ti] = length(intersect(rna_genes, imp_genes))
    }
}
print(pvals)
print(inter_nums)
```

```
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050 0.6588256 0.4303270 0.4039343  0.237694
PM_0.010 0.9603273 0.6207155 0.6392440  0.370203
PM_0.005 0.9261106 1.0000000 1.0000000  1.000000
PM_0.001 0.6570605 1.0000000 1.0000000  1.000000
         IMP_0.050 IMP_0.010 IMP_0.005 IMP_0.001
PM_0.050        64        22        14         7
PM_0.010        11         5         3         2
PM_0.005         5         0         0         0
PM_0.001         2         0         0         0
```

Overlap results are still no good...

## More parameters

While we run the ALL analysis above, let's also run the Assoc function with the
ability to generate the COLOR and permutation inputs.

```bash
# bw
cd ~/data/expression_impute/fusion_twas-master/
module load R
for c in {1..22}; do
    echo $c;
    Rscript FUSION.assoc_test.R \
        --sumstats adhd_eur_jun2017.sumstats.gz \
        --weights ./WEIGHTS/Brain_Anterior_cingulate_cortex_BA24.P01.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr $c \
        --perm 10000 --coloc_P .05 --GWASN 55374 \
        --out ADHD_ACC.P01_p05_10K_${c}.dat;
done
```

And do the same thing for the Caudate:

```bash
for c in {1..22}; do
    echo $c;
    Rscript FUSION.assoc_test.R \
        --sumstats adhd_eur_jun2017.sumstats.gz \
        --weights ./WEIGHTS/Brain_Caudate_basal_ganglia.P01.pos \
        --weights_dir ./WEIGHTS/ \
        --ref_ld_chr ./LDREF/1000G.EUR. \
        --chr $c \
        --perm 10000 --coloc_P .05 --GWASN 55374 \
        --out ADHD_Caudate.P01_p05_10K_${c}.dat;
done
```

# TODO
 * try fizi
 * try permutation test (see their FAQ "How can I validate the TWAS associations)
 * try coloc?
 * try focus?
 * try mesc?
 * SPrediXcan again?
 * improve FUSION results by checking these messages?
  ```WARNING :  Brain_Anterior_cingulate_cortex_BA24 Brain_Anterior_cingulate_cortex_BA24/Brain_Anterior_cingulate_cortex_BA24.ENSG00000186891.9.wgt.RDat TNFRSF18 1 1138888 1142071 109 had 144 / 231 non-overlapping GWAS Z-scores, skipping this gene.
WARNING :  Brain_Anterior_cingulate_cortex_BA24 Brain_Anterior_cingulate_cortex_BA24/Brain_Anterior_cingulate_cortex_BA24.ENSG00000272420.1.wgt.RDat RP4-740C4.7 1 2294500 2295067 109 had mean GWAS Z-score imputation r2 of 0.6148983 at expression weight SNPs, skipping this gene.
```