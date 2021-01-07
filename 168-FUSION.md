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

# TODO
 * why are some predicted genes no in our study? issue with reference genome?
 * try ALL instead of p01?
 * try coloc?
 * try focus?
 * try mesc?