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

# 2020-12-17 07:03:14

Now that we have the sample data, we run:

```bash
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate imlabtools

GWAS_TOOLS=~/data/expression_impute/summary-gwas-imputation/src;
METAXCAN=~/data/expression_impute/MetaXcan/software;
DATA=~/data/expression_impute/MetaXcan_sample_data;
OUTPUT=~/data/expression_impute/output;
python $GWAS_TOOLS/gwas_parsing.py \
    -gwas_file ~/pgc2017/adhd_eur_jun2017 \
    -liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
    -snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
    -output_column_map SNP variant_id \
    -output_column_map A2 non_effect_allele \
    -output_column_map A1 effect_allele \
    -output_column_map OR or \
    -output_column_map P pvalue \
    -output_column_map CHR chromosome \
    -output_column_map SE standard_error \
    --chromosome_format \
    -output_column_map BP position \
    --insert_value sample_size 55374 --insert_value n_cases 20183 \
    -output_order variant_id panel_variant_id chromosome position \
        effect_allele non_effect_allele pvalue or standard_error sample_size \
        n_cases \
    -output $OUTPUT/harmonized_gwas/ADHD_EUR_ADDITIVE.txt.gz
```

I got the two sample numbers from the README file, because without them the
script was crashing.

```bash
python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/ADHD_EUR_ADDITIVE.txt.gz \
    -parquet_genotype $DATA/reference_panel_1000G/chr1.variants.parquet \
    -parquet_genotype_metadata $DATA/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome 1 \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 10 \
    -sub_batch 0 \
    --standardise_dosages \
    --cache_variants \
    -output $OUTPUT/summary_imputation/ADHD_EUR_ADDITIVE_chr1_sb0_reg0.1_ff0.01_by_region.txt.gz
```

It's tripping up on the OR column... let me see if I do the quick step it works
better. In fact, I might be able to just run SPredXcan and let it handle the
harmonization in the background?

Didn't work... let's create the zscore column ourselves. Based on this:

https://huwenboshi.github.io/data%20management/2017/11/23/tips-for-formatting-gwas-summary-stats.html

I can just do:

```bash
echo "zscore" > zcolumn.txt;
tail -n +2 ~/pgc2017/adhd_eur_jun2017 | \
    awk '{a = log($7)/$8; printf("%0.4f\n", a)}' >> zcolumn.txt
paste ~/pgc2017/adhd_eur_jun2017 zcolumn.txt > adhd_eur_gwas_withZ.txt
```

Now, let's try this again:

```bash
GWAS_TOOLS=~/data/expression_impute/summary-gwas-imputation/src;
METAXCAN=~/data/expression_impute/MetaXcan-master/software;
DATA=~/data/expression_impute/MetaXcan_sample_data;
OUTPUT=~/data/expression_impute/output;
python $GWAS_TOOLS/gwas_parsing.py \
    -gwas_file ~/pgc2017/adhd_eur_jun2017 \
    -liftover $DATA/liftover/hg19ToHg38.over.chain.gz \
    -snp_reference_metadata $DATA/reference_panel_1000G/variant_metadata.txt.gz METADATA \
    -output_column_map SNP variant_id \
    -output_column_map A2 non_effect_allele \
    -output_column_map A1 effect_allele \
    -output_column_map OR or \
    -output_column_map P pvalue \
    -output_column_map CHR chromosome \
    -output_column_map SE standard_error \
    --chromosome_format \
    --enforce_numeric_columns \
    -output_column_map BP position \
    --insert_value sample_size 55374 --insert_value n_cases 20183 \
    -output_order variant_id panel_variant_id chromosome position \
        effect_allele non_effect_allele pvalue or zscore standard_error \
        sample_size n_cases \
    -output $OUTPUT/harmonized_gwas/ADHD_EUR_ADDITIVE_Z.txt.gz
```

And try the imputation command above again, but now using the file with Z.

That took less than 1min... I bet I can do a whole chromosome in less than
10min?

```bash
for c in {1..22}; do
python $GWAS_TOOLS/gwas_summary_imputation.py \
    -by_region_file $DATA/eur_ld.bed.gz \
    -gwas_file $OUTPUT/harmonized_gwas/ADHD_EUR_ADDITIVE_Z.txt.gz \
    -parquet_genotype $DATA/reference_panel_1000G/chr${c}.variants.parquet \
    -parquet_genotype_metadata $DATA/reference_panel_1000G/variant_metadata.parquet \
    -window 100000 \
    -parsimony 7 \
    -chromosome $c \
    -regularization 0.1 \
    -frequency_filter 0.01 \
    -sub_batches 1 \
    -sub_batch 0 \
    --standardise_dosages \
    --cache_variants \
    -output $OUTPUT/summary_imputation/ADHD_EUR_ADDITIVE_chr${c}_sb0_reg0.1_ff0.01_by_region.txt.gz;
done
```

Now that I fixed the Zscore issue, I can see that it takes over an hour per
chromosome (if I don't subbatch it), and it's completely overclocking my 32
CPUs. It's at 72 at some point, and taking a little less than 10Gb (if I want to
sbatch it).

So, I put a version of the script above in gwas_impute.sh and swarm ti like
this:

```bash
for c in {1..22}; do
    for b in {0..9}; do
        echo "bash /home/sudregp/data/expression_impute/gwas_impute.sh $c $b" >> swarm.impute;
    done;
done;
swarm -t 2 -g 10 --job-name gwas_imp --time 1:00:00 -f swarm.impute \
    --partition quick --logdir trash_impute -m python
```

# 2020-12-18 06:54:55

Now we gather everything:

```bash
python $GWAS_TOOLS/gwas_summary_imputation_postprocess.py \
    -gwas_file $OUTPUT/harmonized_gwas/ADHD_EUR_ADDITIVE_Z.txt.gz \
    -folder $OUTPUT/summary_imputation \
    -pattern ADHD_EUR_ADDITIVE* \
    -parsimony 7 \
    -output $OUTPUT/processed_summary_imputation/imputed_ADHD_EUR_ADDITIVE_Z.txt.gz
```

Had to make a few changes to the postprocessing script to not look for columns
that are not there. But now it seems to be working. Time to run S-PredXcan:

```bash
python $METAXCAN/SPrediXcan.py \
    --gwas_file  $OUTPUT/processed_summary_imputation/imputed_ADHD_EUR_ADDITIVE_Z.txt.gz \
    --snp_column panel_variant_id \
    --effect_allele_column effect_allele \
    --non_effect_allele_column non_effect_allele \
    --zscore_column zscore \
    --model_db_path $DATA/models/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.db \
    --covariance $DATA/models/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.txt.gz \
    --keep_non_rsid \
    --additional_output \
    --model_db_snp_key varID \
    --throw \
    --output_file $OUTPUT/spredixcan/eqtl/ADHD_ACC_MASHR.csv
```

And repeat for Caudate:

```bash
python $METAXCAN/SPrediXcan.py \
    --gwas_file  $OUTPUT/processed_summary_imputation/imputed_ADHD_EUR_ADDITIVE_Z.txt.gz \
    --snp_column panel_variant_id \
    --effect_allele_column effect_allele \
    --non_effect_allele_column non_effect_allele \
    --zscore_column zscore \
    --model_db_path $DATA/models/eqtl/mashr/mashr_Brain_Caudate_basal_ganglia.db \
    --covariance $DATA/models/eqtl/mashr/mashr_Brain_Caudate_basal_ganglia.txt.gz \
    --keep_non_rsid \
    --additional_output \
    --model_db_snp_key varID \
    --throw \
    --output_file $OUTPUT/spredixcan/eqtl/ADHD_Caudate_MASHR.csv
```

And let's see if this works for the EN models as well:

```bash
mydir=~/data/expression_impute;
python $METAXCAN/SPrediXcan.py \
    --gwas_file  $OUTPUT/processed_summary_imputation/imputed_ADHD_EUR_ADDITIVE_Z.txt.gz \
    --snp_column panel_variant_id \
    --effect_allele_column effect_allele \
    --non_effect_allele_column non_effect_allele \
    --zscore_column zscore \
    --model_db_path $mydir/elastic_net_models/en_Brain_Anterior_cingulate_cortex_BA24.db \
    --covariance $mydir/elastic_net_models/en_Brain_Anterior_cingulate_cortex_BA24.txt.gz \
    --keep_non_rsid \
    --additional_output \
    --model_db_snp_key varID \
    --throw \
    --output_file $OUTPUT/spredixcan/eqtl/ADHD_ACC_EN.csv
```

No RSIDs were found... I'll probably have to debug this one, but let's not do it
now. We have something for MASHR, which is the recommended model anyways. Let's
see how high this flies.

## Association with RNAseq PM

Let's see if the results are at all similar. We could try to check for
correlation, but maybe in the end we'll have to check overlaps under different
cutoffs:

```r
spred = read.csv('~/data/expression_impute/spredixcan/eqtl/ADHD_ACC_MASHR.csv')
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
both_res = merge(rnaseq_acc, spred, by.x='hgnc_symbol', by.y='gene_name',
                 all.x=F, all.y=F)
```

```
r$> cor.test(both_res$pvalue, both_res$P.Value, method='spearman')                          

        Spearman's rank correlation rho

data:  both_res$pvalue and both_res$P.Value
S = 1.9575e+11, p-value = 0.003134
alternative hypothesis: true rho is not equal to 0
sample estimates:
       rho 
-0.0288959 
```

So, this is encouraging. There is a correlation between the PM and the imputed
results. How does it look for Caudate?

```
r$> cor.test(both_res$pvalue, both_res$P.Value, method='spearman')                          

        Spearman's rank correlation rho

data:  both_res$pvalue and both_res$P.Value
S = 2.186e+11, p-value = 0.368
alternative hypothesis: true rho is not equal to 0
sample estimates:
         rho 
-0.008618003 
```

As expected (based on the PRS results), there is no correlation there.

It might also be informative to compute the gene overlap for different
thresholds:

```r
library(GeneOverlap)
spred = read.csv('~/data/expression_impute/spredixcan/eqtl/ADHD_ACC_MASHR.csv')
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
both_res = merge(rnaseq_acc, spred, by.x='hgnc_symbol', by.y='gene_name',
                 all.x=F, all.y=F)
thresh = c(.05, .01, .005, .001)
imp_nums = vector(length=length(thresh), mode='numeric')
rna_nums = vector(length=length(thresh), mode='numeric')
pvals = matrix(data=NA, nrow=length(thresh), ncol=length(thresh))
inter_nums = pvals
for (ti in 1:length(thresh)) {
    imp_genes = both_res[both_res$pvalue < thresh[ti], 'hgnc_symbol']
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

```
r$> inter_nums                                                                              
     [,1] [,2] [,3] [,4]
[1,]   58   15   10    0
[2,]    7    3    2    0
[3,]    2    1    0    0
[4,]    1    1    0    0

r$> pvals                                                                                   
          [,1]      [,2]      [,3] [,4]
[1,] 0.6643450 0.8213633 0.6968390    1
[2,] 0.9914974 0.8129635 0.7511822    1
[3,] 0.9972365 0.9076994 1.0000000    1
[4,] 0.8453085 0.4207225 1.0000000    1

r$> imp_nums                                                                                
[1] 818 249 155  52

r$> rna_nums                                                                                
[1] 779 185 100  23
```

My original idea was to make a nice plot here, but the results are not worth it.

At least the correlation result looks nice. Maybe I could do some WG analysis in
the imputed results? Maybe... I'll see what Philip says.

Let's spend some time in the brain phenotype then.

# 2020-12-21 19:12:00

Philip suggested I should use for rank the sign of the effect as well, similar
to what we did before. Let's see if that changes the results:

```r
spred = read.csv('~/data/expression_impute/spredixcan/eqtl/ADHD_ACC_MASHR.csv')
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
both_res = merge(rnaseq_acc, spred, by.x='hgnc_symbol', by.y='gene_name',
                 all.x=F, all.y=F)
```

There was no relationship if we took into consideration the up and down
regulation (ACC first):

```
r$> cor.test(sign(both_res$t)*both_res$pvalue, sign(both_res$zscore)*both_res$P.Value, method='spearman')                  

        Spearman's rank correlation rho

data:  sign(both_res$t) * both_res$pvalue and sign(both_res$zscore) * both_res$P.Value
S = 1.8963e+11, p-value = 0.7406
alternative hypothesis: true rho is not equal to 0
sample estimates:
        rho 
0.003238165

[...]

r$> cor.test(sign(both_res$t)*both_res$pvalue, sign(both_res$zscore)*both_res$P.Value, method='spearman')                  

        Spearman's rank correlation rho

data:  sign(both_res$t) * both_res$pvalue and sign(both_res$zscore) * both_res$P.Value
S = 2.2199e+11, p-value = 0.5368
alternative hypothesis: true rho is not equal to 0
sample estimates:
         rho 
-0.005894363 
```

