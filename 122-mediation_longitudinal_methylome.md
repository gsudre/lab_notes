# 2020-06-25 16:14:02

Just some notes for running the script. I ran 3 Xs in 4 minutes using a single
thread. If I split the files into groups of 15, that should give me 20 minutes
per Rscript call, which means 480 Xs every 20 minutes if using 32 cores. Even
the HI dataset, with 3.3K Xs, should take tops 3h. 

```bash
cd ~/data/longitudinal_methylome
split -l 15 IN_xs inatt_
ls -1 inatt_* > IN_sets.txt;
cat IN_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        INATT_ROC_methyl_sx_dti_82.csv {} IN_ms IN_ys res_{}.csv;
```

And for HI is very similar:

```bash
cd ~/data/longitudinal_methylome
split -l 15 HI_xs hi_
ls -1 hi_* > HI_sets.txt;
cat HI_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        HI_ROC_methyl_sx_dti_82.csv {} HI_ms HI_ys res_{}.csv;
```

Then, to compile: 

```bash
cd ~/data/longitudinal_methylome/results
head -n 1 res_inatt_bj.csv > inatt_compiled.csv;
for f in `ls -1 res_inatt_??.csv`; do tail -n +2 $f >> inatt_compiled.csv ; done
```

And I added the FDR, within and across Ms, in R:

```r
a = read.csv('~/data/longitudinal_methylome/results/inatt_compiled.csv')
a$tot_p_FDR = p.adjust(a$tot_p, method='fdr')
a$acme_p_FDR = p.adjust(a$acme_p, method='fdr')
a$ade_p_FDR = p.adjust(a$ade_p, method='fdr')
a$prop_p_FDR = p.adjust(a$prop_p, method='fdr')
idx = a$M=='AD_left_unc'
a$tot_p_FDR_withinM = NA
a$acme_p_FDR_withinM = NA
a$ade_p_FDR_withinM = NA
a$prop_p_FDR_withinM = NA
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
idx = a$M=='RD_left_unc'
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
write.csv(a, file='~/data/longitudinal_methylome/results/inatt_compiled_FDR.csv', row.names=F)
```

# 2020-06-29 20:27:08

I had ot make a few changes to the script to make it all run nicely. So, now
here's how to call it. I'll change it to only 8 cores so I can leave all
permutations running over night (HI each direction, plus without baseline)

```bash
cd ~/data/longitudinal_methylome
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} IN_ms_jhu IN_ys F \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_{}.csv;
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} IN_ms_jhu IN_ys T \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_{}.csv;
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv IN_ys IN_ms_jhu {} F \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_xSX_{}.csv;
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv IN_ys IN_ms_jhu {} T \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_xSX_{}.csv;
```

And its similar for HI, just that it takes 3 times longer...

```bash
cd ~/data/longitudinal_methylome
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_{}.csv;
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys T \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_{}.csv;
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv HI_ys HI_ms_jhu {} F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_xSX_{}.csv;
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv HI_ys HI_ms_jhu {} T \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_xSX_{}.csv;
```

And to compile it's similar:

```r
a = read.csv('~/data/longitudinal_methylome/results/hi_JHU_compiled.csv')
a$tot_p_FDR = p.adjust(a$tot_p, method='fdr')
a$acme_p_FDR = p.adjust(a$acme_p, method='fdr')
a$ade_p_FDR = p.adjust(a$ade_p, method='fdr')
a$prop_p_FDR = p.adjust(a$prop_p, method='fdr')
idx = a$M=='rd_11_rate'
a$tot_p_FDR_withinM = NA
a$acme_p_FDR_withinM = NA
a$ade_p_FDR_withinM = NA
a$prop_p_FDR_withinM = NA
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
idx = a$M=='rd_13_rate'
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
write.csv(a, file='~/data/longitudinal_methylome/results/hi_JHU_compiled_FDR.csv', row.names=F)
```

# 2020-06-30 21:33:59

Philip asked me to re-run the mediations for probe (X) to HI (.Y) through the
DTI. Up to 10000, maybe even 100k bootstraps? Based on the results, I think he
meant the option without baseline correction. I'll run that first, only 10K
bootstraps, and if it all looks good I can change the code to run 10K, because
it will likely need to be swarmed.

```bash
head -n 109 HI_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_10K_{}.csv;
# and also
tail -n 110 HI_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_10K_{}.csv;
```