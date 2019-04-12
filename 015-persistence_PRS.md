# 2019-04-12 13:22:32

The idea is to replicate an analysis Philip's collaborators ran to calculate PRS
for persistence of ADHD. They sent us a file with different weights to be used,
so in that sense it's quite similar to calculating PRS for ASD or SCZ. 

I just modified the file they sent to ranem the column BETA_diff to BETA and
P_diff to P. I also added the p-threshold they found to be best in their
results, just for kicks:

```bash
cd data/prs/geno3/imputation_1KG
mkdir persistence_genop05MAFbtp01rsbtp9
module load R
R --file=${HOME}/PRSice_v1.25/PRSice_v1.25.R -q --args  \
  plink  ${HOME}/PRSice_v1.25/plink_1.9_linux_160914 \
  base ~/pgc2017/ToShare_fused-gwas-adult-minus-child-impact-pgc_mod.txt  \
  target ../merged_1KG_biAllelicOnly_genop05MAFbtp01rsbtp9_renamed \
  report.individual.scores T \
  wd ./persistence_genop05MAFbtp01rsbtp9 \
  cleanup F \
  report.best.score.only F \
  covary F \
  fastscore T \
  barchart.levels 1e-5,1e-4,1e-3,1e-2,1e-1,5e-5,5e-4,5e-3,5e-2,2e-1,3e-1,4e-1,5e-1,0.00135
  
Rscript ~/research_code/lab_mgmt/collect_PRS.R
```

