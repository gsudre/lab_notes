# 2019-08-02 11:09:53

Philip asked me to generate PRS for the ABCD data we downloaded. Because we ran
our PRS on the imputed data, let's go ahead and merge that into PLINK format to
make things run faster.

```bash
# bw
cd /data/NCR_SBRB/ABCD/imputed_1KG/
for c in {1..22}; do
    plink --vcf chr${c}.dose.vcf.gz --biallelic-only strict \
        --make-bed --out chr${c};
done
for c in {2..22}; do
    echo chr${c} >> merge_list.txt;
done
plink --bfile chr1 --merge-list merge_list.txt  --make-bed --out merged_1KG
```

I'm actually not sure if I can load all this data into memory, so I might need
to split into chuncks of 1K subjects and then just concatenate them later...
let's see.


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