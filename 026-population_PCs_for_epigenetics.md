# 2019-07-02 10:53:51

The idea here is to construct new population PCs using only our cohort analyzed
for the methylation study.

Philip sent me a file called identified.csv by secure e-mail, so I'll need to
find which NSB was used for those MRNs to filter them out from the original
file, and then contruct the PCs using KING. The file from Philip had 248
samples.

I'm looking at
shaw/prs/FINAL_FILES_08242018/REGRESSION/kingpc_geno3_1KG_genop05MAFbtp01rsbtp9.csv
to get the MRN to NSB map. 

```bash
# bw sin
cd /lscratch/$SLURM_JOBID
cp -v ~/data/prs/geno3/merged_noDups_clean.* .
for n in `cat keep_ids.txt`; do
    grep " ${n} " merged_noDups_clean.fam >> keep_ids_with_fam.txt;
done;
plink --bfile merged_noDups_clean --keep keep_ids_with_fam.txt --make-bed --out methylated
/data/NCR_SBRB/software/KING/king -b methylated.bed --mds
```

Then, it's just a matter of merging them with the MRNs in R:

```r
a = read.csv('~/tmp/mrn_longnsb_n232.csv')
b = read.table('/lscratch/30548903/kingpc.ped')
colnames(b) = c('j1', 'long_NSB', 'j2', 'j3', 'j4', 'j5',
                sapply(1:20, function(x) sprintf('PC%02d', x)))
m = merge(a, b, by='long_NSB')
m[, 3:7] = NULL
write.csv(m, file='~/tmp/methylation_PCs_n232.csv', row.names=F)
```

