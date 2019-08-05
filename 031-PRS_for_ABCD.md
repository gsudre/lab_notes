# 2019-08-02 11:09:53

Philip asked me to generate PRS for the ABCD data we downloaded. Because we ran
our PRS on the imputed data, let's go ahead and merge that into PLINK format to
make things run faster.

```bash
# bw
module load plink
cd /data/NCR_SBRB/ABCD/imputed_1KG/
for c in {1..22}; do
    plink --vcf chr${c}.dose.vcf.gz --biallelic-only strict \
        --double-id --make-bed --out chr${c};
done
for c in {2..22}; do
    echo chr${c} >> merge_list.txt;
done
plink --bfile chr1 --merge-list merge_list.txt  --make-bed --out merged_1KG
```

I'm actually not sure if I can load all this data into memory, so I might need
to split into chuncks of 1K subjects and then just concatenate them later...
let's see.

But it was actually taking a long time just to convert from dosage files to
bim/bam. So, let's do this without imputation first:


```bash
cd /data/NCR_SBRB/ABCD/
mkdir prs_pgc2017
module load R
R --file=${HOME}/PRSice_v1.25/PRSice_v1.25.R -q --args  \
  plink  ${HOME}/PRSice_v1.25/plink_1.9_linux_160914 \
  base ~/pgc2017/adhd_jun2017  \
  target ../ABCD_release_2_genomics_smokescreen/ABCD_release_2.0_r2 \
  report.individual.scores T \
  wd ./prs_pgc2017 \
  cleanup F \
  report.best.score.only F \
  covary F \
  fastscore T \
  barchart.levels 1e-5,1e-4,1e-3,1e-2,1e-1,5e-5,5e-4,5e-3,5e-2,2e-1,3e-1,4e-1,5e-1
  
Rscript ~/research_code/lab_mgmt/collect_PRS.R
```

# 2019-08-05 12:21:59

This didn't quite work, and I'll need to reduce the files to see what's going
on:

```bash
# bw
cd /data/NCR_SBRB/ABCD/ABCD_release_2_genomics_smokescreen
head -n 100 ABCD_release_2.0_r2.fam > first100.txt
module load plink
plink -bfile ABCD_release_2.0_r2 --keep first100.txt --make-bed --out ABCD_release_2.0_r2_first100
module load R

cd /data/NCR_SBRB/ABCD/
R --file=${HOME}/PRSice_v1.25/PRSice_v1.25.R -q --args  \
  plink  ${HOME}/PRSice_v1.25/plink_1.9_linux_160914 \
  base ~/pgc2017/adhd_jun2017_fast  \
  target ../ABCD_release_2_genomics_smokescreen/ABCD_release_2.0_r2_first100 \
  report.individual.scores T \
  wd ./prs_pgc2017 \
  cleanup F \
  report.best.score.only F \
  covary F \
  fastscore T \
  barchart.levels 1e-5,1e-4,1e-3,1e-2,1e-1,5e-5,5e-4,5e-3,5e-2,2e-1,3e-1,4e-1,5e-1
```

But since we're doing this from scratch, let's use the new version of PRSice
from http://www.prsice.info:

```bash
cd /data/NCR_SBRB/ABCD/
Rscript /data/NCR_SBRB/software/PRSice_2.2.5/PRSice.R  \
    --prsice /data/NCR_SBRB/software/PRSice_2.2.5/PRSice_linux \
    --base ~/pgc2017/adhd_jun2017_fast  \
    --target /data/NCR_SBRB/ABCD/ABCD_release_2_genomics_smokescreen/ABCD_release_2.0_r2_first100 \
    --all-score \
    --lower 5e-08 --upper .5 --interval 5e-05 \
    --no-regress
```

Note that in this new version we get all the values in .all.score file that are
specified between lower and upper. The bar-levels variable is only there for
plotting!

So, now we should be able to run the entire thing:

```bash
cd /data/NCR_SBRB/ABCD/
Rscript /data/NCR_SBRB/software/PRSice_2.2.5/PRSice.R  \
    --prsice /data/NCR_SBRB/software/PRSice_2.2.5/PRSice_linux \
    --base ~/pgc2017/adhd_jun2017  \
    --target /data/NCR_SBRB/ABCD/ABCD_release_2_genomics_smokescreen/ABCD_release_2.0_r2 \
    --all-score \
    --lower 5e-08 --upper .5 --interval 5e-05 \
    --no-regress \
    --out ABCD_rel2_PRS_adhd_jun2017
```

Since we're here, we might as well compute the PRS for SCZ and ASD using rps.txt
and daner_AUT_meta14_WW_all.hg19.Mar2016_info_0.60_maf_0.05_release_Jun2017.tsv
(respectively), following the original Evernote note. Note that I used --stat OR
to make the ASD file run, as it had both BETA and OR columns.

I also did the PRS for WNH only, using adhd_eur_jun2017.

Then, we compute the usual population components (MDS, actually) using KING:

```bash
#bw
cd /data/NCR_SBRB/ABCD/
/data/NCR_SBRB/software/KING/king --mds --cpus 30 \
    -b /data/NCR_SBRB/ABCD/ABCD_release_2_genomics_smokescreen/ABCD_release_2.0_r2.bed
```

Finally, merge everything:

```r
# this takes a while because we're reading in TXT files!

# keep only some of the PRs columns that were created
mycols = c('FID', 'IID', 'X0.00010005', 'X0.00100005', 'X0.01', 'X0.1',
            'X5.005e.05', 'X0.00050005', 'X0.00500005', 'X0.0500001',
            'X0.5', 'X0.4', 'X0.3', 'X0.2')
a = read.table('/data/NCR_SBRB/ABCD/ABCD_rel2_PRS_adhd_jun2017.all.score', header=1)
af = a[, mycols]
b = read.table('/data/NCR_SBRB/ABCD/ABCD_rel2_PRS_adhd_eur_jun2017.all.score', header=1)
bf = b[, mycols]

m = merge(af, bf, by='IID', suffixes = c('.ADHD', '.ADHDeur'))

```