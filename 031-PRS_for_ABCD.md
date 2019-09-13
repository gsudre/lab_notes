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

*This took 8h for the SVD step!!!*

Finally, merge everything:

```r
# this takes a while because we're reading in TXT files!
a = read.table('/data/NCR_SBRB/ABCD/ABCD_rel2_PRS_adhd_jun2017.all.score', header=1)
b = read.table('/data/NCR_SBRB/ABCD/ABCD_rel2_PRS_adhd_eur_jun2017.all.score', header=1)
s = read.table('/data/NCR_SBRB/ABCD/SCZ.all.score', header=1)
d = read.table('/data/NCR_SBRB/ABCD/ASD.all.score', header=1)
pcs = read.table('/data/NCR_SBRB/ABCD/kingpc.ped')

# keep only some of the PRs columns that were created
mycols = c('IID', 'X0.00010005', 'X0.00100005', 'X0.01', 'X0.1',
            'X5.005e.05', 'X0.00050005', 'X0.00500005', 'X0.0500001',
            'X0.5', 'X0.4', 'X0.3', 'X0.2')
new_names = c('IID', sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                              .5, .4, .3, .2),
                            function(x) sprintf('ADHD_PRS%f', x)))
af = a[, mycols]
colnames(af) = new_names
bf = b[, mycols]
new_names = gsub('ADHD', x=new_names, 'ADHDeur')
colnames(bf) = new_names
df = d[, mycols]
new_names = gsub('ADHDeur', x=new_names, 'ASD')
colnames(df) = new_names
mycols = c('IID', 'X0.00010005', 'X0.00100005', 'X0.01', 'X0.1',
            'X5.005e.05', 'X0.00050005', 'X0.00500005', 'X0.0500001',
            'X0.5', 'X0.4001', 'X0.3002', 'X0.2002')
sf = s[,mycols]
new_names = gsub('ASD', x=new_names, 'SCZ')
colnames(sf) = new_names

m = merge(af, bf, by='IID')
m = merge(m, df, by='IID')
m = merge(m, sf, by='IID')

pcsf = pcs[, c(2, 7:26)]
new_names = c('IID', sapply(1:20, function(x) sprintf('PC%.2d', x)))
colnames(pcsf) = new_names
m = merge(m, pcsf, by='IID')

m = merge(pcs[, 1:2], m, by.x='V2', by.y='IID')
colnames(m)[1:2] = c('IID', 'FID')
write.csv(m, file='/data/NCR_SBRB/ABCD/merged_PRS_08062019.csv', row.names=F)
```

# 2019-09-13 13:08:02

Sam is having issues defining the ethnic groups based on the MDS components. So,
I'm running the pca through KING as well, just in case:

I'm also running MDS through plink, following the ENIGMA pipeline, to see if we
get anything better:

http://enigma.ini.usc.edu/wp-content/uploads/2012/07/ENIGMA2_1KGP_cookbook_v3.pdf

So, it goes like this:

```bash
cd /data/NCR_SBRB/ABCD/
wget "http://enigma.ini.usc.edu/wp-content/uploads/2012/07/HM3.bed.gz"
wget "http://enigma.ini.usc.edu/wp-content/uploads/2012/07/HM3.bim.gz"
wget "http://enigma.ini.usc.edu/wp-content/uploads/2012/07/HM3.fam.gz"
export datafileraw=/data/NCR_SBRB/ABCD/ABCD_release_2_genomics_smokescreen/ABCD_release_2.0_r2
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb --make-bed --out ${datafileraw}_filtered
gunzip HM3*.gz
export datafile=${datafileraw}_filtered
awk '{print $2}' HM3.bim > HM3.snplist.txt
plink --bfile ${datafile} --extract HM3.snplist.txt --make-bed --noweb --out local
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep -v ambig > local.snplist.txt
plink --bfile HM3 --extract local.snplist.txt --make-bed --noweb --out external
plink --bfile local --bmerge external.bed external.bim external.fam --make-bed --noweb --out HM3merge
plink --bfile local --flip HM3merge-merge.missnp --make-bed --noweb --out flipped
plink --bfile flipped --bmerge external.bed external.bim external.fam --make-bed --noweb --out HM3merge
plink --bfile HM3merge --cluster --mind .05 --mds-plot 10 --extract local.snplist.txt --noweb --out HM3mds
```