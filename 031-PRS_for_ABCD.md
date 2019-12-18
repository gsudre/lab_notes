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

BTW, using straight up PCA in KING is not converging, even after 24h.

# 2019-09-18 10:26:38

Philip asked to join the PNC and our cohort with the ABCD data before generating
the MDS components according to the ENIGMA paradigm. First, we need to figure
out what's going on with PNC genomics. The most recent version of the dbGap link is this one:
https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000607.v3.p2,
, according to the original PNC paper:
https://www.sciencedirect.com/science/article/pii/S1053811915002529?via%3Dihub#s0015

I downloaded the data a while ago, and it looks like all samples the first
version are in the second version. 

```bash
HG-02113362-DM4:GenotypeFiles sudregp$ wc -l */_*
    8715 phg000381.v2.NIMH_NeurodevelopmentalGenomics.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_file_7removed
    9269 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_use-nosra
   17984 total
HG-02113362-DM4:GenotypeFiles sudregp$ head phg000381.v2.NIMH_NeurodevelopmentalGenomics.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_file_7removed
600001676724    600001676724    0       0       600001676724    2       1       GO_Affy60
600003245643    600003245643    0       0       600003245643    1       1       GO_Affy60
600004963801    600004963801    0       0       600004963801    2       1       GO_Affy60
600005394890    600005394890    0       0       600005394890    1       1       GO_Affy60
600005726384    600005726384    0       0       600005726384    1       1       GO_Affy60
600008688531    600008688531    0       0       600008688531    1       1       GO_Affy60
600009963128    600009963128    0       0       600009963128    2       1       GO_Affy60
600010814166    600010814166    0       0       600010814166    2       1       GO_Affy60
600012815174    600012815174    0       0       600012815174    1       1       GO_Affy60
600013511285    600013511285    0       0       600013511285    1       1       GO_Affy60
HG-02113362-DM4:GenotypeFiles sudregp$ grep 600001676724 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_use-nosra
0       600001676724    0       0       600001676724    2       1       Array_SNP
HG-02113362-DM4:GenotypeFiles sudregp$ grep 600010814166 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.sample-info.MULTI/_fam_sub_mot_fat_sam_sex_con_use-nosra
0       600010814166    0       0       600010814166    2       1       Array_SNP
```

I didn't test all of them, but the numbers are somewhat indicative, and the 9269
number is a closer match to what the wbesite says they have. Also, closer to 500
or so subjects they say they added from version 1 to 2.

Now, we just need to get the PLINK files and start playing. But it does seem
like there were many different genotyping arrays, which will make things quite
interesting in merging the datasets.

```
HG-02113362-DM4:GenotypeFiles sudregp$ wc -l */*fam
     722 phg000381.v2.NIMH_NeurodevelopmentalGenomics.genotype-calls-matrixfmt.Axiom.c1.GRU-NPU/GO_Axiom.fam
      66 phg000381.v2.NIMH_NeurodevelopmentalGenomics.genotype-calls-matrixfmt.Genome-Wide_Human_SNP_Array_6.0.c1.GRU-NPU/GO_Affy60.fam
    3802 phg000381.v2.NIMH_NeurodevelopmentalGenomics.genotype-calls-matrixfmt.Human610-Quadv1_B.c1.GRU-NPU/GO_Quad_5removed.fam
     555 phg000381.v2.NIMH_NeurodevelopmentalGenomics.genotype-calls-matrixfmt.HumanHap550_v1.c1.GRU-NPU/GO_v1_1removed.fam
    1913 phg000381.v2.NIMH_NeurodevelopmentalGenomics.genotype-calls-matrixfmt.HumanHap550_v3.c1.GRU-NPU/GO_v3_1removed.fam
    1657 phg000381.v2.NIMH_NeurodevelopmentalGenomics.genotype-calls-matrixfmt.HumanOmniExpress.c1.GRU-NPU/GO_Omni.fam
     225 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.Axiom.c1.GRU-NPU/GO_AxiomTx.fam
      40 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.Axiom.c1.GRU-NPU/GO_Axiom_set2.fam
      17 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.BDCHP-1X10-HUMANHAP550.c1.GRU-NPU/GO_v1set2.fam
     141 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.Human1M-Duov3_B.c1.GRU-NPU/GO_1MDuo.fam
      40 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.Human610-Quadv1_B.c1.GRU-NPU/GO_Quadset2.fam
      31 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.HumanHap550_v3.c1.GRU-NPU/GO_v3set2.fam
      37 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.HumanOmniExpress-12v1_A.c1.GRU-NPU/GO_Omniset2.fam
      18 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.HumanOmniExpress-12v1_B.c1.GRU-NPU/GO_OMNI12v11.fam
       3 phg000661.v1.NIMH_NeurodevelopmentalGenomics_v2.genotype-calls-matrixfmt.HumanOmniExpressExome-8v1_A.c1.GRU-NPU/GO_OEE.fam
    9267 total
```

My plan is to modify the ENIGMA protocol slightly to account for the different
datasets. HEre, they will be ours, ABCD, and all of PNC. I'll extract just the
HM3 SNPS from each of them, and merge only afterwards.

http://enigma.ini.usc.edu/wp-content/uploads/2012/07/ENIGMA2_1KGP_cookbook_v3.pdf

```bash
# interactive
cd /data/NCR_SBRB/combined_genetics
module load plink

[sudregp@cn1741 combined_genetics]$ ls -1 *bed
ABCD_release_2.0_r2.bed
CIDR_1to22.bed
GO_1MDuo.bed
GO_Affy60.bed
GO_Axiom.bed
GO_Axiom_set2.bed
GO_AxiomTx.bed
GO_OEE.bed
GO_OMNI12v11.bed
GO_Omni.bed
GO_Omniset2.bed
GO_Quad_5removed.bed
GO_Quadset2.bed
GO_v1_1removed.bed
GO_v1set2.bed
GO_v3_1removed.bed
GO_v3set2.bed
HM3.bed
Shaw01_2017_1to22.bed
Shaw02_2017_1to22.bed
Shaw03_1to22.bed
Shaw03_2017_1to22.bed
Shaw04_1to22.bed
Shaw04_2017_1to22.bed
Shaw05_2017_1to22.bed
Shaw06_2017_1to22.bed
Shaw07_2017_1to22.bed
twins_1to22.bed

awk '{print $2}' HM3.bim > HM3.snplist.txt;
for data in `ls *bed`; do
    datafileraw=`basename -s .bed $data`;
    if [ $datafileraw != 'HM3' ]; then
        plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb \
            --make-bed --out ${datafileraw}_filtered
        datafile=${datafileraw}_filtered;
        plink --bfile ${datafile} --extract HM3.snplist.txt --make-bed --noweb --out ${datafile}_local;
        awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' ${datafile}.bim | grep -v ambig > ${datafile}_local.snplist.txt;
    fi;
done
[sudregp@cn1741 combined_genetics]$ wc -l *local.snplist.txt
   334195 ABCD_release_2.0_r2_filtered_local.snplist.txt
   673148 CIDR_1to22_filtered_local.snplist.txt
  1029421 GO_1MDuo_filtered_local.snplist.txt
   705032 GO_Affy60_filtered_local.snplist.txt
   459738 GO_Axiom_filtered_local.snplist.txt
   330049 GO_Axiom_set2_filtered_local.snplist.txt
   549524 GO_AxiomTx_filtered_local.snplist.txt
   511493 GO_OEE_filtered_local.snplist.txt
   628580 GO_OMNI12v11_filtered_local.snplist.txt
   653198 GO_Omni_filtered_local.snplist.txt
   687381 GO_Omniset2_filtered_local.snplist.txt
   459395 GO_Quad_5removed_filtered_local.snplist.txt
   577678 GO_Quadset2_filtered_local.snplist.txt
   516440 GO_v1_1removed_filtered_local.snplist.txt
   517539 GO_v1set2_filtered_local.snplist.txt
   475132 GO_v3_1removed_filtered_local.snplist.txt
   539374 GO_v3set2_filtered_local.snplist.txt
   680975 Shaw01_2017_1to22_filtered_local.snplist.txt
   637883 Shaw02_2017_1to22_filtered_local.snplist.txt
   633280 Shaw03_1to22_filtered_local.snplist.txt
   675513 Shaw03_2017_1to22_filtered_local.snplist.txt
   694192 Shaw04_1to22_filtered_local.snplist.txt
   676186 Shaw04_2017_1to22_filtered_local.snplist.txt
   661454 Shaw05_2017_1to22_filtered_local.snplist.txt
   421349 Shaw06_2017_1to22_filtered_local.snplist.txt
   666876 Shaw07_2017_1to22_filtered_local.snplist.txt
   645735 twins_1to22_filtered_local.snplist.txt

```

There's lots of variability in the number of SNPs in each dataset, even after
extracting only the ones that are in HM3. So, let's grab only the intersection
of the lists, and then re-extract everything, including HM3.

```bash
cp ABCD_release_2.0_r2_filtered_local.snplist.txt intersect.snplist.txt
for f in `ls *_filtered_local.snplist.txt`; do 
    echo $f;
    # combine the two files, sort them, and then display only the lines that appear more than once: that is, the ones that appear in both files.
    sort intersect.snplist.txt $f | uniq -d > junk.txt;
    cp junk.txt intersect.snplist.txt;
    wc -l intersect.snplist.txt;
done
```

```
ABCD_release_2.0_r2_filtered_local.snplist.txt
331992 intersect.snplist.txt
CIDR_1to22_filtered_local.snplist.txt
68279 intersect.snplist.txt
GO_1MDuo_filtered_local.snplist.txt
57286 intersect.snplist.txt
GO_Affy60_filtered_local.snplist.txt
22407 intersect.snplist.txt
GO_Axiom_filtered_local.snplist.txt
10295 intersect.snplist.txt
GO_Axiom_set2_filtered_local.snplist.txt
7729 intersect.snplist.txt
GO_AxiomTx_filtered_local.snplist.txt
6028 intersect.snplist.txt
GO_OEE_filtered_local.snplist.txt
4383 intersect.snplist.txt
GO_OMNI12v11_filtered_local.snplist.txt
4302 intersect.snplist.txt
GO_Omni_filtered_local.snplist.txt
4297 intersect.snplist.txt
GO_Omniset2_filtered_local.snplist.txt
4269 intersect.snplist.txt
GO_Quad_5removed_filtered_local.snplist.txt
2792 intersect.snplist.txt
GO_Quadset2_filtered_local.snplist.txt
2788 intersect.snplist.txt
GO_v1_1removed_filtered_local.snplist.txt
2625 intersect.snplist.txt
GO_v1set2_filtered_local.snplist.txt
2549 intersect.snplist.txt
GO_v3_1removed_filtered_local.snplist.txt
2523 intersect.snplist.txt
GO_v3set2_filtered_local.snplist.txt
2497 intersect.snplist.txt
Shaw01_2017_1to22_filtered_local.snplist.txt
2481 intersect.snplist.txt
Shaw02_2017_1to22_filtered_local.snplist.txt
2476 intersect.snplist.txt
Shaw03_1to22_filtered_local.snplist.txt
2470 intersect.snplist.txt
Shaw03_2017_1to22_filtered_local.snplist.txt
2469 intersect.snplist.txt
Shaw04_1to22_filtered_local.snplist.txt
2465 intersect.snplist.txt
Shaw04_2017_1to22_filtered_local.snplist.txt
2463 intersect.snplist.txt
Shaw05_2017_1to22_filtered_local.snplist.txt
2406 intersect.snplist.txt
Shaw06_2017_1to22_filtered_local.snplist.txt
1533 intersect.snplist.txt
Shaw07_2017_1to22_filtered_local.snplist.txt
1528 intersect.snplist.txt
twins_1to22_filtered_local.snplist.txt
1528 intersect.snplist.txt
```

The number of SNPs went down really fast. Not sure if we can run anything with
this... but let's see.

```bash
for data in `ls *_filtered_local.bed`; do
    datafile=`basename -s .bed $data`;
    plink --bfile ${datafile} --extract intersect.snplist.txt --make-bed --noweb --out ${datafile}_tiny;
    echo "${datafile}_tiny.bed ${datafile}_tiny.bim ${datafile}_tiny.fam" >> merge_list.txt
done
plink --bfile HM3 --extract intersect.snplist.txt --make-bed --noweb --out HM3_tiny;
```

I'm running into several issues of flipping variants, so I'll have to merge them
slowly...

```
datafile=ABCD_release_2.0_r2_filtered_local_tiny;
cnt=1;
plink --bfile HM3_tiny -bmerge ${datafile}.bed ${datafile}.bim ${datafile}.fam --make-bed --noweb --out HM3merge${cnt};

datafile=CIDR_1to22_filtered_local_tiny;
plink --bfile HM3merge${cnt} -bmerge ${datafile}.bed ${datafile}.bim ${datafile}.fam --make-bed --noweb --out HM3merge$(( cnt + 1))
plink --bfile $datafile --flip HM3merge$(( cnt + 1))-merge.missnp --make-bed --noweb --out ${datafile}_flipped
plink --bfile HM3merge${cnt} -bmerge ${datafile}_flipped.bed ${datafile}_flipped.bim ${datafile}_flipped.fam --make-bed --noweb --out HM3merge$(( cnt + 1))
let cnt=$cnt+1;

datafile=GO_1MDuo_filtered_local_tiny;
plink --bfile HM3merge${cnt} -bmerge ${datafile}.bed ${datafile}.bim ${datafile}.fam --make-bed --noweb --out HM3merge$(( cnt + 1))
plink --bfile $datafile --flip HM3merge$(( cnt + 1))-merge.missnp --make-bed --noweb --out ${datafile}_flipped
plink --bfile HM3merge${cnt} -bmerge ${datafile}_flipped.bed ${datafile}_flipped.bim ${datafile}_flipped.fam --make-bed --noweb --out HM3merge$(( cnt + 1))
let cnt=$cnt+1;

```

I don't like where this is going... maybe the best approach will be to do the
intersection of each one with HM3?

At least for within PNC, I'll need to find a way to merge them. I could impute
each one first, and then merge? That would take care of the flipping. What's the
intersection of just PNC?

```bash
cd /data/NCR_SBRB/PNC_genetics
for f in `ls *bim`; do
    fname=`basename -s .bim $f`;
    awk '{print $2}' $f > $fname.snplist.txt;
done
cp GO_1MDuo.snplist.txt intersect.snplist.txt;
for f in `ls *.snplist.txt`; do
    echo $f;
    # combine the two files, sort them, and then display only the lines that appear more than once: that is, the ones that appear in both files.
    sort intersect.snplist.txt $f | uniq -d > junk.txt;
    cp junk.txt intersect.snplist.txt;
    wc -l intersect.snplist.txt;
done
```

```
GO_1MDuo.snplist.txt
1199187 intersect.snplist.txt
GO_Affy60.snplist.txt
308753 intersect.snplist.txt
GO_Axiom_set2.snplist.txt
92294 intersect.snplist.txt
GO_Axiom.snplist.txt
92294 intersect.snplist.txt
GO_AxiomTx.snplist.txt
23952 intersect.snplist.txt
GO_OEE.snplist.txt
16221 intersect.snplist.txt
GO_OMNI12v11.snplist.txt
16056 intersect.snplist.txt
GO_Omniset2.snplist.txt
16056 intersect.snplist.txt
GO_Omni.snplist.txt
16056 intersect.snplist.txt
GO_Quad_5removed.snplist.txt
10547 intersect.snplist.txt
GO_Quadset2.snplist.txt
10547 intersect.snplist.txt
GO_v1_1removed.snplist.txt
10061 intersect.snplist.txt
GO_v1set2.snplist.txt
10061 intersect.snplist.txt
GO_v3_1removed.snplist.txt
10059 intersect.snplist.txt
GO_v3set2.snplist.txt
10059 intersect.snplist.txt
```

Even within PNC we're losing lots of SNPs. I'll go ahead and run the ENIGMA
protocol just on these, and we can see what we get in the end. Merging will be a
nightmare, but let's give it a try. 

```bash
cd /data/NCR_SBRB/PNC_genetics
cp ../ABCD/HM3.??? .
awk '{print $2}' HM3.bim > HM3.snplist.txt;
for data in `ls *bed`; do
    datafileraw=`basename -s .bed $data`;
    if [ $datafileraw != 'HM3' ]; then
        plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb \
            --make-bed --out ${datafileraw}_filtered
        datafile=${datafileraw}_filtered;
        plink --bfile ${datafile} --extract HM3.snplist.txt --make-bed --noweb --out ${datafile}_local;
        awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' ${datafile}.bim | grep -v ambig > ${datafile}_local.snplist.txt;
    fi;
done

cp GO_1MDuo_filtered_local.snplist.txt intersect.snplist.txt
for f in `ls *_filtered_local.snplist.txt`; do
    echo $f;
    # combine the two files, sort them, and then display only the lines that appear more than once: that is, the ones that appear in both files.
    sort intersect.snplist.txt $f | uniq -d > junk.txt;
    cp junk.txt intersect.snplist.txt;
    wc -l intersect.snplist.txt;
done
```

```
GO_1MDuo_filtered_local.snplist.txt
1029421 intersect.snplist.txt
GO_Affy60_filtered_local.snplist.txt
264439 intersect.snplist.txt
GO_Axiom_filtered_local.snplist.txt
75443 intersect.snplist.txt
GO_Axiom_set2_filtered_local.snplist.txt
53395 intersect.snplist.txt
GO_AxiomTx_filtered_local.snplist.txt
15144 intersect.snplist.txt
GO_OEE_filtered_local.snplist.txt
8051 intersect.snplist.txt
GO_OMNI12v11_filtered_local.snplist.txt
7793 intersect.snplist.txt
GO_Omni_filtered_local.snplist.txt
7721 intersect.snplist.txt
GO_Omniset2_filtered_local.snplist.txt
7672 intersect.snplist.txt
GO_Quad_5removed_filtered_local.snplist.txt
4324 intersect.snplist.txt
GO_Quadset2_filtered_local.snplist.txt
4319 intersect.snplist.txt
GO_v1_1removed_filtered_local.snplist.txt
4081 intersect.snplist.txt
GO_v1set2_filtered_local.snplist.txt
3964 intersect.snplist.txt
GO_v3_1removed_filtered_local.snplist.txt
3909 intersect.snplist.txt
GO_v3set2_filtered_local.snplist.txt
3873 intersect.snplist.txt
```

Again, not much overlap. Actually mirroring what we had when combining all
datasets. I actually re-ran the numbers and without any filtering the
intersection between PNC and HM# goes from 10: to 8K SNPs. So, the other 5K we
are losing above is because of QC. But I'm not sure if using QC in those small
sets is really that fair. So, let's merge the datasets first, 

OK, let's do this slowly then:

```bash
[sudregp@cn1741 PNC_genetics]$ ls -1 *bim
GO_1MDuo.bim
GO_Affy60.bim
GO_Axiom.bim
GO_Axiom_set2.bim
GO_AxiomTx.bim
GO_OEE.bim
GO_OMNI12v11.bim
GO_Omni.bim
GO_Omniset2.bim
GO_Quad_5removed.bim
GO_Quadset2.bim
GO_v1_1removed.bim
GO_v1set2.bim
GO_v3_1removed.bim
GO_v3set2.bim
```

```bash
#starting with Illumina arrays
plink --bfile GO_Omni -bmerge GO_Omniset2.bed GO_Omniset2.bim GO_Omniset2.fam --make-bed --noweb --out mergeOmni
plink --bfile GO_v1_1removed -bmerge GO_v1set2.bed GO_v1set2.bim GO_v1set2.fam --make-bed --noweb --out mergev1
plink --bfile GO_v3_1removed -bmerge GO_v3set2.bed GO_v3set2.bim GO_v3set2.fam --make-bed --noweb --out mergev3
plink --bfile GO_Quad_5removed -bmerge GO_Quadset2.bed GO_Quadset2.bim GO_Quadset2.fam --make-bed --noweb --out mergeQuad
plink --bfile mergeQuad -bmerge mergeOmni.bed mergeOmni.bim mergeOmni.fam --make-bed --noweb --out merge1
plink --bfile merge1 -bmerge mergev1.bed mergev1.bim mergev1.fam --make-bed --noweb --out merge2
plink --bfile merge2 -bmerge GO_OMNI12v11.bed GO_OMNI12v11.bim GO_OMNI12v11.fam --make-bed --noweb --out merge3
plink --bfile merge3 -bmerge GO_OEE.bed GO_OEE.bim GO_OEE.fam --make-bed --noweb --out merge4
plink --bfile merge4 -bmerge GO_1MDuo.bed GO_1MDuo.bim GO_1MDuo.fam --make-bed --noweb --out merge5
plink --bfile merge5 -bmerge mergev3.bed mergev3.bim mergev3.fam --make-bed --noweb --out merge6

# this takes care of all Illumina datasets from what I can see... now, to Affymetrix
plink --bfile GO_Axiom -bmerge GO_Axiom_set2.bed GO_Axiom_set2.bim GO_Axiom_set2.fam --make-bed --noweb --out mergeAxiom
plink --bfile mergeAxiom -bmerge GO_Affy60.bed GO_Affy60.bim GO_Affy60.fam --make-bed --noweb --out merge7
plink --bfile merge7 -bmerge GO_AxiomTx.bed GO_AxiomTx.bim GO_AxiomTx.fam --make-bed --noweb --out merge8
plink --bfile merge7 --flip merge8-merge.missnp --make-bed --noweb --out merge7_flipped
plink --bfile merge7_flipped -bmerge GO_AxiomTx.bed GO_AxiomTx.bim GO_AxiomTx.fam --make-bed --noweb --out merge8

# the big merge across platforms
plink --bfile merge8 -bmerge merge6.bed merge6.bim merge6.fam --make-bed --noweb --out merged
plink --bfile merge8 --flip merged-merge.missnp --make-bed --noweb --out merge8_flipped
plink --bfile merge8_flipped -bmerge merge6.bed merge6.bim merge6.fam --make-bed --noweb --out merged
```

The big merge still didn't work... maybe if I do the HM3 part first?

```bash
cd /data/NCR_SBRB/PNC_genetics
cp ../ABCD/HM3.??? .
awk '{print $2}' HM3.bim > HM3.snplist.txt;

export datafileraw=merge6
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb --make-bed --out ${datafileraw}_filtered
export datafile=${datafileraw}_filtered
plink --bfile ${datafile} --extract HM3.snplist.txt --make-bed --noweb --out illumina_local
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep -v ambig > illumina_local.snplist.txt

export datafileraw=merge8
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb --make-bed --out ${datafileraw}_filtered
export datafile=${datafileraw}_filtered
plink --bfile ${datafile} --extract HM3.snplist.txt --make-bed --noweb --out affy_local
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep -v ambig > affy_local.snplist.txt

plink --bfile HM3 --extract illumina_local.snplist.txt --make-bed --noweb --out external_illumina
plink --bfile illumina_local --bmerge external_illumina.bed external_illumina.bim external_illumina.fam --make-bed --noweb --out HM3merge1
plink --bfile illumina_local --flip HM3merge1-merge.missnp --make-bed --noweb --out illumina_flipped
plink --bfile illumina_flipped --bmerge external_illumina.bed external_illumina.bim external_illumina.fam --make-bed --noweb --out HM3merge1
plink --bfile HM3merge --cluster --mind .05 --mds-plot 10 --extract local.snplist.txt --noweb --out HM3mds
```

Having problems here again. I think the best solution would be to impute within
platform, and then do this. But for now we'll just use Veera's PCs for PNC, and
I'll recompute ours, doing age restriction.

```bash
cd /data/NCR_SBRB/NCR_genetics
cp -v ~/data/prs/geno3/*_1to22.??? .
[sudregp@cn1741 NCR_genetics]$ ls -1 *_1to22.bed
CIDR_1to22.bed
Shaw01_2017_1to22.bed
Shaw02_2017_1to22.bed
Shaw03_1to22.bed
Shaw03_2017_1to22.bed
Shaw04_1to22.bed
Shaw04_2017_1to22.bed
Shaw05_2017_1to22.bed
Shaw06_2017_1to22.bed
Shaw07_2017_1to22.bed
twins_1to22.bed
```

```bash
plink --bfile Shaw01_2017_1to22 -bmerge Shaw02_2017_1to22.bed Shaw02_2017_1to22.bim Shaw02_2017_1to22.fam --make-bed --noweb --out mergeShaw2017_1
plink --bfile mergeShaw2017_1 -bmerge Shaw03_2017_1to22.bed Shaw03_2017_1to22.bim Shaw03_2017_1to22.fam --make-bed --noweb --out mergeShaw2017_2
plink --bfile mergeShaw2017_2 -bmerge Shaw04_2017_1to22.bed Shaw04_2017_1to22.bim Shaw04_2017_1to22.fam --make-bed --noweb --out mergeShaw2017_3
plink --bfile mergeShaw2017_3 -bmerge Shaw05_2017_1to22.bed Shaw05_2017_1to22.bim Shaw05_2017_1to22.fam --make-bed --noweb --out mergeShaw2017_4
plink --bfile mergeShaw2017_4 -bmerge Shaw06_2017_1to22.bed Shaw06_2017_1to22.bim Shaw06_2017_1to22.fam --make-bed --noweb --out mergeShaw2017_5
plink --bfile mergeShaw2017_5 -bmerge Shaw07_2017_1to22.bed Shaw07_2017_1to22.bim Shaw07_2017_1to22.fam --make-bed --noweb --out mergeShaw2017_6
plink --bfile Shaw04_1to22 -bmerge Shaw03_1to22.bed Shaw03_1to22.bim Shaw03_1to22.fam --make-bed --noweb --out mergeShaw
plink --bfile mergeShaw -bmerge mergeShaw2017_6.bed mergeShaw2017_6.bim mergeShaw2017_6.fam --make-bed --noweb --out merge1
plink --bfile merge1 -bmerge twins_1to22.bed twins_1to22.bim twins_1to22.fam --make-bed --noweb --out merge2
plink --bfile merge2 -bmerge CIDR_1to22.bed CIDR_1to22.bim CIDR_1to22.fam --make-bed --noweb --out merged

export datafileraw=merged
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb --make-bed --out ${datafileraw}_filtered
cp ../ABCD/HM3.??? .
export datafile=${datafileraw}_filtered
awk '{print $2}' HM3.bim > HM3.snplist.txt
plink --bfile ${datafile} --extract HM3.snplist.txt --make-bed --noweb --out local
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep -v ambig > local.snplist.txt
plink --bfile HM3 --extract local.snplist.txt --make-bed --noweb --out external
plink --bfile local --bmerge external.bed external.bim external.fam --make-bed --noweb --out HM3merge
plink --bfile local --flip HM3merge-merge.missnp --make-bed --noweb --out flipped
plink --bfile flipped --bmerge external.bed external.bim external.fam --make-bed --noweb --out HM3merge

# at this stage I have everyone in H3merge. Now it's just a matter of removing any NSBs for people with age above 22
plink --bfile HM3merge --keep keep_younger_22.txt --make-bed --out HM3merge_LT22yo
plink --bfile HM3merge_LT22yo --cluster --mind .05 --mds-plot 20 --extract local.snplist.txt --noweb --out HM3_LT22yo_mds
plink --bfile HM3merge --keep keep_younger_22_noDups.txt --make-bed --out HM3merge_LT22yo_noDups
plink --bfile HM3merge_LT22yo_noDups --cluster --mind .05 --mds-plot 20 --extract local.snplist.txt --noweb --out HM3_LT22yo_noDups_mds
```

# 2019-10-17 15:17:53

Generating similar file for Kathryn's analysis, but here we can just pick up
from where we left off in Sam's analysis:

```bash
module load plink
cd /data/NCR_SBRB/NCR_genetics

plink --bfile HM3merge --keep keep_younger_22.txt --make-bed --out HM3merge_LT22yo
plink --bfile HM3merge_LT22yo --cluster --mind .05 --mds-plot 20 --extract local.snplist.txt --noweb --out HM3_LT22yo_mds
plink --bfile HM3merge --keep keep_younger_22_noDups.txt --make-bed --out HM3merge_LT22yo_noDups
plink --bfile HM3merge_LT22yo_noDups --cluster --mind .05 --mds-plot 20 --extract local.snplist.txt --noweb --out HM3_LT22yo_noDups_mds
```

# 2019-11-27 20:09:30

Going back to this issue with PNC, I'll clearly have to understand a bit more
about how Illumina and Affymetrix data are coding their alleles, before solving
this. For example, here I continued the analysis from above, but tried matching
only the SNPs that were common to both platforms:

```bash
# interactive
cd /data/NCR_SBRB/PNC_genetics
plink --bfile merge6 --write-snplist
plink --bfile merge8 --extract plink.snplist --make-bed --out merge8_filtered
plink --bfile merge8 --write-snplist
plink --bfile merge6 --extract plink.snplist --make-bed --out merge6_filtered
plink --bfile merge8_filtered -bmerge merge6_filtered.bed merge6_filtered.bim merge6_filtered.fam --make-bed --noweb --out merged
```

Doing the --flip version didn't help either, so that's not it. For example,
here's the output:

```
596900 markers loaded from merge8_filtered.bim.
596900 markers to be merged from merge6_filtered.bim.
Of these, 0 are new, while 596900 are present in the base dataset.
102929 more multiple-position warnings: see log file.
Error: 595116 variants with 3+ alleles present.
```

So, almost all SNPs are coded wrong. For example:

```
[sudregp@cn1864 PNC_genetics]$ head merged-merge.missnp 
rs1000002
rs1000003
rs10000037
rs10000041
rs10000042
rs10000073
rs10000081
rs10000085
rs10000091
rs10000092
[sudregp@cn1864 PNC_genetics]$ grep rs10000037 merge6_filtered.bim
4       rs10000037      0       38600725        1       2
[sudregp@cn1864 PNC_genetics]$ grep rs10000037 merge8_filtered.bim
4       rs10000037      0       38600725        A       G
```

But it's not just an issue of numbers vs letters, because if I try converting it
using --aleleACGT they don't match either. Need to read more about it.

# 2019-12-18 15:22:33

I downloaded the new ABCD data, so let's impute and run PRS again on it. Before
we do any of that, let's check on sgender assignments of this new dataset, and
take the usual ShapeIt steps:

```bash
# sinteractive
[sudregp@cn2167 v201]$ pwd
/data/NCR_SBRB/ABCD/v201
[sudregp@cn2167 v201]$ wc -l ABCD_release_2.0.1_r1.fam
10627 ABCD_release_2.0.1_r1.fam
```

OK, doing good in terms of using the correct file, based on the README for
release 2.0.1. Now, we need to update sex and remove problematic IDs. In fact,
let's just go ahead and use the exact same pipeline we used for our own data, in
note 59:

```bash
# used abcddemo01.txt to create the update file
plink --bfile ABCD_release_2.0.1_r1 --update-sex update_sex.txt --make-bed \
  --out ABCD_2.0.1_sex;
plink --bfile ABCD_2.0.1_sex --check-sex;
```

PLINK found problems determining the sex of 216 samples.  Since I cannot tell
where the error actually is, I'll just go ahead and remove those IDs. 

```bash
grep "PROBLEM" plink.sexcheck | awk '{print $1, $2}' >> failed_sex_check.txt;
plink --bfile ABCD_2.0.1_sex --remove failed_sex_check.txt --make-bed \
  --out ABCD_2.0.1_sexClean;
```

Time to check for identical samples. Here I'll need further investigation to see
if they're indeed twins, as I don't expect subjects to have more than one
sample. In any case, whether we'll use twins later will depend on the analysis.

```bash
plink --bfile ABCD_2.0.1_sexClean --genome
awk '{if ($10 > .95) print $1, $2, $3, $4}' plink.genome | wc -l
```

**OK, so there are 364 pairs of identical samples. Yes, they could just be twins,
so we'll need to be careful when doing further analysis.**

From now on, it's just a matter of following the ENIGMA imputation protocol. But
I'll just go ahead and copy what I did in note 59 because ShapeIt seems to be
quite necessary for the alignment.

```bash
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.bed.gz"
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.bim.gz"
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.fam.gz"
# Filter SNPs out from your dataset which do not meet Quality Control criteria
# (Minor Allele Frequency < 0.01; Genotype Call Rate < 95%; Hardy­Weinberg
# Equilibrium < 1x10­6)
export datafileraw=ABCD_2.0.1_sexClean
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb \
      --make-bed --out ${datafileraw}_filtered
# Unzip the HM3 genotypes. Prepare the HM3 and the raw genotype data by
# extracting only snps that are in common between the two genotype data sets
# this avoids exhausting the system memory. We are also removing the strand
# ambiguous snps from the genotyped data set to avoid strand mismatch among
# these snps. Your genotype files should be filtered to remove markers which
# do not satisfy the quality control criteria above.
gunzip HM3_b37*.gz
export datafile=${datafileraw}_filtered
awk '{print $2}' HM3_b37.bim > HM3_b37.snplist.txt
plink --bfile ${datafile} --extract HM3_b37.snplist.txt --make-bed --noweb --out local
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep -v ambig > local.snplist.txt
plink --bfile HM3_b37 --extract local.snplist.txt --make-bed --noweb --out external
# Merge the two sets of plink files. In merging the two files plink will check
# for strand differences. If any strand differences are found plink will crash
# with the following error (ERROR: Stopping due to mis­matching SNPs - check +/­
# strand?). Ignore warnings regarding different physical positions
plink --bfile local --bmerge external.bed external.bim external.fam \
  --make-bed --noweb --out HM3_b37merge
# got the error
plink --bfile local --flip HM3_b37merge-merge.missnp --make-bed --noweb \
  --out flipped
plink --bfile flipped --bmerge external.bed external.bim external.fam \
  --make-bed --noweb --out HM3_b37merge
# running MDS analysis... switching to 10 dimensions to conform to old analysis
plink --bfile HM3_b37merge --cluster --mind .05 --mds-plot 10 \
  --extract local.snplist.txt --noweb --out HM3_b37mds
# making the MDS plot
awk 'BEGIN{OFS=","};{print $1, $2, $3, $4, $5, $6, $7}' HM3_b37mds.mds >> HM3_b37mds2R.mds.csv
```

<!-- Then, I made the plot locally:

```R
library(calibrate)
mds.cluster = read.csv("~/data/tmp/HM3_b37mds2R.mds.csv", header=T);
colors=rep("red",length(mds.cluster$C1));
colors[which(mds.cluster$FID == "CEU")] <- "lightblue";
colors[which(mds.cluster$FID == "CHB")] <- "brown";
colors[which(mds.cluster$FID == "YRI")] <- "yellow";
colors[which(mds.cluster$FID == "TSI")] <- "green";
colors[which(mds.cluster$FID == "JPT")] <- "purple";
colors[which(mds.cluster$FID == "CHD")] <- "orange";
colors[which(mds.cluster$FID == "MEX")] <- "grey50";
colors[which(mds.cluster$FID == "GIH")] <- "black";
colors[which(mds.cluster$FID == "ASW")] <- "darkolivegreen";
colors[which(mds.cluster$FID == "LWK")] <- "magenta";
colors[which(mds.cluster$FID == "MKK")] <- "darkblue";
# pdf(file="mdsplot.pdf",width=7,height=7)
plot(rev(mds.cluster$C2), rev(mds.cluster$C1), col=rev(colors),
         ylab="Dimension 1", xlab="Dimension 2",pch=20)
legend("topright", c("My Sample", "CEU", "CHB", "YRI", "TSI", "JPT", "CHD",
                     "MEX", "GIH", "ASW","LWK", "MKK"),
       fill=c("red", "lightblue", "brown", "yellow", "green", "purple",
              "orange", "grey50", "black", "darkolivegreen", "magenta",
              "darkblue"))
# if you want to know the subject ID label of each sample on the graph,
# uncomment the value below
# FIDlabels <- c("CEU", "CHB", "YRI", "TSI", "JPT", "CHD", "MEX", "GIH", "ASW",
#                "LWK", "MKK");
# textxy(mds.cluster[which(!(mds.cluster$FID %in% FIDlabels)), "C2"],
#        mds.cluster[which(!(mds.cluster$FID %in% FIDlabels)), "C1"],
#        mds.cluster[which(!(mds.cluster$FID %in% FIDlabels)), "IID"])
# dev.off();
```

![](images/2019-12-06-18-06-07.png)

Now, for the imputation:

```bash
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep ambig | awk '{print $1}' > ambig.list
plink --bfile $datafile --exclude ambig.list --make-founders --out lastQC \
  --maf 0.01 --hwe 0.000001 --make-bed --noweb
awk '{print $2, $1":"$4}' lastQC.bim > updateSNPs.txt
plink --bfile lastQC --update-name updateSNPs.txt --make-bed --out lastQCb37 \
  --noweb --list-duplicate-vars
plink --bfile lastQCb37 --exclude lastQCb37.dupvar --out lastQCb37_noduplicates \
  --make-bed --noweb
module load vcftools
for i in {1..22}; do
  plink --bfile lastQCb37_noduplicates --chr $i --recode vcf --out NCR_chr"$i";
  vcf-sort NCR_chr"$i".vcf | bgzip -c > NCR_chr"$i".vcf.gz
done
```

Then, uploading to Michigan Imputation Server.

![](images/2019-12-06-17-57-52.png)


But I got an error from the imputation server:

```
Warning: 1 Chunk(s) excluded: < 3 SNPs (see chunks-excluded.txt for details).
Warning: 145 Chunk(s) excluded: at least one sample has a call rate < 50.0% (see chunks-excluded.txt for details).
Remaining chunk(s): 8
Error: More than 100 obvious strand flips have been detected. Please check strand. Imputation cannot be started!
```

# 2019-12-09 09:06:13

I certainly need to understand the strand issue a bit better. But since
imputation might take a while, for now I'll follow the steps I've taken before
to fix the strand issue, and try that:

```bash
cd /data/NCR_SBRB/NCR_genetics/v2
module load shapeit
refdir=/fdb/impute2/1000Genomes_Phase3_integrated_haplotypes_Oct2014/1000GP_Phase3/
for c in {1..22}; do
    shapeit -check -T 16 -V NCR_chr${c}.vcf.gz \
        --input-ref $refdir/1000GP_Phase3_chr${c}.hap.gz \
        $refdir/1000GP_Phase3_chr${c}.legend.gz $refdir/1000GP_Phase3.sample \
        --output-log chr${c}.alignments;
done
# format the files:
for c in {1..22}; do
    grep Strand chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> flip_snps.txt;
    grep Missing chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> missing_snps.txt;
done
# I'm having some issues with duplicate ID that even the list-duplicate command
# in ENIGMA's protocol is not finding, because they have different alleles. So,
# let's remove them completely from the analysis, before we flip it using ShapeIt results:
module load plink
plink --bfile lastQCb37 --write-snplist --out all_snps
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink --bfile lastQCb37 --exclude duplicated_snps.snplist --make-bed --out lastQCb37_noduplicates
# flip and remove all bad ids
plink --bfile lastQCb37_noduplicates --flip flip_snps.txt \
    --exclude missing_snps.txt --make-bed --out lastQCb37_noduplicates_flipped
#reconstruct the VCFs as above to send it to the imputation server.
module load vcftools
for i in {1..22}; do
    plink --bfile lastQCb37_noduplicates_flipped --chr ${i} --recode-vcf \
        --out NCR_chr${i}_flipped;
    vcf-sort NCR_chr"$i"_flipped.vcf | bgzip -c > NCR_chr"$i"_flipped.vcf.gz
done
```

And then we try the imputation server again.

![](images/2019-12-09-09-51-09.png)

Even though it now survives the server QC, only 8 chromossomes are working
ebcause we have a few samples with low call rate. Let's see if we can remove
that in PLINK or if we have to go back to GenomeStudio.

```bash
export datafileraw=merged_inter_noCtrl_sexClean_noDups
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb \
      --make-bed --out ${datafileraw}_filtered
export datafile=${datafileraw}_filtered
export datafile=merged_inter_noCtrl_sexClean_noDups_filtered
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep ambig | awk '{print $1}' > ambig.list
plink --bfile $datafile --exclude ambig.list --make-founders --out lastQC \
    --maf 0.01 --hwe 0.000001 --mind .05 --make-bed --noweb
awk '{print $2, $1":"$4}' lastQC.bim > updateSNPs.txt
plink --bfile lastQC --update-name updateSNPs.txt --make-bed --out lastQCb37 \
    --noweb --list-duplicate-vars -->

 -->
