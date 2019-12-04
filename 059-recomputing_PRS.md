# 2019-12-04 15:00:37

Let's recompute PRS now that we have a few more arrays. We can also transcribe
the Evernote note here.

Note that I can potentially still do some further cleaning inside GenomeStudio.
For example, following this paper: 

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6171493/pdf/bbx012.pdf

But since I cannot find the raw data anymore (dealing with IT to retrieve
backups because it looks like it was moved accidentally), I'll go from PLINK
files from now. At least those look like they are still in the cluster. When I
get the raw data back I'll check if this will make any difference.

I'll do many of the same checks in PLINK anyways (like call rate, and MAF), but
I might be able to do better after fixing some clusters, for example.

```bash
cd /data/NCR_SBRB/NCR_genetics
mkdir v2
cd v2
module load plink
# converting to PLINK binary format
plink --file ~/data/prs/PLINK_160517_0102/Shaw03_2017 --out ./Shaw03_2017
plink --file ~/data/prs/PLINK_160517_0102/Shaw03_2017 --make-bed --out ./Shaw03_2017
plink --file ~/data/prs/PLINK_160517_0143/Shaw04_2017 --make-bed --out ./Shaw04_2017
plink --file ~/data/prs/PLINK_160517_0338/Shaw07_2017 --make-bed --out ./Shaw07_2017
plink --file ~/data/prs/PLINK_160517_0412/Shaw06_2017 --make-bed --out ./Shaw06_2017
plink --file ~/data/prs/PLINK_160517_0533/Shaw05_2017 --make-bed --out ./Shaw05_2017
plink --file ~/data/prs/PLINK_160517_0711/twins --make-bed --out ./twins
plink --file ~/data/prs/PLINK_160517_0852/Shaw03 --make-bed --out ./Shaw03
plink --file ~/data/prs/PLINK_160517_0953/Shaw04 --make-bed --out ./Shaw04
plink --file ~/data/prs/PLINK_160517_1052/Shaw01_2017 --make-bed --out ./Shaw01_2017
plink --file ~/data/prs/PLINK_160517_1143/Shaw02_2017 --make-bed --out ./Shaw02_2017
plink --file ~/data/prs/PLINK_220517_0838/Shaw_Final_Release_05112015 --make-bed --out ./CIDR
plink --file PLINK_231019_0742/Shaw01_2019 --make-bed --out ./Shaw01_2019
plink --file PLINK_231019_0821/Shaw02_2019 --make-bed --out ./Shaw02_2019
plink --file PLINK_231019_0905/Shaw03_2019 --make-bed --out ./Shaw03_2019
```

Maybe we can reduce this to only the SNPs present in all boxes?

```bash
merge_intersection() {
    plink --bfile $1 --write-snplist
    plink --bfile $2 --extract plink.snplist --make-bed --out ${2}_filtered
    plink --bfile $2 --write-snplist
    plink --bfile $1 --extract plink.snplist --make-bed --out ${1}_filtered
    plink --bfile $1_filtered -bmerge ${2}_filtered.bed ${2}_filtered.bim \
        ${2}_filtered.fam --make-bed --noweb --out $3
}
merge_intersection Shaw01_2017 Shaw02_2017 mergeShaw2017_inter1
merge_intersection mergeShaw2017_inter1 Shaw03_2017 mergeShaw2017_inter2
merge_intersection mergeShaw2017_inter2 Shaw04_2017 mergeShaw2017_inter3
merge_intersection mergeShaw2017_inter3 Shaw05_2017 mergeShaw2017_inter4
merge_intersection mergeShaw2017_inter4 Shaw06_2017 mergeShaw2017_inter5
merge_intersection mergeShaw2017_inter5 Shaw07_2017 mergeShaw2017_inter6

merge_intersection Shaw01_2019 Shaw02_2019 mergeShaw2019_inter1
merge_intersection mergeShaw2019_inter1 Shaw03_2019 mergeShaw2019_inter2

merge_intersection Shaw04 Shaw03 mergeShaw_inter

merge_intersection mergeShaw2019_inter2 mergeShaw2017_inter6 merge1_inter
merge_intersection merge1_inter mergeShaw_inter merge2_inter
merge_intersection merge2_inter twins merge3_inter
merge_intersection merge3_inter CIDR merged_inter
```















I decided to re-run everything in geno3:

awk '{print $2}' ~/pgc2017/adhd_jun2017 > pgc_snps.txt
module load plink
plink --file ~/data/prs/PLINK_160517_0102/Shaw03_2017 --out Shaw03_2017
plink --file ~/data/prs/PLINK_160517_0102/Shaw03_2017 --make-bed --out Shaw03_2017
plink --file ~/data/prs/PLINK_160517_0143/Shaw04_2017 --make-bed --out Shaw04_2017
plink --file ~/data/prs/PLINK_160517_0338/Shaw07_2017 --make-bed --out Shaw07_2017
plink --file ~/data/prs/PLINK_160517_0412/Shaw06_2017 --make-bed --out Shaw06_2017
plink --file ~/data/prs/PLINK_160517_0533/Shaw05_2017 --make-bed --out Shaw05_2017
plink --file ~/data/prs/PLINK_160517_0711/twins --make-bed --out twins
plink --file ~/data/prs/PLINK_160517_0852/Shaw03 --make-bed --out Shaw03
plink --file ~/data/prs/PLINK_160517_0953/Shaw04 --make-bed --out Shaw04
plink --file ~/data/prs/PLINK_160517_1052/Shaw01_2017 --make-bed --out Shaw01_2017
plink --file ~/data/prs/PLINK_160517_1143/Shaw02_2017 --make-bed --out Shaw02_2017
plink --file ~/data/prs/PLINK_220517_0838/Shaw_Final_Release_05112015 --make-bed --out CIDR

ls -1 *bim | sed -e 's/\.bim//g' > boxes.txt
rm failed_sexcheck_ids;

# make family and sex update regardless of family
n=`cat update_sex.txt | wc -l`;
rm inds_ids.txt
for i in `seq 1 $n`; do echo 1 >> inds_ids.txt; done;
cut -f 2,3 update_sex.txt > tmp.txt;
paste inds_ids.txt tmp.txt > noid_update_sex.txt
n=`cat update_ids.txt | wc -l`;
rm inds_ids.txt
for i in `seq 1 $n`; do echo 1 >> inds_ids.txt; done;
cut -f 2,3,4 update_ids.txt > tmp.txt;
paste inds_ids.txt tmp.txt > noid_update_ids.txt

while read box; do
     # remove family info
     n=`cat ${box}.fam | wc -l`;
     rm inds_ids.txt;
     for i in `seq 1 $n`; do echo 1 >> inds_ids.txt; done;
     cut -d' ' -f 2-6 ${box}.fam > fam_no_famids.txt;
     cp ${box}.fam ${box}.old_fam;
     paste -d ' ' inds_ids.txt fam_no_famids.txt > ${box}.fam

     plink --bfile ${box} --update-sex noid_update_sex.txt --make-bed --out $box;
     plink --bfile ${box} --update-ids noid_update_ids.txt --make-bed --out $box;
     plink --bfile ${box} --check-sex --out ${box};   
     grep "PROBLEM" ${box}.sexcheck | awk '{print $1, $2}' >> failed_sexcheck_ids;
     plink --bfile ${box} -chr 1-22 --extract intersect_snps.txt --make-bed --out ${box}_1to22;
     plink --bfile ${box}_1to22 --extract pgc_snps.txt --make-bed --out ${box}_1to22_PGC;
done < boxes.txt

plink --bfile twins_1to22_PGC --merge-list merge_list.txt --make-bed --out merged

Then added the following to failed_sex_ids to be removed as well, due to double sampling or twins:

10278 10230
871 223
871 332
8 1116
813 1253
10489 10202
10490 10248
10296 8679
413 580
471 623
528 678
813 1254
10491 2635
10495   10098
1036    10095
1036    WG1023567-DNAC01-10094@0175461097
911     10566  
911     10603
1036    WG1023567-DNAC01-10094@0175461097
1036     10095
10306   10067
10306     10313
10429   2704
10429     2705
10483   10288  
10483     10294
10484   10300  
10484     10306
10485   2746   
10485     2747
10486   10270  
10486     10276
10487   10261  
10487     10264
10488   2354   
10488     2355
10489   WG1023563-DNAD12-2332@0175460578       
10498     10206
10490   WG1023563-DNAC09-2220@0175460566       
10490     10244
10491   WG1023567-DNAC02-2635@0175460548       
10491     10342
10492   1247   
10492     1264
10498   864    
10498     865
10499   883    
10499     884
10501   969    
10501 970


box=merged
plink --bfile ${box} --remove failed_sexcheck_ids --make-bed --out ${box}_noDups
plink --bfile ${box}_noDups --maf .01 --geno .01 --hwe .00001 --mind .01 --make-bed --out ${box}_noDups_clean

Using KING to double check we got rid of twins:

/data/NCR_SBRB/software/KING/king -b merged_noDups_clean.bed

And now we regenerate PRS:

R --file=${HOME}/PRSice_v1.25/PRSice_v1.25.R -q --args  \
  plink  ${HOME}/PRSice_v1.25/plink_1.9_linux_160914 \
  base ~/pgc2017/adhd_jun2017  \
  target ../merged_noDups_clean \
  report.individual.scores T \
  wd ./noFlip \
  cleanup F \
  report.best.score.only F \
  covary F \
  fastscore T \
  barchart.levels 1e-5,1e-4,1e-3,1e-2,1e-1,5e-5,5e-4,5e-3,5e-2,2e-1,3e-1,4e-1,5e-1

And also for the EUR GWAS:

R --file=${HOME}/PRSice_v1.25/PRSice_v1.25.R -q --args  \
  plink  ${HOME}/PRSice_v1.25/plink_1.9_linux_160914 \
  base ~/pgc2017/adhd_eur_jun2017  \
  target ../merged_noDups_clean \
  report.individual.scores T \
  wd ./noFlip_eur \
  cleanup F \
  report.best.score.only F \
  covary F \
  fastscore T \
  barchart.levels 1e-5,1e-4,1e-3,1e-2,1e-1,5e-5,5e-4,5e-3,5e-2,2e-1,3e-1,4e-1,5e-1

Now, we work on imputations:

for i in {1..22}; do   plink --bfile merged_noDups_clean --chr ${i} --recode-vcf --out vcf_chr${i}; done
module load vcftools
for f in `ls vcf*.vcf`; do echo $f; vcf-sort $f | bgzip -c > ${f}.gz; done

And we should do a sneak peak using shapeit to see what kind of flips and removals we'll need to do:

module load shapeit
for c in {1..22}; do shapeit -check -T 16 -V vcf_chr${c}.vcf.gz --input-ref ../1000GP_Phase3/1000GP_Phase3_chr${c}.hap.gz ../1000GP_Phase3/1000GP_Phase3_chr${c}.legend.gz ../1000GP_Phase3/1000GP_Phase3.sample --output-log chr${c}.alignments; done

I downloaded the reference panel from https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html

Now we need to properly flip and remove all bad ids. First, format the files:

for c in {1..22}; do grep Strand chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> flip_snps.txt; done
for c in {1..22}; do grep Missing chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> missing_snps.txt; done

plink --bfile merged_noDups_clean --flip flip_snps.txt --exclude missing_snps.txt --make-bed --out merged_noDups_clean_flipped

Let's just make sure the set is still clean:

plink --bfile merged_noDups_clean_flipped --missing --out ibd
awk '$6 > .01 {print $0}' ibd.imiss | wc -l
awk '$5 > .01 {print $0}' ibd.lmiss | wc -l

And construct the VCFs as above to send it to the imputation server.

for i in {1..22}; do   plink --bfile merged_noDups_clean_flipped --chr ${i} --recode-vcf --out vcf_chr${i}_flipped; done
for f in `ls vcf*flipped.vcf`; do echo $f; vcf-sort $f | bgzip -c > ${f}.gz; done

We also do the population analysis, first in the non-imputed data as we wait on the imputation server:

/data/NCR_SBRB/software/KING/king -b merged_noDups_clean.bed --mds

KING 2.1 - (c) 2010-2018 Wei-Min Chen

The following parameters are in effect:
                   Binary File : merged_noDups_clean.bed (-bname)

Additional Options
         Close Relative Inference : --related, --duplicate
   Pairwise Relatedness Inference : --kinship, --ibdseg, --ibs, --homo
              Inference Parameter : --degree
         Relationship Application : --unrelated, --cluster, --build
                        QC Report : --bysample, --bySNP, --roh, --autoQC
                     QC Parameter : --callrateN, --callrateM
             Population Structure : --pca, --mds [ON]
              Structure Parameter : --projection
              Disease Association : --tdt
   Quantitative Trait Association : --mtscore
                Association Model : --trait [], --covariate []
            Association Parameter : --invnorm, --maxP
               Genetic Risk Score : --risk, --model [], --prevalence
              Computing Parameter : --cpus
                           Output : --prefix [king], --rplot

KING starts at Wed Oct 18 23:52:26 2017
Loading genotype data in PLINK binary format...
Read in PLINK fam file merged_noDups_clean.fam...
  PLINK pedigrees loaded: 996 samples
Read in PLINK bim file merged_noDups_clean.bim...
  Genotype data consist of 544897 autosome SNPs
  PLINK maps loaded: 544897 SNPs
Read in PLINK bed file merged_noDups_clean.bed...
  PLINK binary genotypes loaded.
  129 MB memory allocated for KING format genotype data.
  1 CPU cores are used to convert data from SNP-major to individual-major...
    KING format genotype data successfully converted.
MDS starts at Wed Oct 18 23:52:30 2017
Genotypes stored in 8515 words for each of 996 individuals.
1 CPU cores are used.
SVD starts at Wed Oct 18 23:53:01 2017
done
Largest 20 eigenvalues: 10.43 2.76 1.72 1.52 1.49 1.30 1.26 1.25 1.23 1.21 1.18 1.16 1.14 1.11 1.10 1.07 1.03 1.02 1.00 0.99
MDS ends at Wed Oct 18 23:53:14 2017
20 principal components saved in files kingpc.dat and kingpc.ped

And copied them to Jen.
