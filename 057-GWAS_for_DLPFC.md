# 2019-11-29 08:52:14

Philip said we are running GWAS for out cohort and then for ABCD, not PGC. That
makes life easier, although I still want to read more and understand all those
things about allele representation. For now, just to get the ball rolling, let's
run GWAS with our cohort and send out the code to Louk. Then I can work with
ABCD data, because it will take some work. Specifically, I'll need to download
the newer release, as detailed in Argyris' student e-mail from 10/7/2019:

```
Hi Argyris, Hi Dr. Shaw,

Our collaborator, Dr. Ric Anney, at Cardiff did the qc on the release 2 data and found sex mismatches between SNP sex and sex from the demographic file. The mismatches seem to contain in some plates (i.e., batches). 
 
However, the ABCD seems to fix it in the release 2.0.1. So far, we found only a few mismatches (see below), and these mismatches are not specific to certain plates. See attached for the qc result and log file. The ABCD still suggests we remove one plate (plate 461). Note also that the sex information in the plink file in the release 2 is wrong--they removed that information in the release 2.0.1. 
 
SNPSEX 1= MALE 2= FEMALE 3=OUTSIDE THRESHOLDS
 
                   |              SNPSEX
       gender_demo |         0          1          2 |     Total
-------------------+---------------------------------+----------
                 F |       183          9      4,871 |     5,063
                 M |         2      5,541         15 |     5,558
-------------------+---------------------------------+----------
             Total |       185      5,550      4,886 |    10,621
 
Note: Some of the anomalies are pure errors male = female (very small (n= 9,15), and feasibly explained aware as identifying other that sex at birth) others are thresholds issues (n=185).
```

So, let's get started. I'm going to re-create our merged file using only the
intersection SNPs:

```bash
cd /data/NCR_SBRB/NCR_genetics

merge_intersection() {
    plink --bfile $1 --write-snplist
    plink --bfile $2 --extract plink.snplist --make-bed --out ${2}_filtered
    plink --bfile $2 --write-snplist
    plink --bfile $1 --extract plink.snplist --make-bed --out ${1}_filtered
    plink --bfile $1_filtered -bmerge ${2}_filtered.bed ${2}_filtered.bim \
        ${2}_filtered.fam --make-bed --noweb --out $3
}
merge_intersection Shaw01_2017_1to22 Shaw02_2017_1to22 mergeShaw2017_inter1
merge_intersection mergeShaw2017_inter1 Shaw03_2017_1to22 mergeShaw2017_inter2
merge_intersection mergeShaw2017_inter2 Shaw04_2017_1to22 mergeShaw2017_inter3
merge_intersection mergeShaw2017_inter3 Shaw05_2017_1to22 mergeShaw2017_inter4
merge_intersection mergeShaw2017_inter4 Shaw06_2017_1to22 mergeShaw2017_inter5
merge_intersection mergeShaw2017_inter5 Shaw07_2017_1to22 mergeShaw2017_inter6
merge_intersection Shaw04_1to22 Shaw03_1to22 mergeShaw_inter
merge_intersection mergeShaw_inter mergeShaw2017_inter6 merge1_inter
merge_intersection merge1_inter twins_1to22 merge2_inter
merge_intersection merge2_inter CIDR_1to22 merged_inter
```

