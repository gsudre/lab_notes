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

# 2019-12-31 10:15:40

Going back to this analysis, we now have the imputed data for our cohort and
also for ABCD. So, let's run the GWAS. Veera said in his e-mail he'd be using
the genotype data in PNC, but the differences from running that and imputed data
were trivial as long as we have beta and se for all the SNPs. Also he said we
should just run everything, regardless of any filtering, and just make sure to
include the INFO score column in the summary statistics, so that we can exclude
the low info variants during metaxcan. 

For GWAS, the plan was to use Plink (version 1.90), --linear option (if
phenotype quantitative) and –logistic (if phenotype binary). And supply the
covariates using –cov. Covariates are age, gender, PCs (10).

I grabbed the files from Sam in GREX_common/PHENO_SENDING_11132019, and running
scaled_ADHD, scaled_inatt, scaled_hi, and dx_six_or_more columns. And
scaled_age, sex, and PCs for covariates.

```bash
cd /data/NCR_SBRB/NCR_genetics/v2/1KG
mkdir GWAS
cd GWAS
# convert from NSB to the project's mergeID (N+MRN)
plink --bfile ../NCR_1KG --update-ids to_mergeids.txt --make-bed --out NCR_1KG_ids
```

Note that I had to merge famids to the phenotype files Sam sent so PLINK would
see them.

```r
fam = read.table('NCR_1KG_ids.fam')
colnames(fam)[1:2] = c('FID', 'mergeID')
phen = read.csv('nhgri_pheno_1132019.csv')
phen2 = merge(fam[, 1:2], phen, by='mergeID', all.x=F, all.y=T)
phen3 = phen2[, c(2, 1, 3:25)]
colnames(phen3)[2] = 'IID'
# we only need to update sex for whom we have phenotypes
sex = phen3[, c('FID', 'IID', 'Sex')]
sex$binSex = as.numeric(sex$Sex)-1
write.table(phen3, file='nhgri_pheno_1132019_withFAMIDs.txt', row.names=F, quote=F)
write.table(sex, file='update_sex.txt', row.names=F, col.names=F, quote=F)
```

Then, we go back to PLINK:

```bash
plink --bfile NCR_1KG_ids --update-sex update_sex.txt --make-bed --out NCR_1KG_ids_sex
# I double-checked, and only the people that have sex in the .fam are being used when sex is specified below!
# it's also not using any NAs, which is good
plink --bfile NCR_1KG_ids_sex --chr 1 --linear sex hide-covar \
    --pheno nhgri_pheno_1132019_withFAMIDs.txt --pheno-name inatt_9_scale \
    --covar nhgri_pheno_1132019_withFAMIDs.txt \
    --covar-name scaled_age, PC1 \
    --out tmpassoc1
```

So, that's working. Now, I don't think I need to know the associations with all
covariates. Also, the columns I'll need to merge with the results are all in the
info.gz, so that should be straight forward. Let's do it per chromossome first
just to make more manageable files.

Also, the operation only takes one thread, so we can tmux it in a big machine to
make use of all the memory we have.

```bash
phen=scaled_inatt;
for chr in {1..22}; do
    plink --bfile NCR_1KG_ids_sex --chr $chr --linear sex hide-covar \
        --pheno nhgri_pheno_1132019_withFAMIDs.txt --pheno-name $phen \
        --covar nhgri_pheno_1132019_withFAMIDs.txt \
        --covar-name scaled_age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr${chr}_${phen};
done
```

For --logistic, I'll also need to add the --1 flag so that 0s are unnafected and
1s are cases! Otherwise it thinks affected are 2s, and we'd have any...

Then, I imagine we can assume the rows in the same order to do a straight up
merge of the columns? Otherwise we'd have to use R, which will be so slow...

```bash
zcat ../chr1.info.gz | awk '{print $1}' - > snps.info
cat tmpassoc1.assoc.linear | awk '{print $2}' - > snps.res
wc -l snps.res
wc -l snps.info
```

No, they're different :(

```r
library(data.table)
res = fread('tmpassoc1.assoc.linear', header = T, sep = ' ')
info = fread('../chr1.info.gz', header = T, sep = '\t')
m = merge(res, info, by='SNP', all.x=T, all.y=F)
fwrite(m, file='test.gz', row.names=F, quote=F, compress='gzip', sep=' ')
```

But that runs a bit faster than the actual analysis, so we need to wait until
all chromosomes are done before we can gather the results.

```r
library(data.table)
phen = 'scaled_inatt'
for (chr in 1:22) {
    print(chr)
    res = fread(sprintf('assoc_chr%d_%s.assoc.linear', chr, phen),
                header = T, sep = ' ')
    info = fread(sprintf('../chr%d.info.gz', chr), header = T, sep = '\t')
    m = merge(res, info, by='SNP', all.x=T, all.y=F)
    out_fname = sprintf('assoc_chr%d_%s.gz', chr, phen)
    fwrite(m, file=out_fname, row.names=F, quote=F, compress='gzip', sep=' ')
}
```



# TODO
