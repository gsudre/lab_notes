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
1s are cases! Otherwise it thinks affected are 2s, and we wouldn't have any...

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
    fwrite(m, file=out_fname, row.names=F, quote=F, compress='gzip', sep=' ', na='NA')
}
```

# 2020-01-02 10:52:58

I should have all the results by now. So, let's compile them and do some quick
check so see what survives. Just running the code above in 4 different tmuxs.

```bash
cd /data/NCR_SBRB/NCR_genetics/v2/1KG/GWAS
# checking how many hits under .05 we had in chr1
zcat assoc_chr1_scaled_inatt.gz | awk '{ if ($9 < .05) { print $_ } }' - | wc -l
```

I can also get started with the ABCD GWAS, which will follow the exact same
format as the one above. The main issue is that I cannot combine the different
chromossomes in PLINK, because the process is being killed. Let's see how much
of that I can do with the broken down files...

```bash
cd /data/NCR_SBRB/ABCD/v201/1KG/GWAS
```

Note that I had to merge famids to the phenotype files Sam sent so PLINK would
see them.

```r
library(gdata)
a = read.xls('KEY_ab_subject_mergeID_09102019.xlsx')
phen = read.csv('abcd_pheno_1132019.csv')
m = merge(phen, a, by='mergeID', all.x=F, all.y=F)
m$FID = 0
m$IID = sapply(1:nrow(m),
               function(x) paste0(strsplit(as.character(m$ab_key[x]),
                                                        '_')[[1]][1:3],
                                  collapse='_'))
m2 = m[, c(27, 28, 1:26)]
# we only need to update sex for whom we have phenotypes
sex = m2[, c('FID', 'IID', 'Sex')]
sex$binSex = as.numeric(sex$Sex)-1
write.table(m2, file='abcd_pheno_1132019_withFAMIDs.txt', row.names=F, quote=F)
write.table(sex, file='update_sex.txt', row.names=F, col.names=F, quote=F)
```

```bash
cd /data/NCR_SBRB/ABCD/v201/1KG/GWAS
for chr in {1..22}; do
    plink --bfile ../chr${chr} --update-sex update_sex.txt --make-bed --out chr${chr}_sex;
done
```

Now that all sex variables are updated, let's run the association per
chromosome:

```bash
cd /data/NCR_SBRB/ABCD/v201/1KG/GWAS
phen=scaled_inatt;
for chr in {1..22}; do
    plink --bfile chr${chr}_sex --linear sex hide-covar \
        --pheno abcd_pheno_1132019_withFAMIDs.txt --pheno-name $phen \
        --covar abcd_pheno_1132019_withFAMIDs.txt \
        --covar-name scaled_age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr${chr}_${phen};
done
```

This is taking a long time, so I might have to swarm it and run locally... we'll
see. Let me generate the updated sex files above first, and then write the
instructions to Louk. That should give it some time to get some associations
run.

## Instructions for Louk

There are many ways you can get the same results in PLINK. This is what I did
for our dataset.

I started from PLINK binary files with all our data after imputation using the
Michigan Server. We get dosage files from the imputation server, so I imported
them with:

```bash
cd /data/NCR_SBRB/NCR_genetics/v2/1KG
for c in {1..22}; do
    plink --vcf chr${c}.dose.vcf.gz --make-bed --out chr${c};
done
rm -rf merge_list.txt;
for c in {2..22}; do
    echo "chr${c}" >> merge_list.txt;
done
plink --bfile chr1 --merge-list merge_list.txt  --make-bed --out NCR_1KG
```

Then, make sure that the sex assignment in the .fam file is correct, as it will
be one of our covariates. Below I specified to remove anyone who does not have
sex assigned in the .fam file. If you need to update it, the .txt file should
have 3 columns: FAMID, IID, sex. I also had to change the IDs in my .fam to
match what was in my phenotype.

In your case, you can download the phenotypes Sam created from:

```
https://hpc.nih.gov/~sudregp/genr_pheno_1132019.csv
```

and see what needs to be done.

The phenotype file will need to have the first two columns as FID and IID,
matching the same ones in the .fam (not in the same order, though). I added the
extra column (FID) in R:

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

So, this step will be very site specific. Here's what I did on my end:

```bash
mkdir GWAS
cd GWAS
# convert from NSB to the project's mergeID (N+MRN)
plink --bfile ../NCR_1KG --update-ids to_mergeids.txt --make-bed --out NCR_1KG_ids
# update the sex in .fam:
plink --bfile NCR_1KG_ids --update-sex update_sex.txt --make-bed --out NCR_1KG_ids_sex
```

Finally, I decided to run the GWAS by chromosomes, because it's easier to
modularize it that way. Especially in the ABCD dataset, this becomes quite
important. There are 3 continuous phenotypes and one binary. So, in our dataset, I did:

```bash
phen=scaled_inatt;  # or scaled_hi, scaled_ADHD
for chr in {1..22}; do
    plink --bfile NCR_1KG_ids_sex --chr $chr --linear sex hide-covar \
        --pheno nhgri_pheno_1132019_withFAMIDs.txt --pheno-name $phen \
        --covar nhgri_pheno_1132019_withFAMIDs.txt \
        --covar-name scaled_age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr${chr}_${phen};
done
```

For --logistic, I'll also need to add the --1 flag so that 0s are unnafected and
1s are cases! Otherwise it thinks affected are 2s, and we wouldn't have any...

```bash
phen=dx_six_or_more;
for chr in {1..22}; do
    plink --bfile NCR_1KG_ids_sex --chr $chr --1 --logistic sex hide-covar \
        --pheno nhgri_pheno_1132019_withFAMIDs.txt --pheno-name $phen \
        --covar nhgri_pheno_1132019_withFAMIDs.txt \
        --covar-name scaled_age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr${chr}_${phen};
done
```

Once you have results, I merged them with the INFO columns from the imputation,
as Veera suggested in his e-mail, so that we can easily filter what SNPs we want
to play with later. I did it in R:

```r
library(data.table)
for (phen in c('scaled_inatt', 'scaled_hi', 'scaled_ADHD', 'dx_six_or_more')) {
    if (phen == 'dx_six_or_more') {
        suf = 'logistic'
    } else {
        suf = 'linear'
    }
    for (chr in 1:22) {
        print(sprintf('%s, chr %d', phen, chr))
        res = fread(sprintf('assoc_chr%d_%s.assoc.%s', chr, phen, suf),
                    header = T, sep = ' ')
        info = fread(sprintf('../chr%d.info.gz', chr), header = T, sep = '\t')
        m = merge(res, info, by='SNP', all.x=T, all.y=F)
        out_fname = sprintf('assoc_chr%d_%s.gz', chr, phen)
        fwrite(m, file=out_fname, row.names=F, quote=F, compress='gzip',
               sep=' ', na='NA')
    }
}
```

Now the results are just big compressed text tables. It's fast to do operations
in the comand line, like counting the number of variables at some random
significance:

```bash
zcat assoc_chr1_scaled_inatt.gz | awk '{ if ($9 < .05) { print $_ } }' - | wc -l
```

The ABCD data was too big to be combined into one big file, so I kept the per-chromosome binary
files in PLINK. Not much else changes in the analysis, except for the calls for
the actual GWAS command that don't take the chr number anymore, but a changing
binary file name.

## Back to GWAS in ABCD

So, my tests of running it locally didn't help much. I might as well just do one
tmux per chromosome and go from there. I can load about 10 at once per big
machine, if not more... I'm actually running all 22 at once without any issues.
So, let's just do that

# 2020-01-03 09:11:03

Chatting with Philip, I should try running the WNH only, and then I saw Veera's
e-mail about using the different batches and arrays as covariates. Since I'll be
doing this a few times, let's set up parallel.

I should likely also remove the text results from the analysis, as they're
already incorporated in the .gz. It's a 20x decrease in size...

```bash
cd /data/NCR_SBRB/NCR_genetics/v2/1KG/GWAS;
head -n 1 nhgri_pheno_1132019_withFAMIDs.txt > nhgri_pheno_1132019_withFAMIDs_WNHonly.txt;
grep WNH nhgri_pheno_1132019_withFAMIDs.txt >> nhgri_pheno_1132019_withFAMIDs_WNHonly.txt;
phen=scaled_inatt;  # or scaled_hi, scaled_ADHD

for i in {1..22}; do echo $i; done | parallel --max-args=1 \
    plink --bfile NCR_1KG_ids_sex --chr {} --linear sex hide-covar \
        --pheno nhgri_pheno_1132019_withFAMIDs_WNHonly.txt --pheno-name $phen \
        --covar nhgri_pheno_1132019_withFAMIDs_WNHonly.txt \
        --covar-name scaled_age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr{}_${phen}_WNHonly;
```

And to merge the results we do:

```r
library(data.table)
for (phen in c('scaled_inatt', 'scaled_hi', 'scaled_ADHD', 'dx_six_or_more')) {
    if (phen == 'dx_six_or_more') {
        suf = 'logistic'
    } else {
        suf = 'linear'
    }
    for (chr in 1:22) {
        print(sprintf('%s, chr %d', phen, chr))
        res = fread(sprintf('assoc_chr%d_%s_WNHonly.assoc.%s', chr, phen, suf),
                    header = T, sep = ' ')
        info = fread(sprintf('../chr%d.info.gz', chr), header = T, sep = '\t')
        m = merge(res, info, by='SNP', all.x=T, all.y=F)
        out_fname = sprintf('assoc_chr%d_%s_WNHonly.gz', chr, phen)
        fwrite(m, file=out_fname, row.names=F, quote=F, compress='gzip',
               sep=' ', na='NA')
    }
}
```

And we can do the same thing for ABCD:

```bash
cd /data/NCR_SBRB/ABCD/v201/1KG/GWAS;
head -n 1 abcd_pheno_1132019_withFAMIDs.txt > abcd_pheno_1132019_withFAMIDs_WNHonly.txt;
grep WNH abcd_pheno_1132019_withFAMIDs.txt >> abcd_pheno_1132019_withFAMIDs_WNHonly.txt;

phen=scaled_inatt;  # or scaled_hi, scaled_ADHD
for i in {1..22}; do echo $i; done | parallel --max-args=1 \
    plink --bfile chr{}_sex --linear sex hide-covar \
        --pheno abcd_pheno_1132019_withFAMIDs_WNHonly.txt --pheno-name $phen \
        --covar abcd_pheno_1132019_withFAMIDs_WNHonly.txt \
        --covar-name scaled_age, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr{}_${phen}_WNHonly;
```

According to Veera we should always keep the PCs, even if running WNHonly
analysis. So, to save time, I'm not even going to run the noPCcovariates version
for ABCD.

Sam made new files for ABCD including the plate and batch as covariates. There
is one including the problematic plate 461 and one without it. There are only 56
samples in that plate, so I won't use them. If the analysis breaks because of 56
samples it's wea to begin with. I'll use plate and batch as covariates as well.

```r
library(gdata)
a = read.xls('KEY_ab_subject_mergeID_09102019.xlsx')
phen = read.csv('abcd_pheno_01022020.csv')
m = merge(phen, a, by='mergeID', all.x=F, all.y=F)
m$FID = 0
m$IID = sapply(1:nrow(m),
               function(x) paste0(strsplit(as.character(m$ab_key[x]),
                                                        '_')[[1]][1:3],
                                  collapse='_'))
m2 = m[, c(29, 30, 1:28)]
write.table(m2, file='abcd_pheno_01022020_withFAMIDs.txt', row.names=F, quote=F)

# had to do it here because there's one ABCD ID with WNH in the name
idx = m2$PC_DEFINED_SUBGROUPS=='WNH'
write.table(m2[idx,], file='abcd_pheno_01022020_withFAMIDs_WNHonly.txt',
            row.names=F, quote=F)
```

```bash
cd /data/NCR_SBRB/ABCD/v201/1KG/GWAS;

phen=scaled_inatt;  # or scaled_hi, scaled_ADHD
for i in {1..22}; do echo $i; done | parallel --max-args=1 \
    plink --bfile chr{}_sex --linear sex hide-covar \
        --pheno abcd_pheno_01022020_withFAMIDs_WNHonly.txt --pheno-name $phen \
        --covar abcd_pheno_01022020_withFAMIDs_WNHonly.txt \
        --covar-name scaled_age, Plate, batch_number, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr{}_${phen}_WNHonly_agePlateBatchPCs;

phen=dx_six_or_more;
for i in {1..22}; do echo $i; done | parallel --max-args=1 \
    plink --bfile chr{}_sex --1 --logistic sex hide-covar \
        --pheno abcd_pheno_01022020_withFAMIDs_WNHonly.txt --pheno-name $phen \
        --covar abcd_pheno_01022020_withFAMIDs_WNHonly.txt \
        --covar-name scaled_age, Plate, batch_number, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr{}_${phen}_WNHonly_agePlateBatchPCs;
```

And Sam added the batch numbers to our cohort. For our data, they all have the
same chip, like ABCD. Some have different versions of the chip, but not like
PNC that have actual different chips. We do have different batches, though.
There are five batches, so just number them sequentially: twins, CIDR,
Shaw03+Shaw04, Shaw0?_2017, and Shaw0?_2019.

```r
fam = read.table('NCR_1KG_ids.fam')
colnames(fam)[1:2] = c('FID', 'mergeID')
phen = read.csv('nhgri_pheno_01032020.csv')
phen2 = merge(fam[, 1:2], phen, by='mergeID', all.x=F, all.y=T)
phen3 = phen2[, c(2, 1, 3:26)]
colnames(phen3)[2] = 'IID'
write.table(phen3, file='nhgri_pheno_01032020_withFAMIDs.txt', row.names=F, quote=F)
idx = phen3$PC_DEFINED_SUBGROUPS=='WNH'
write.table(phen3[idx,], file='nhgri_pheno_01032020_withFAMIDs_WNHonly.txt',
            row.names=F, quote=F)
```

```bash
cd /data/NCR_SBRB/NCR_genetics/v2/1KG/GWAS;

# our dataset runs quite fast, so no need to parallel
phen=scaled_inatt;  # or scaled_hi, scaled_ADHD
for i in {1..22}; do
    plink --bfile NCR_1KG_ids_sex --chr ${i} --linear sex hide-covar \
        --pheno nhgri_pheno_01032020_withFAMIDs_WNHonly.txt --pheno-name $phen \
        --covar nhgri_pheno_01032020_withFAMIDs_WNHonly.txt \
        --covar-name scaled_age, batch_number, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr${i}_${phen}_WNHonly_ageBatchPCs;
done

phen=dx_six_or_more;
for i in {1..22}; do
    plink --bfile NCR_1KG_ids_sex --chr ${i} --1 --logistic sex hide-covar \
        --pheno nhgri_pheno_01032020_withFAMIDs_WNHonly.txt --pheno-name $phen \
        --covar nhgri_pheno_01032020_withFAMIDs_WNHonly.txt \
        --covar-name scaled_age, batch_number, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10 \
        --out assoc_chr${i}_${phen}_WNHonly_ageBatchPCs;
done
```

```r
library(data.table)
setwd('/data/NCR_SBRB/NCR_genetics/v2/1KG/GWAS')
for (phen in c('scaled_inatt', 'scaled_hi', 'scaled_ADHD', 'dx_six_or_more')) {
    if (phen == 'dx_six_or_more') {
        suf = 'logistic'
    } else {
        suf = 'linear'
    }
    for (chr in 1:22) {
        print(sprintf('%s, chr %d', phen, chr))
        res = fread(sprintf('assoc_chr%d_%s_WNHonly_ageBatchPCs.assoc.%s', chr, phen, suf),
                    header = T, sep = ' ')
        info = fread(sprintf('../chr%d.info.gz', chr), header = T, sep = '\t')
        m = merge(res, info, by='SNP', all.x=T, all.y=F)
        out_fname = sprintf('assoc_chr%d_%s_WNHonly_ageBatchPCs.gz', chr, phen)
        fwrite(m, file=out_fname, row.names=F, quote=F, compress='gzip',
               sep=' ', na='NA')
    }
}
```

```r
library(data.table)
setwd('/data/NCR_SBRB/ABCD/v201/1KG/GWAS/')
for (phen in c('scaled_inatt', 'scaled_hi', 'scaled_ADHD', 'dx_six_or_more')) {
    if (phen == 'dx_six_or_more') {
        suf = 'logistic'
    } else {
        suf = 'linear'
    }
    for (chr in 1:22) {
        print(sprintf('%s, chr %d', phen, chr))
        res = fread(sprintf('assoc_chr%d_%s_WNHonly_agePlateBatchPCs.assoc.%s', chr, phen, suf),
                    header = T, sep = ' ')
        info = fread(sprintf('../chr%d.info.gz', chr), header = T, sep = '\t')
        m = merge(res, info, by='SNP', all.x=T, all.y=F)
        out_fname = sprintf('assoc_chr%d_%s_WNHonly_agePlateBatchPCs.gz', chr, phen)
        fwrite(m, file=out_fname, row.names=F, quote=F, compress='gzip',
               sep=' ', na='NA')
    }
}
```

# TODO
