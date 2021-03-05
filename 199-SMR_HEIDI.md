# 2021-03-01 15:09:26

After reading a bit more in the FUSION webpage, I found two different software
that could help with the mediation analysis. COLOC and SMR/HEIDI. I'll go for
the latter first, because thoeir papers looked a bit more like what we want to
do:

 * https://www.nature.com/articles/s41467-018-04558-1
 * https://www.nature.com/articles/s41467-018-03371-0
 * https://www.nature.com/articles/ng.3538
 * https://www.um.edu.mt/__data/assets/pdf_file/0005/289427/eQTL_intro.pdf

The first step is to create the eQTL file. I'll go with MatrixeQTL, because it
looks fast and easy to use. But anything in that last link could do. And
SMR/HEIDI (https://cnsgenomics.com/software/smr/#DataManagement) can easily
convert from MatrixeQTL format to their own.

FastQTL is also an option (http://fastqtl.sourceforge.net), as it can read the
dosage VCFs directly. I might need to rename things though, so maybe easier just
to use MatrixEQTL?

```bash
cd ~/data/post_mortem/genotyping/1KG
module load plink
plink --bfile PM_1KG_genop05MAFbtp01rsbtp9_renamed_hbcc \
    --export A --out genotyped_out
```

Then, in R:

```r
library(data.table)
dread = fread('~/data/post_mortem/genotyping/1KG/genotyped_out.raw',
              header = T, sep = ' ')
d = as.data.frame(dread)
d = d[, c(2, 7:ncol(d))]
d2 = t(d)
saveRDS(d2, file='~/data/post_mortem/genotyping/1KG/genotype_raw.rds')
```

# 2021-03-02 11:34:36

So, since MatrixEQTL is in R, it makes sense to not write out the matrices and
just do it all in R:

```r
library("MatrixEQTL")

# SNPs
d = readRDS('~/data/post_mortem/genotyping/1KG/genotype_raw.rds')
colnames(d) = d[1, ]
d = d[2:nrow(d), ]
# d = readRDS('~/data/post_mortem/genotyping/genotype_raw.rds')

# this file has results for all gene subtypes
load('~/data/post_mortem/DGE_03022021.RData')
res = dge_acc[['all']]
covars = t(res$design[, 3:ncol(res$design)])
library(DESeq2)
clean_count = counts(res$dds)
colnames(clean_count) = colnames(covars)

# making sure all subjects are present
d = d[, colnames(d) %in% colnames(covars)]
covars = covars[, colnames(covars) %in% colnames(d)]
clean_count = clean_count[, colnames(clean_count) %in% colnames(d)]

snps = MatrixEQTL::SlicedData(d)
cvrt = MatrixEQTL::SlicedData(covars)

# computing eQTLs
useModel = modelLINEAR;
errorCovariance = numeric();
pvOutputThreshold = 1e-2;

# couldn't house everything in memory
cnt = 1
step = 99
while (cnt < nrow(clean_count)) {
    cat(cnt, '\n')
    mymax = min(cnt + step, nrow(clean_count))
    
    gene = MatrixEQTL::SlicedData(clean_count[cnt:mymax,])
    output_file_name = sprintf('/scratch/sudregp/matrixeqtl_%d.txt', cnt);
    cnt = cnt + step + 1

    me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = TRUE,
    noFDRsaveMemory = TRUE);
}
```

That's breaking... lack of memory.

Let me use the non-imputed data to see if it helps.

```bash
cd ~/data/post_mortem/genotyping/Strandscript-master/test
module load plink
plink --file flipped_PM_test --export A --out ../../genotyped_out
```

```r
library(data.table)
dread = fread('~/data/post_mortem/genotyping/genotyped_out.raw',
              header = T, sep = ' ')
d = as.data.frame(dread)
d = d[, c(2, 7:ncol(d))]
d = t(d)
colnames(d) = d[1, ]
d2 = d[2:nrow(d), ]
a = sapply(colnames(d2), function(x) { br = strsplit(x, '_')[[1]][1];
                                       as.numeric(gsub(br, pattern='BR',
                                                      replacement=''))})
colnames(d2) = a
d3 = matrix(as.numeric(d2), nrow=nrow(d2), ncol=ncol(d2))
rownames(d3) = rownames(d2)
colnames(d3) = colnames(d2)
saveRDS(d3, file='~/data/post_mortem/genotyping/genotype_raw.rds')
```

And let's try the same code as above, but using this non-imputed genotype. Still
getting killed... Let me see if I can get a bigger machine, otherwise I'll need
to somewhow parallelize it.

Now we combine all files and construct the initial BESD:

```bash
#bw
module load SMR
cd ~/data/tmp
head -n +1 matrixeqtl_1.txt > mateQTL.txt;
for f in `/bin/ls matrixeqtl_*txt`; do
    echo $f;
    tail -n +2 $f >> mateQTL.txt;
done
smr --eqtl-summary mateQTL.txt --matrix-eqtl-format --make-besd --out mybesd
```

```
[sudregp@cn3511 tmp]$ smr --eqtl-summary mateQTL.txt --matrix-eqtl-format --make-besd --out mybesd
*******************************************************************
* Summary-data-based Mendelian Randomization (SMR)
* version 1.03
* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang
* The University of Queensland
* MIT License
*******************************************************************
Analysis started: 20:19:1,Tue Mar 2,2021

Options:
--eqtl-summary mateQTL.txt
--matrix-eqtl-format 
--make-besd 
--out mybesd

Reading eQTL summary data from mateQTL.txt ...
Total 5 columns in the Matrix eQTL output.
189024009 rows to be included from mateQTL.txt.

Generating the .epi file...
25943 probes have been saved in the file mybesd.epi.

Generating the .esi file...
646900 SNPs have been saved in the file mybesd.esi.

Generating the .besd file...
Effect sizes (beta) and SE for 25943 probes and 646900 SNPs have been saved in a binary file [mybesd.besd].

Analysis completed: 20:37:21,Tue Mar 2,2021
Computational time: 0:18:20
```

# 2021-03-03 19:18:07

Since we're using the non-imputed data, we're still in grch38, for both genotype
and expression data. Maybe we can just lift over the GWAS? Let's plan on that.

```r
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"),
               columns=c('GENEID','TXCHROM', 'TXSTART', 'TXEND'),
               "GENEID")
txdf = txdf[!duplicated(txdf$GENEID),] 

epi = read.table('~/data/tmp/mybesd.epi', header=0)
epi$GENEID = substr(epi$V2, 1, 15)
tx_meta = merge(epi, txdf, by='GENEID', sort=F)
epi2 = tx_meta[, c('TXCHROM', 'V2', 'V3', 'TXSTART', 'GENEID')]
epi2$orient = '+'
write.table(epi2, row.names=F, col.names=F, quote=F,
            file='~/data/tmp/mybesd.epi.new')
```

Looking at the esi file, the rsids have _Allele in their names. When googling
for those rsids, the allele seems to be the alt allele... not always.

This is getting very confusing, and I think it's because of the whole Ilumina
standard. Let's then do this using the imputed data. It'll take a bit longer,
but we can do it for one and then script it out. I'll need to use the liftOver
expression as well, but for now let's just make sure all script work all the way
to the end, just using the first 100 transcripts. FYI, 100 transcripts is taking
about 33Gb of memory with the imputed genotype, and finished in 375sec.

I think I figured it out.

```bash
cd /scratch/sudregp
module load plink
plink \
    --file ~/data/post_mortem/genotyping/Strandscript-master/test/flipped_PM_test
    --make-bed --out tmp
```

Now my .bim has the correct format. Let's go back to the original plan:

```bash
cd ~/data/tmp
cp /scratch/sudregp/tmp.bim mybesd.esi.new
cp mybesd.esi mybesd.esi.old
cp mybesd.esi.new mybesd.esi
cp mybesd.epi mybesd.epi.old
cp mybesd.epi.new mybesd.epi
```

Now let's keep on going with the data SMRS needs. First we need a LD file. I'll
calculate it from 1KG Phase 3, which I have already converted to binary PLINK in
note 188. Let's then just concatenate them all:

```bash
cd /scratch/sudregp
mydir=~/data/post_mortem/genotyping/1KG
for i in {2..22}; do
    echo "$mydir/chr$i.bed $mydir/chr$i.bim $mydir/chr$i.fam" >> merge_list.txt;
done
plink --bfile $mydir/chr1 --merge-list merge_list.txt --make-bed --out all_1KG
```

That crashed. Let's use the reduce dataset, which might be enough:

```bash
cd /scratch/sudregp
plink --vcf /fdb/1000genomes/release/20130502/reduced.ALL.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz --make-bed --out ALL_reduced
module load SMR
smr --bfile ALL_reduced --make-bld --r --ld-wind 4000 --out ~/data/tmp/mybld
```

Not sure if I need to liftOver my GWAS file. Do rsids change between reference
genomes? The file SMRS takes as input doesn't take in positions, simply the rsid
name. I'm waiting on further information about MAF, so I might need to wait ot
run this anyways. Or, I can see if it runs without MAF?

# 2021-03-04 07:33:48

This took the whole night and still nothing. Let me try to do it per chromosome
then.

```bash
cd /scratch/sudregp
module load SMR
smr --bfile ~/data/post_mortem/genotyping/1KG/ch1_1KG --make-bld --r --ld-wind 4000 --chr 1 --out chr1
```

Again, this is taking forever to save the file. And it's taking almost 100Gb of
memory with a single chromossome... 

Maybe I could just subset the entire analysis by the GWAS SNPs? The GWAS has 8M
SNPs. Just chr1 in the ALL set has 6.5M, and the reduced version of the entire
dataset has 47M. So, maybe I could extract the SNPs I want from GWAS in each chr
from the ALL set, and then merge them all together. Then, I can use the bfile
directly in SMR...

```bash
cd /scratch/sudregp
awk '{ print $2 }' ~/pgc2017/adhd_eur_jun2017 | tail -n +2 > keep_snps.txt;
for c in {1..22}; do
    plink --bfile ~/data/post_mortem/genotyping/1KG/ch${c}_1KG \
        --extract keep_snps.txt --make-bed --out ch${c}_1KG_gwasOnly;
done
for c in {1..22}; do
    plink --bfile ch${c}_1KG_gwasOnly --list-duplicate-vars;
    plink --bfile ch${c}_1KG_gwasOnly --exclude plink.dupvar --make-bed \
        -out ch${c}_1KG_gwasOnly_noDups;
done
rm -f merge_list.txt;
for i in {2..22}; do
    echo "ch${i}_1KG_gwasOnly_noDups.bed ch${i}_1KG_gwasOnly_noDups.bim ch${i}_1KG_gwasOnly_noDups.fam" >> merge_list.txt;
done
plink --bfile ch1_1KG_gwasOnly_noDups --merge-list merge_list.txt \
    --make-bed --out all_1KG_gwasOnly_noDups

```

This is not working, because many times the same variant have two different
positions. So, we need to keep only the good ones in the VCF after filtering by
chromosome and position in the VCF. Just a simple update-name won't work because
PLINK won't know which one to update. So, we'll need to recreate the .bed using
chr:pos instead of rsid.

```bash
awk '{ print $1":"$3":"$4":"$5 }' ~/pgc2017/adhd_eur_jun2017 | tail -n +2  > keep_ids_chrpos.txt;
cat keep_ids_chrpos.txt | sort -u -k 1,1 | uniq > unique_keep_ids_chrpos.txt;
# this takes a long time, so best to parallelize it
module load bcftools
for c in {1..22}; do
    echo "bcftools annotate -Ob -x 'ID' -I +'%CHROM:%POS:%REF:%ALT' /fdb/1000genomes/release/20130502/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > chr${c}_ren.bcf && bcftools index chr${c}_ren.bcf && plink --bcf chr${c}_ren.bcf --allow-extra-chr --extract unique_keep_ids_chrpos.txt --make-bed --out chr${c}_byPos;" >> my_cmds.txt
done
cat my_cmds.txt | parallel -j 22 --max-args=1 {};

rm -f merge_list.txt;
for c in {2..22}; do
    echo "chr${c}_byPos.bed chr${c}_byPos.bim chr${c}_byPos.fam" >> merge_list.txt;
done
plink --bfile chr${c}_byPos --merge-list merge_list.txt \
    --make-bed --out all_1KG_byPos
```

This might work, and it might even be a better way to go about it in the future
because it's less variable than simply the rsids. For now, let's try a
reconversion keeping only the bilallelic ones:

```bash
for i in {4..22}; do echo $i >> junk.txt; done
cat junk.txt | parallel -j 18 plink --vcf /fdb/1000genomes/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --make-bed --biallelic-only strict list --extract keep_snps.txt \
        --out ch{}_1KG_biAllelicOnly_gwasOnly;

for c in {1..22}; do
    plink --bfile ch${c}_1KG_biAllelicOnly_gwasOnly --list-duplicate-vars;
    # wiping out any remaining IDs... might be removing too much here
    tail -n +2 plink.dupvar | awk '{ print $4 }' > rm_ids.txt;
    tail -n +2 plink.dupvar | awk '{ print $5 }' >> rm_ids.txt;
    plink --bfile ch${c}_1KG_biAllelicOnly_gwasOnly --write-snplist --out all_snps;
    cat all_snps.snplist | sort | uniq -d >> rm_ids.txt;
    plink --bfile ch${c}_1KG_biAllelicOnly_gwasOnly --exclude rm_ids.txt \
        --make-bed -out ch${c}_1KG_biAllelicOnly_gwasOnly_noDups;
done
rm -f merge_list.txt;
for i in {2..22}; do
    echo "ch${i}_1KG_biAllelicOnly_gwasOnly_noDups.bed ch${i}_1KG_biAllelicOnly_gwasOnly_noDups.bim ch${i}_1KG_biAllelicOnly_gwasOnly_noDups.fam" >> merge_list.txt;
done
plink --bfile ch1_1KG_biAllelicOnly_gwasOnly_noDups --merge-list merge_list.txt \
    --make-bed --out all_1KG_biAllelicOnly_gwasOnly_noDups
```

At this stage I have all_1KG_biAllelicOnly_gwasOnly_noDups and all_1KG_byPos as
two possible candidates to run SMRS on. They are already filtered to GWAS SNPs
and I can just calculate LD on the fly. I'll still need to change back to rsids
in the second file because of the way the eQTL was run. But I could also try
lifting those position and re-renaming the eQTL file? We'll see. Let's start
with the first file, which might be less of a headache. The biAllelic version
has more than twice as many snps as byPos though.

Let's see how this runs, and we can double check all inputs later:

```bash
module load SMR
echo "SNP    A1  A2  freq    b   se  p   n" > mygwas.ma;
tail -n +2 ~/pgc2017/adhd_eur_jun2017 | \
    awk ' { print $2, $4, $5, "NA", $6, $7, $8, 53293 } ' >> mygwas.ma;
smr --bfile /scratch/sudregp/all_1KG_biAllelicOnly_gwasOnly_noDups \
    --gwas-summary mygwas.ma --beqtl-summary mybesd --out mysmr --thread-num 10
```

Great. This ran all the way to the end. So, we finally have a pipeline. Now,
it's just making sure all is are dotted. I can also try using the allel
frequency from 1KG, under the assumption I won't get it from PGC.

```bash
cd /scratch/sudregp/
plink --bfile all_1KG_biAllelicOnly_gwasOnly_noDups --freq
```

```r
library(data.table)
dread = fread('~/pgc2017/adhd_eur_jun2017', header = T, sep = '\t')
gwas = as.data.frame(dread)
dread = fread('/scratch/sudregp/plink.frq', header = T, sep = ' ')
freqs = as.data.frame(dread)
m = merge(gwas, freqs, by='SNP')
m2 = m[, c('SNP', 'A1.x', 'A2.x', 'MAF', 'OR', 'SE', 'P')]
m2$n = 53293
colnames(m2) = c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'n')
write.table(m2, row.names=F, quote=F, file='~/data/tmp/mygwas.ma')
```

That actually broke because of the allele frequency discrepancy in the datasets.
Maybe we won't do it that way then.

OK, so let's redo the entire analysis.

## A new beginning...

Let's focus only in the GWAS SNPs. For that, let's make sure that, when we
recode them, everything still looks good:

```bash
cd /scratch/sudregp/
gwas=~/pgc2017/adhd_eur_jun2017
tail -n +2 $gwas | awk '{ print $1":"$3 }' > gwas_renamed_snps.txt;
# I won't specify the alleles because the GWAS file doesn't mean minor or major
# only reference

# keep only the GWAS SNPs in the imputed PM data using the converted and 
# filtered data prior to renaming using the HRC file
plink_file=~/data/post_mortem/genotyping/1KG/PM_1KG_genop05MAFbtp01rsbtp9
awk '{ print $1 ":" $4 }' ${plink_file}.bim > new_name.txt;
awk '{ print $2 }' ${plink_file}.bim > old_name.txt;
paste old_name.txt new_name.txt > update_snps.txt
plink --bfile $plink_file --update-name update_snps.txt --make-bed --out tmp
plink --bfile tmp -extract gwas_renamed_snps.txt --make-bed \
    --out PM_1KG_gwasOnly
plink --bfile PM_1KG_gwasOnly --list-duplicate-vars;
# wiping out any remaining IDs... might be removing too much here
tail -n +2 plink.dupvar | awk '{ print $4 }' > rm_ids.txt;
tail -n +2 plink.dupvar | awk '{ print $5 }' >> rm_ids.txt;
plink --bfile PM_1KG_gwasOnly --write-snplist --out all_snps;
cat all_snps.snplist | sort | uniq -d >> rm_ids.txt;
plink --bfile PM_1KG_gwasOnly --exclude rm_ids.txt --make-bed \
    -out PM_1KG_gwasOnly_noDups;
plink --bfile PM_1KG_gwasOnly_noDups --export A --out genotyped_out
```

Now we need to run the EQTL, but only for the transcripts that remained when
lifted over to hg19:

```r
library(data.table)
dread = fread('/scratch/sudregp/genotyped_out.raw',
              header = T, sep = ' ')
d = as.data.frame(dread)
d = d[, c(2, 7:ncol(d))]
d2 = t(d)
colnames(d2) = d2[1, ]
d2 = d2[2:nrow(d2), ]
a = sapply(colnames(d2), function(x) { br = strsplit(x, '_')[[1]][2];
                                       as.numeric(gsub(br, pattern='BR',
                                                      replacement=''))})
colnames(d2) = a
d3 = matrix(as.numeric(d2), nrow=nrow(d2), ncol=ncol(d2))
rownames(d3) = rownames(d2)
colnames(d3) = colnames(d2)
saveRDS(d3, file='/scratch/sudregp/genotype_raw.rds')
```

Now that the genotype data is a bit cleaner, let's match it to the rest of the
data we'll need:

```r
# SNPs
d = readRDS('/scratch/sudregp/genotype_raw.rds')

# this file has results for all gene subtypes
load('~/data/post_mortem/DGE_03022021.RData')
res = dge_acc[['all']]
covars = t(res$design[, 3:ncol(res$design)])
library(DESeq2)
clean_count = counts(res$dds)
colnames(clean_count) = colnames(covars)

# making sure all subjects are present
d = d[, colnames(d) %in% colnames(covars)]
covars = covars[, colnames(covars) %in% colnames(d)]
clean_count = clean_count[, colnames(clean_count) %in% colnames(d)]

# trimming to only lifted transcripts
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"),
               columns=c('GENEID','TXCHROM', 'TXSTART', 'TXEND'),
               "GENEID")
txdf = txdf[!duplicated(txdf$GENEID),] 
tx_meta = data.frame(GENEID = substr(rownames(clean_count), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
tx_meta$bed = sapply(1:nrow(tx_meta),
                     function(x) sprintf('chr%s:%d-%d', tx_meta[x, 'TXCHROM'],
                                         tx_meta[x, 'TXSTART'],
                                         tx_meta[x, 'TXEND']))
write.table(tx_meta$bed, file='/scratch/sudregp/accCleanAllBED_hg38.txt',
            row.names=F, quote=F, col.names=F)

# liftOver at https://genome.ucsc.edu/cgi-bin/hgLiftOver
errs = read.table('/scratch/sudregp/hglft_genome_67f17_1670a0.err.txt')
rm_me = tx_meta$bed %in% errs[,1]
tx_meta_clean = tx_meta[!rm_me, ]
count_matrix_clean = clean_count[!rm_me,]
new_pos = read.table('/scratch/sudregp/hglft_genome_67f17_1670a0.bed')
new_pos = gsub(x=new_pos[, 1], pattern='chr', replacement='')
tx_meta_clean$CHR = sapply(new_pos, function(x) strsplit(x, ':')[[1]][1])
nochr = sapply(new_pos, function(x) strsplit(x, ':')[[1]][2])
tx_meta_clean$GENE_COORD = sapply(nochr, function(x) strsplit(x, '-')[[1]][1])
write.table(tx_meta_clean, row.names=F, quote=F,
            file='/scratch/sudregp/lifted_transcripts_accAllClean.txt')
save(count_matrix_clean, tx_meta_clean, d, covars,
     file='/scratch/sudregp/accAllClean_forMatrixEQTL.rData')
```

Let's code it now so that we can parallelize the eQTL computations, because they
take quite a while.

```r
library("MatrixEQTL")
load('/scratch/sudregp/accAllClean_forMatrixEQTL.rData')

snps = MatrixEQTL::SlicedData(d)
cvrt = MatrixEQTL::SlicedData(covars)

# computing eQTLs
useModel = modelLINEAR;
errorCovariance = numeric();
pvOutputThreshold = 1e-2;

# couldn't house everything in memory
args <- commandArgs(trailingOnly = TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])
# from = 1
# to = 250

mymax = min(to, nrow(count_matrix_clean))
gene = MatrixEQTL::SlicedData(count_matrix_clean[from:mymax,])
output_file_name = sprintf('/scratch/sudregp/matrixeqtl_%05dto%05d.txt',
                            from, mymax);
cat(output_file_name, '\n')
me = Matrix_eQTL_engine(snps = snps,
                        gene = gene,
                        cvrt = cvrt,
                        output_file_name = output_file_name,
                        pvOutputThreshold = pvOutputThreshold,
                        useModel = useModel, 
                        errorCovariance = errorCovariance, 
                        verbose = TRUE,
                        pvalue.hist = TRUE,
                        min.pv.by.genesnp = TRUE,
                        noFDRsaveMemory = TRUE);
```

I wrote that as /scratch/sudregp/tmp.R, and then it's just a matter of swarming
it:

```bash
cd /scratch/sudregp/tmp.R
rm swarm.eqtl
cnt=1;
step=250;
while [ $cnt -lt 25818 ]; do
    let to=$cnt+$step-1;
    echo "Rscript --vanilla ./tmp.R $cnt $to" >> swarm.eqtl;
    let cnt=$cnt+$step;
done
swarm -t 1 -g 60 -f swarm.eqtl --job-name eqtl --time 30:00 \
        --logdir trash -m R --gres=lscratch:10 --partition quick,norm;
```

Now we combine all files and construct the initial BESD:

```bash
#bw
module load SMR
cd /scratch/sudregp/
head -n +1 matrixeqtl_00001to00250.txt > mateQTL.txt;
for f in `/bin/ls matrixeqtl_*to*txt`; do
    echo $f;
    tail -n +2 $f >> mateQTL.txt;
done
smr --eqtl-summary mateQTL.txt --matrix-eqtl-format --make-besd --out mybesd
```

While this is running, I'll go ahead and correct the GWAS file:

```bash
module load SMR
cd /scratch/sudregp/
echo "SNP    A1  A2  freq    b   se  p   n" > mygwas.ma;
tail -n +2 ~/pgc2017/adhd_eur_jun2017 | \
    awk ' { print $1":"$3, $4, $5, "NA", $6, $7, $8, 53293 } ' >> mygwas.ma;
```

Note that I'm still waiting on the right MAFs from PGC, but this should do it
for now.

Also in parallel we can construct a PLINK binary file with only the SNPS we care
about:

```bash
cd /scratch/sudregp/
# this takes a long time, so best to parallelize it
module load bcftools
rm -f my_cmds.txt;
for c in {1..22}; do
    echo "bcftools annotate -Ob -x 'ID' -I +'%CHROM:%POS' /fdb/1000genomes/release/20130502/ALL.chr${c}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz > chr${c}_chrPos.bcf && bcftools index chr${c}_chrPos.bcf && plink --bcf chr${c}_chrPos.bcf --allow-extra-chr --extract gwas_renamed_snps.txt --make-bed --out chr${c}_chrPos;" >> my_cmds.txt
done
cat my_cmds.txt | parallel -j 22 --max-args=1 {};

# removing duplicate SNPs
for c in {1..22}; do
    plink --bfile chr${c}_chrPos --list-duplicate-vars;
    # wiping out any remaining IDs... might be removing too much here
    tail -n +2 plink.dupvar | awk '{ print $4 }' > rm_ids.txt;
    tail -n +2 plink.dupvar | awk '{ print $5 }' >> rm_ids.txt;
    plink --bfile chr${c}_chrPos --write-snplist --out all_snps;
    cat all_snps.snplist | sort | uniq -d >> rm_ids.txt;
    plink --bfile chr${c}_chrPos --exclude rm_ids.txt --make-bed \
        -out chr${c}_chrPos_noDups;
done

rm -f merge_list.txt;
for c in {2..22}; do
    echo "chr${c}_chrPos.bed chr${c}_chrPos.bim chr${c}_chrPos.fam" >> merge_list.txt;
done
plink --bfile chr${c}_chrPos --merge-list merge_list.txt \
    --make-bed --out all_1KG_chrPos
```

Getting lots of duplicate errors here too...

It turns out that creating the BESD from the imputed data is taking too long.
Actually, it's getting killed even in the bigger machine... I'll have to go back
to the non-imputed data.

## Back to the non-imputed data

Let's then first liftOver the PM genotypes to hg19:

```bash
module load ucsc
module load python/2.7
python ~/data/tmp/liftOverPlink-master/liftOverPlink.py \
    -m ~/data/post_mortem/genotyping/Strandscript-master/test/flipped_PM_test.map \
    -p ~/data/post_mortem/genotyping/Strandscript-master/test/flipped_PM_test.ped \
    -o PM_hg19 -c ~/data/tmp/hg38ToHg19.over.chain.gz
# create an output file with just the GWAS snps. GWAs is already in hg19
awk '{ print $2 }' ~/pgc2017/adhd_eur_jun2017 | tail -n +2 > keep_snps.txt;
module load plink
plink --file PM_hg19 --extract keep_snps.txt --allow-extra-chr \
    --export A --out PM_hg19_GWASonly
```

Let's compute eQTL now:

```r
library(data.table)
dread = fread('/scratch/sudregp/PM_hg19_GWASonly.raw', header = T, sep = ' ')
d = as.data.frame(dread)
d = d[, c(2, 7:ncol(d))]
d2 = t(d)
colnames(d2) = d2[1, ]
d2 = d2[2:nrow(d2), ]
a = sapply(colnames(d2), function(x) { br = strsplit(x, '_')[[1]][1];
                                       as.numeric(gsub(br, pattern='BR',
                                                      replacement=''))})
colnames(d2) = a
d3 = matrix(as.numeric(d2), nrow=nrow(d2), ncol=ncol(d2))
rownames(d3) = rownames(d2)
colnames(d3) = colnames(d2)
saveRDS(d3, file='/scratch/sudregp/PM_hg19_GWASonly.rds')
```

Now that the genotype data is a bit cleaner, let's match it to the rest of the
data we'll need:

```r
# SNPs
d = readRDS('/scratch/sudregp/PM_hg19_GWASonly.rds')

# this file has results for all gene subtypes
load('~/data/post_mortem/DGE_03022021.RData')
res = dge_acc[['all']]
covars = t(res$design[, 3:ncol(res$design)])
library(DESeq2)
clean_count = counts(res$dds)
colnames(clean_count) = colnames(covars)

# making sure all subjects are present
d = d[, colnames(d) %in% colnames(covars)]
covars = covars[, colnames(covars) %in% colnames(d)]
clean_count = clean_count[, colnames(clean_count) %in% colnames(d)]

# trimming to only lifted transcripts
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"),
               columns=c('GENEID','TXCHROM', 'TXSTART', 'TXEND'),
               "GENEID")
txdf = txdf[!duplicated(txdf$GENEID),] 
tx_meta = data.frame(GENEID = substr(rownames(clean_count), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
tx_meta$bed = sapply(1:nrow(tx_meta),
                     function(x) sprintf('chr%s:%d-%d', tx_meta[x, 'TXCHROM'],
                                         tx_meta[x, 'TXSTART'],
                                         tx_meta[x, 'TXEND']))
write.table(tx_meta$bed, file='/scratch/sudregp/accCleanAllBED_hg38.txt',
            row.names=F, quote=F, col.names=F)

# liftOver at https://genome.ucsc.edu/cgi-bin/hgLiftOver
errs = read.table('/scratch/sudregp/hglft_genome_67f17_1670a0.err.txt')
rm_me = tx_meta$bed %in% errs[,1]
tx_meta_clean = tx_meta[!rm_me, ]
count_matrix_clean = clean_count[!rm_me,]
new_pos = read.table('/scratch/sudregp/hglft_genome_67f17_1670a0.bed')
new_pos = gsub(x=new_pos[, 1], pattern='chr', replacement='')
tx_meta_clean$CHR = sapply(new_pos, function(x) strsplit(x, ':')[[1]][1])
nochr = sapply(new_pos, function(x) strsplit(x, ':')[[1]][2])
tx_meta_clean$GENE_COORD = sapply(nochr, function(x) strsplit(x, '-')[[1]][1])
write.table(tx_meta_clean, row.names=F, quote=F,
            file='/scratch/sudregp/lifted_transcripts_accAllClean.txt')
save(count_matrix_clean, tx_meta_clean, d, covars,
     file='/scratch/sudregp/accAllCleanPM_forMatrixEQTL.rData')
```

tmp.R is very similar, except that now it'll take much faster as we have less
SNPs (still the same number of transcripts, though):

```r
library("MatrixEQTL")
load('/scratch/sudregp/accAllCleanPM_forMatrixEQTL.rData')

snps = MatrixEQTL::SlicedData(d)
cvrt = MatrixEQTL::SlicedData(covars)

# computing eQTLs
useModel = modelLINEAR;
errorCovariance = numeric();
pvOutputThreshold = 1e-2;

# couldn't house everything in memory
args <- commandArgs(trailingOnly = TRUE)
from = as.numeric(args[1])
to = as.numeric(args[2])
# from = 1
# to = 250

mymax = min(to, nrow(count_matrix_clean))
gene = MatrixEQTL::SlicedData(count_matrix_clean[from:mymax,])
output_file_name = sprintf('/scratch/sudregp/matrixeqtlPM_%05dto%05d.txt',
                            from, mymax);
cat(output_file_name, '\n')
me = Matrix_eQTL_engine(snps = snps,
                        gene = gene,
                        cvrt = cvrt,
                        output_file_name = output_file_name,
                        pvOutputThreshold = pvOutputThreshold,
                        useModel = useModel, 
                        errorCovariance = errorCovariance, 
                        verbose = TRUE,
                        pvalue.hist = TRUE,
                        min.pv.by.genesnp = TRUE,
                        noFDRsaveMemory = TRUE);
```

```bash
cd /scratch/sudregp/
rm swarm.eqtl
cnt=1;
step=250;
while [ $cnt -lt 25818 ]; do
    let to=$cnt+$step-1;
    echo "Rscript --vanilla ./tmp.R $cnt $to" >> swarm.eqtl;
    let cnt=$cnt+$step;
done
swarm -t 1 -g 60 -f swarm.eqtl --job-name eqtlPM --time 30:00 \
        --logdir trash -m R --gres=lscratch:10 --partition quick,norm;

module load SMR
cd /scratch/sudregp/
head -n +1 matrixeqtlPM_00001to00250.txt > mateQTLPM.txt;
for f in `/bin/ls matrixeqtlPM_*to*txt`; do
    echo $f;
    tail -n +2 $f >> mateQTLPM.txt;
done
smr --eqtl-summary mateQTLPM.txt --matrix-eqtl-format --make-besd --out mybesdPM
```

```
*******************************************************************
* Summary-data-based Mendelian Randomization (SMR)
* version 1.03
* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang
* The University of Queensland
* MIT License
*******************************************************************
Analysis started: 18:13:42,Fri Mar 5,2021

Options:--eqtl-summary mateQTLPM.txt
--matrix-eqtl-format 
--make-besd 
--out mybesdPM

Reading eQTL summary data from mateQTLPM.txt ...
Total 5 columns in the Matrix eQTL output.
159825577 rows to be included from mateQTLPM.txt.

Generating the .epi file...
25818 probes have been saved in the file mybesdPM.epi.

Generating the .esi file...
553991 SNPs have been saved in the file mybesdPM.esi.

Generating the .besd file...
Effect sizes (beta) and SE for 25818 probes and 553991 SNPs have been saved in a binary file [mybesdPM.besd].

Analysis completed: 18:35:22,Fri Mar 5,2021
Computational time: 0:21:40
```

Now we need to create our binary PLINK file for LD computation, but again we can
extract only the GWAS SNPs to make it slimmer:

```bash
cd /scratch/sudregp/
rm -f junk.txt;
for i in {1..22}; do echo $i >> junk.txt; done
cat junk.txt | parallel -j 22 plink --vcf /fdb/1000genomes/release/20130502/ALL.chr{}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
        --make-bed --biallelic-only strict list --extract keep_snps.txt \
        --out ch{}_1KG_biAllelicOnly_gwasOnly;

# some extra cleaning to remove duplicate SNPs
for c in {1..22}; do
    plink --bfile ch${c}_1KG_biAllelicOnly_gwasOnly --list-duplicate-vars;
    # wiping out any remaining IDs... might be removing too much here
    tail -n +2 plink.dupvar | awk '{ print $4 }' > rm_ids.txt;
    tail -n +2 plink.dupvar | awk '{ print $5 }' >> rm_ids.txt;
    plink --bfile ch${c}_1KG_biAllelicOnly_gwasOnly --write-snplist --out all_snps;
    cat all_snps.snplist | sort | uniq -d >> rm_ids.txt;
    plink --bfile ch${c}_1KG_biAllelicOnly_gwasOnly --exclude rm_ids.txt \
        --make-bed -out ch${c}_1KG_biAllelicOnly_gwasOnly_noDups;
done

rm -f merge_list.txt;
for i in {2..22}; do
    echo "ch${i}_1KG_biAllelicOnly_gwasOnly_noDups.bed ch${i}_1KG_biAllelicOnly_gwasOnly_noDups.bim ch${i}_1KG_biAllelicOnly_gwasOnly_noDups.fam" >> merge_list.txt;
done
plink --bfile ch1_1KG_biAllelicOnly_gwasOnly_noDups \
    --merge-list merge_list.txt \
    --make-bed --out all_1KG_biAllelicOnly_gwasOnly_noDups
```

And let's clean up the esi and epi files:

```r
txdf = read.table('/scratch/sudregp/lifted_transcripts_accAllClean.txt',
                  header=T)
epi = read.table('/scratch/sudregp/mybesdPM.epi', header=0)
epi$GENEID = substr(epi$V2, 1, 15)
tx_meta = merge(epi, txdf, by='GENEID', sort=F)
epi2 = tx_meta[, c('TXCHROM', 'V2', 'V3', 'GENE_COORD', 'GENEID')]
epi2$orient = '+'
write.table(epi2, row.names=F, col.names=F, quote=F,
            file='/scratch/sudregp/mybesdPM.epi.new')
```

And we copy the esi from bim:

```bash
cd /scratch/sudregp
module load plink
plink --file PM_hg19 --extract keep_snps.txt --allow-extra-chr \
    --make-bed --out PM_hg19_GWASonly

cp PM_hg19_GWASonly.bim mybesdPM.esi.new
cp mybesdPM.esi mybesdPM.esi.old
cp mybesdPM.esi.new mybesdPM.esi
cp mybesdPM.epi mybesdPM.epi.old
cp mybesdPM.epi.new mybesdPM.epi
```

And finally fix the GWAS to run the whole thing:

```bash
cd /scratch/sudregp
echo "SNP    A1  A2  freq    b   se  p   n" > mygwas.ma;
tail -n +2 ~/pgc2017/adhd_eur_jun2017 | \
    awk ' { print $2, $4, $5, "NA", $6, $7, $8, 53293 } ' >> mygwas.ma;

module load SMR

smr --bfile all_1KG_biAllelicOnly_gwasOnly_noDups \
    --gwas-summary mygwas.ma --beqtl-summary mybesdPM \
    --out mysmrPM --thread-num 10
```

```
*******************************************************************
* Summary-data-based Mendelian Randomization (SMR)
* version 1.03
* (C) 2015 Futao Zhang, Zhihong Zhu and Jian Yang
* The University of Queensland
* MIT License
*******************************************************************
Analysis started: 19:17:13,Fri Mar 5,2021

Options:
--bfile all_1KG_biAllelicOnly_gwasOnly_noDups
--gwas-summary mygwas.ma
--beqtl-summary mybesdPM
--out mysmrPM
--thread-num 10

Reading GWAS summary data from [mygwas.ma].
WARNING: frequency is 'NA' in one or more rows.
GWAS summary data of 8094094 SNPs to be included from [mygwas.ma].
Reading eQTL SNP information from [mybesdPM.esi].
553991 SNPs to be included from [mybesdPM.esi].
Reading PLINK FAM file from [all_1KG_biAllelicOnly_gwasOnly_noDups.fam].
2504 individuals to be included from [all_1KG_biAllelicOnly_gwasOnly_noDups.fam].
Reading PLINK BIM file from [all_1KG_biAllelicOnly_gwasOnly_noDups.bim].
7622495 SNPs to be included from [all_1KG_biAllelicOnly_gwasOnly_noDups.bim].
Checking the consistency of the alleles of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).
551289 SNPs are included after allele checking. 
Reading PLINK BED file from [all_1KG_biAllelicOnly_gwasOnly_noDups.bed] in SNP-major format ...
Genotype data for 2504 individuals and 551289 SNPs to be included from [all_1KG_biAllelicOnly_gwasOnly_noDups.bed].
Calculating allele frequencies ...
Checking the consistency of allele frequency of each SNP between pairwise data sets (including the GWAS summary data, the eQTL summary data and the LD reference data).
0 SNPs (0.00% <= 5.00%) with allele frequency differences > 0.20 between any pair of the data sets are excluded from the analysis.
Reading eQTL summary data...
Reading eQTL probe information from [mybesdPM.epi].
25818 Probes to be included from [mybesdPM.epi].
Reading eQTL summary data from [mybesdPM.besd].
eQTL summary data of 25818 Probes to be included from [mybesdPM.besd].

Performing SMR analysis (SMR and HEIDI tests) for 10 probes (with at least a cis-eQTL at p < 5.00e-08)... 
For each probe, the analysis will only include SNPs with eQTL p-values < 1.565400e-03,
then exclude SNPs with LD r-squared between top-SNP > 0.90 or < 0.05, and further exclude one of each pair of the remaining SNPs with LD r-squared > 0.90.
Results of 10 probes have been saved in file mysmrPM.smr.

Analysis completed: 19:22:54,Fri Mar 5,2021
Computational time: 0:5:41
```

We can also try the multi-SNP analysis and also trans:

```bash
smr --bfile all_1KG_biAllelicOnly_gwasOnly_noDups \
    --gwas-summary mygwas.ma --beqtl-summary mybesdPM \
    --out mytransPM --thread-num 10 --trans --trans-wind 1000 

smr --bfile all_1KG_biAllelicOnly_gwasOnly_noDups \
    --gwas-summary mygwas.ma --beqtl-summary mybesdPM \
    --out mymultiPM --smr-multi
```

Not much going on for any of the 3 analysis. Could there be something wrong with
how we are specifying the reference SNPs?

# TODO
 * still unclear on where we should have effect and reference alleles
 * look at including methylation data








plink --bfile ~/data/post_mortem/genotyping/1KG/PM_1KG_genop05MAFbtp01rsbtp9_ren \
    -extract kee_ids.txt --make-bed --out PM_1KG_gwasOnly2






plink --bfile lastQCb37 --write-snplist --out all_snps
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink --bfile lastQCb37 --exclude duplicated_snps.snplist --make-bed --out lastQCb37_noduplicates
# flip and remove all bad ids
plink --bfile lastQCb37_noduplicates --flip flip_snps.txt \
    --exclude missing_snps.txt --make-bed --out lastQCb37_noduplicates_flipped










LEt's remove the probes that didn't make it when lifting from grch38 to hg19:



```r
res = run_DGE(count_matrix, data, tx_meta, myregion, NA, .05)
covars = data.frame(FID=rownames(res$design), IID=rownames(res$design))
covars = cbind(covars, res$design[, 3:ncol(res$design)])
write.table(covars, quote=F, row.names=F, sep='\t',
            file='~/tmp/caudate_covars_mesc.txt')

clean_count = counts(res$dds)

library(GenomicFeatures)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
clean_count = clean_count[imautosome, ]
tx_meta = tx_meta[imautosome, ]
tx_meta$bed = sapply(1:nrow(tx_meta),
                     function(x) sprintf('chr%s:%d-%d', tx_meta[x, 'TXCHROM'],
                                         tx_meta[x, 'TXSTART'],
                                         tx_meta[x, 'TXEND']))
write.table(tx_meta$bed, file='~/tmp/tx_meta_caudateClean_hg38.txt',
            row.names=F, quote=F, col.names=F)

errs = read.table('~/tmp/hglft_genome_26baf_8035f0.err.txt')
rm_me = tx_meta$bed %in% errs[,1]
tx_meta_clean = tx_meta[!rm_me, ]
count_matrix_clean = clean_count[!rm_me,]
colnames(count_matrix_clean) = data$hbcc_brain_id
new_pos = read.table('~/tmp/hglft_genome_26baf_8035f0.bed')
new_pos = gsub(x=new_pos[, 1], pattern='chr', replacement='')
tx_meta_clean$CHR = sapply(new_pos, function(x) strsplit(x, ':')[[1]][1])
nochr = sapply(new_pos, function(x) strsplit(x, ':')[[1]][2])
tx_meta_clean$GENE_COORD = sapply(nochr, function(x) strsplit(x, '-')[[1]][1])
out_df = cbind(tx_meta_clean[, c(1, 6, 7)], count_matrix_clean)
colnames(out_df)[1] = 'GENE'
write.table(out_df, quote=F, row.names=F, sep='\t',
            file='~/tmp/caudateClean_mesc.txt')
```