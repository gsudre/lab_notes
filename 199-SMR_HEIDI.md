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