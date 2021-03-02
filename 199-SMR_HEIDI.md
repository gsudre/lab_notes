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
# d = readRDS('~/data/post_mortem/genotyping/1KG/genotype_raw.rds')
# colnames(d) = d[1, ]
# d = d[2:nrow(d), ]
d = readRDS('~/data/post_mortem/genotyping/genotype_raw.rds')

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
step = 999
while (cnt < nrow(clean_count)) {
    cat(cnt, '\n')
    mymax = min(cnt + step, nrow(clean_count))
    
    gene = MatrixEQTL::SlicedData(clean_count[cnt:mymax,])
    output_file_name = sprintf('~/data/tmp/matrixeqtl_%d.txt', cnt);
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

Now we combine all files into a .tar.gz and construct the initial BESD:

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