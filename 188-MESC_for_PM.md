# 2021-02-04 11:45:42

Let's take a stab at running MESC:

https://github.com/douglasyao/mesc/wiki/Estimating-overall-expression-scores

First, let's make the expression matrix. Maybe I can just re-read Derek's
original file, remove the outliers, and add whatever info is still remaining.

The one caveat here is that our RNAseq is in GRCh38, and everything else (plink
files, reference 1KG files, and GWAS) is in hg19. It might be easier to just
liftOver the RNAseq to hg19, and do everything there:

```bash
# laptop
conda activate rnaseq
```

```python
cd ~/data/post_mortem/gtfparse
from gtfparse import read_gtf
fn = '/Users/sudregp/data/post_mortem/Homo_sapiens.GRCh38.97.gtf'
df = read_gtf(fn)
out_fn = '/Users/sudregp/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv'
my_cols = ['transcript_id', 'transcript_biotype', 'gene_id', 'gene_biotype']
df[my_cols].to_csv(out_fn)
```

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]
# data only contains sample metadata, and count_matrix has actual counts

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"),
               columns=c('GENEID','TXCHROM', 'TXSTART', 'TXEND'),
               "GENEID")
txdf = txdf[!duplicated(txdf$GENEID),] 
# store gene names in geneCounts without version in end of name
tx_meta = data.frame(GENEID = substr(rownames(count_matrix), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
count_matrix = count_matrix[imautosome, ]
tx_meta = tx_meta[imautosome, ]
```

# 2021-02-23 17:17:16

Let's go back to this analysis. I'll try to liftOver the genes, so I'll need to
get the gene coordinates. They can be start location or midpoints.

```r
tx_meta$bed = sapply(1:nrow(tx_meta),
                     function(x) sprintf('chr%s:%d-%d', tx_meta[x, 'TXCHROM'],
                                         tx_meta[x, 'TXSTART'],
                                         tx_meta[x, 'TXEND']))
write.table(tx_meta$bed, file='~/tmp/tx_meta_hg38.txt', row.names=F, quote=F,
            col.names=F)
```

Then I used the liftOver website to convert to hg19. There were a few errors, so
I'll remove those first, and then add in the converted coordinates.

```r
errs = read.table('~/tmp/hglft_genome_29a45_586cc0.err.txt')
rm_me = tx_meta$bed %in% errs[,1]
tx_meta_clean = tx_meta[!rm_me, ]
cout_matrix_clean = count_matrix[!rm_me,]
new_pos = read.table('~/Downloads/hglft_genome_29a45_586cc0.bed')
new_pos = gsub(x=new_pos[, 1], pattern='chr', replacement='')
tx_meta_clean$CHR = sapply(new_pos, function(x) strsplit(x, ':')[[1]][1])
nochr = sapply(new_pos, function(x) strsplit(x, ':')[[1]][2])
tx_meta_clean$GENE_COORD = sapply(nochr, function(x) strsplit(x, '-')[[1]][1])
```

# TODO
 * might need to lift over PLINK files to GRCh38 because gene counts are there: https://github.com/sritchie73/liftOverPlink
 * add covariates, especially batch
 * use only clean genes?