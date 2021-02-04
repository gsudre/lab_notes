# 2021-02-04 11:45:42

Let's take a stab at running MESC:

https://github.com/douglasyao/mesc/wiki/Estimating-overall-expression-scores

First, let's make the expression matrix. Maybe I can just re-read Derek's
original file, remove the outliers, and add whatever info is still remaining.

The one caveta here is that our RNAseq is in GRCh38, and everything else (plink
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
               columns=c('GENEID','TXCHROM','EXONSTART'),
               "GENEID")
# store gene names in geneCounts without version in end of name
tx_meta = data.frame(GENEID = substr(rownames(count_matrix), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
count_matrix = count_matrix[imautosome, ]
tx_meta = tx_meta[imautosome, ]
```

# TODO
 * might need to lift over PLINK files to GRCh38 because gene counts are there: https://github.com/sritchie73/liftOverPlink
 * add covariates, especially batch