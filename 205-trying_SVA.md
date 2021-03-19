# 2021-03-19 17:33:06

A quick try with SVA, based on:

http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#using-sva-with-deseq2

and

https://biodatascience.github.io/compbio/dist/sva.html


```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]
# data only contains sample metadata, and count_matrix has actual counts

# cleaning up some variables
data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))
data$Diagnosis = factor(data$Diagnosis, levels=c('Control', 'Case'))
data$substance_group = factor(data$substance_group)
data$comorbid_group = factor(data$comorbid_group_update)
data$evidence_level = factor(data$evidence_level)
data$brainbank = factor(data$bainbank)

# removing everything but autosomes
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
# store gene names in geneCounts without version in end of name
tx_meta = data.frame(GENEID = substr(rownames(count_matrix), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
count_matrix = count_matrix[imautosome, ]
tx_meta = tx_meta[imautosome, ]
```

Now we start the part specific to SVA:

```r
keep_me = rep(TRUE, nrow(count_matrix))
my_count_matrix = count_matrix[keep_me, ]
my_tx_meta = tx_meta[keep_me, ]
data.pm = data
rownames(data.pm) = data.pm$hbcc_brain_id
colnames(my_count_matrix) = rownames(data.pm)
rownames(my_count_matrix) = substr(rownames(my_count_matrix), 1, 15)

# removing variables with zero or near-zero variance
library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(my_count_matrix), method = pp_order)
X = t(predict(pp, t(my_count_matrix)))
cat('Keeping', nrow(X), 'after NZ and NZV filtering\n')

# replace the one subject missing population PCs by the median of their
# self-declared race and ethnicity
idx = (data.pm$Race.x=='White' & data.pm$Ethnicity.x=='Non-Hispanic' &
    !is.na(data.pm$C1))
pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
med_pop = apply(data[idx, pop_pcs], 2, median)
data.pm[which(is.na(data.pm$C1)), pop_pcs] = med_pop

num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5')
categ_vars = c('batch', 'MoD', 'substance_group', 'brainbank',
                'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level',
                'Diagnosis')
covars = data.pm[, c(num_vars, categ_vars)]
for (var in num_vars) {
    covars[, var] = scale(data.pm[, var])
}

countdata = round(X)
library(DESeq2)

# from https://biodatascience.github.io/compbio/dist/sva.html
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                    colData = covars,
                                    design = ~Diagnosis)
dds <- estimateSizeFactors(dds)
# I get the same value whether I do this after DESeq or just estimateSizeFactors
dat  <- counts(dds, normalized = TRUE)

library(sva)
mod  <- model.matrix(~ Diagnosis, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
dat2 <- dat[rowSums(dat) > 0,]
svseq <- svaseq(dat2, mod, mod0, n.sv = 2)
```

So, this works. Let's now try splitting the data, and doing the correlation
analysis.

```r
idx1 = 1:23
idx2 = setdiff(1:nrow(covars), idx1)
dds <- DESeqDataSetFromMatrix(countData = countdata[, idx1],
                                    colData = covars[idx1, ],
                                    design = ~Diagnosis)
dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
mod  <- model.matrix(~ Diagnosis, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
dat2 <- dat[rowSums(dat) > 0,]
svseq <- svaseq(dat2, mod, mod0, n.sv = 2)
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + Diagnosis
ddssva <- DESeq(ddssva)
res1 <- results(ddssva, name = "Diagnosis_Case_vs_Control", alpha = .05)

dds <- DESeqDataSetFromMatrix(countData = countdata[, idx2],
                                    colData = covars[idx2, ],
                                    design = ~Diagnosis)
dds <- estimateSizeFactors(dds)
dat  <- counts(dds, normalized = TRUE)
mod  <- model.matrix(~ Diagnosis, colData(dds))
mod0 <- model.matrix(~   1, colData(dds))
dat2 <- dat[rowSums(dat) > 0,]
svseq <- svaseq(dat2, mod, mod0, n.sv = 2)
ddssva <- dds
ddssva$SV1 <- svseq$sv[,1]
ddssva$SV2 <- svseq$sv[,2]
design(ddssva) <- ~ SV1 + SV2 + Diagnosis
ddssva <- DESeq(ddssva)
res2 <- results(ddssva, name = "Diagnosis_Case_vs_Control", alpha = .05)

both_res = merge(data.frame(res1), data.frame(res2), by=0)
mycorr = cor.test(both_res$log2FoldChange.x, both_res$log2FoldChange.y,
                  method='spearman')
```

OK, so this does the job. Now it's just a matter of scripting it for multiple
permutations and different nmbers of SVs. Also, start with the correct
covariates based on the earth analysis.