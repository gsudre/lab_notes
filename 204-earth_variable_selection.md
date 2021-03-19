# 2021-03-18 19:42:30

Let's see if we can do our covariate selection similar to:

https://science.sciencemag.org/content/sci/suppl/2018/12/12/362.6420.eaat8127.DC1/aat8127_Gandal_SM.pdf

They first used earth to do variable selection, and then added some SVs to it.
Let's do the earth part first.

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

That was the data loading part, now the actual variable selection:

```r
keep_me = rep(TRUE, nrow(count_matrix))
my_count_matrix = count_matrix[keep_me, ]
my_tx_meta = tx_meta[keep_me, ]
data.pm = data
rownames(data.pm) = data.pm$hbcc_brain_id
colnames(my_count_matrix) = rownames(data.pm)
rownames(my_count_matrix) = substr(rownames(my_count_matrix), 1, 15)

# # removing variables where more than half of the subjects have zero counts
# keep_me = rowSums(my_count_matrix==0) < .25*ncol(my_count_matrix)
# my_count_matrix = my_count_matrix[keep_me, ]
# cat('Keeping', nrow(my_count_matrix), 'after zero removal\n')

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

countdata = round(X)
library(DESeq2)

# from https://biodatascience.github.io/compbio/dist/sva.html
dds <- DESeqDataSetFromMatrix(countData = countdata,
                                    colData = data.pm,
                                    design = ~Diagnosis)
dds <- estimateSizeFactors(dds)
# I get the same value whether I do this after DESeq or just estimateSizeFactors
dat  <- counts(dds, normalized = TRUE)
dat2 = scale(t(dat))
# dat2 = t(dat)
bad_gene = colSums(is.na(dat2)) > 0
dat2 = dat2[, !bad_gene]

# removed pH because of too many NAs, RIN because we have RINe for everyone
num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5')
categ_vars = c('batch', 'MoD', 'substance_group', 'brainbank',
                'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level')
covars = data.pm[, c(num_vars, categ_vars)]
squared = data.pm[, num_vars] ** 2
for (var in num_vars) {
    covars[, var] = scale(data.pm[, var])
    squared[, var] = scale(squared[, var])
}
colnames(squared) = sapply(num_vars, function(x) sprintf('%s_2', x))
covars = cbind(covars, squared)

library(earth)
scores = rep(0, ncol(covars))
names(scores) = colnames(covars)
nperms = 50
set.seed(1234)
for (p in 1:nperms) {
    cat(p, '\n')
    var_idx = sample(1:ncol(dat2), 1, replace=F) 
    a = earth(x=covars, y=dat2[, var_idx])
    print(evimp(a))
    scores[rownames(evimp(a))] = scores[rownames(evimp(a))] + 1
}
```

This is not quite working... it's not picking up anything! It could simply be
because they have 1695 subject, and we only have 60! Maybe something to do later
with our big datasets...

But maybe we could do a different version of it? Say, go sequentially through
all genes, and select the variables that most genes are related to?

```r
library(MASS)
mydat = cbind(dat2, covars)
scores = rep(0, ncol(covars))
names(scores) = colnames(covars)

for (p in 1:ncol(dat2)) {
    if (p %% 100 == 0) {
        cat(p, '\n')
    }
    fm_str = sprintf('%s ~ %s', colnames(dat2)[p], paste0(colnames(covars), collapse='+'))
    res.lm <- lm(as.formula(fm_str), data = mydat)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    selected_vars = as.character(attr(terms(step), 'variables'))
    # remove list and gene name
    selected_vars = selected_vars[3:length(selected_vars)]
    scores[selected_vars] = scores[selected_vars] + 1
}
```




    