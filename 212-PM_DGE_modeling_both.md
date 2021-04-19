# 2021-04-19 15:25:45

I was reading the main DESeq2 tutorial, as well as the man page for results(),
and it made me wonder if I cannot get a better result if I model both Caudate
and ACC together.

```r
add_cov=NA


data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name
# this is a ACC outlier
data = data[, ! colnames(data) %in% c('68080')]
# this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# its other replicate
data = data[, ! colnames(data) %in% c('66552')]

library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]

# cleaning up some variables
df$Individual = factor(df$hbcc_brain_id)
df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
df$MoD = factor(df$Manner.of.Death)
df$Sex = factor(df$Sex)
df$batch = factor(df$batch)
df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
df$substance_group = factor(df$substance_group)
df$comorbid_group = factor(df$comorbid_group_update)
df$evidence_level = factor(df$evidence_level)
df$brainbank = factor(df$bainbank)
# replace the one subject missing population PCs by the median of their
# self-declared race and ethnicity
idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
med_pop = apply(df[idx, pop_pcs], 2, median)
df[which(is.na(df$C1)), pop_pcs] = med_pop
                                            
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
            "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                tx_meta$TXCHROM != 'Y' &
                tx_meta$TXCHROM != 'MT')
data = data[imautosome, ]
tx_meta = tx_meta[imautosome, ]

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = data,
                            colData = df,
                            design = ~ Diagnosis)
# don't allow genes with more zeros than the subjects in smallest group
min_subjs = min(table(df$Diagnosis))
keep <- rowSums(counts(dds) == 0) <= min_subjs
dds <- dds[keep,]

dds = DESeq(dds)
norm.cts <- counts(dds, normalized=TRUE)
library(sva)
mm <- model.matrix(~ Region + Diagnosis + Region:Diagnosis, colData(dds))
mm0 <- model.matrix(~ 1, colData(dds))
fit <- svaseq(norm.cts, mod=mm, mod0=mm0, n.sv=2)

ddssva <- dds
fm_str = '~ '
for (s in 1:ncol(fit$sv)) {
    eval(parse(text=sprintf('ddssva$SV%d <- fit$sv[,%d]', s, s)))
    fm_str = sprintf('%s + SV%d', fm_str, s)
}

if (is.na(add_cov)) {
    fm_str = sprintf('%s + Region + Diagnosis + Region:Diagnosis', fm_str)
} else {
    fm_str = sprintf('%s + %s + Region + Diagnosis + Region:Diagnosis',
    fm_str, add_cov)
}
design(ddssva) = as.formula(fm_str)
ddssva = DESeq(ddssva)

library(edgeR)
design = model.matrix(as.formula(fm_str), data=colData(ddssva))
isexpr <- filterByExpr(ddssva, design=design)
ddsexpr = ddssva[isexpr,]

nOutliers = Inf
mydds = ddsexpr
while (nOutliers > 0) {
    cat('Processing', nrow(mydds), 'variables.\n')
    mydds <- DESeq(mydds)
    maxCooks <- apply(assays(mydds)[["cooks"]], 1, max)
    # outlier cut-off uses the 99% quantile of the F(p,m-p) distribution (with 
    # p the number of parameters including the intercept and m number of
    # samples).
    m <- ncol(mydds)
    # number or parameters (SVs + Diagnosis + intercept)
    p <- ncol(fit$sv) + 1 + 1
    if (! is.na(add_cov)) {
        p = p + 1 # add one more to the terms above
    }
    co = qf(.99, p, m - p)
    keep_me = which(maxCooks < co)
    nOutliers = nrow(mydds) - length(keep_me)
    cat('Found', nOutliers, 'outliers.\n')
    mydds = mydds[keep_me, ]
}

resultsNames(mydds)
# the condition effect for Caudate (the main effect, reference level)
results(mydds, contrast=c("Diagnosis","Case","Control"))

# the condition effect for ACC
# this is, by definition, the main effect *plus* the interaction term # (the extra condition effect in genotype II compared to genotype I). 
results(mydds, list( c("Diagnosis_Case_vs_Control", "RegionACC.DiagnosisCase") ))

# the interaction term, answering: is the Case effect *different* across regions?
results(mydds, name="RegionACC.DiagnosisCase")

dds_2sv = mydds

# ...

dds_autosv = mydds
```

I then created the same as above, but with the auto SV, which gave me 17 SVs.
That might be too much, but we'll see how the results look.

