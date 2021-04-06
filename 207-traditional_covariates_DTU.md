# 2021-03-25 13:42:12

Let's see if DTU can helps us decide with which model to go with. The main
difference here is that I could not find any references for using SVA with DTU,
so let's just run the difference between batch and BBB models:

```r
plotProportion <- function(expData = NULL, geneID = NULL, samps = NULL) {
    colnames(expData)[1:2] = c("gid","tid")
    sub = subset(expData, gid == geneID)
    sub = reshape2::melt(sub, id = c("gid", "tid"))
    sub = merge(samps, sub, by.x = "sample_id", by.y = "variable")

    clrs = c("dodgerblue3", "maroon2",  "forestgreen", "darkorange1", "blueviolet", "firebrick2",
"deepskyblue", "orchid2", "chartreuse3", "gold", "slateblue1", "tomato" , "blue", "magenta", "green3",
"yellow", "purple3", "red" ,"darkslategray1", "lightpink1", "lightgreen", "khaki1", "plum3", "salmon")

    p = ggplot(sub, aes(tid, value, color = group, fill = group)) +
    geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8, lwd = 0.5) +
    stat_summary(fun = mean, geom = "point", color = "black", shape = 5, size = 3, position=position_dodge(width = 0.8)) +
    scale_color_manual(values = clrs) + scale_fill_manual(values = clrs) +
    geom_quasirandom(color = "black", size = 1, dodge.width = 0.8) + theme_bw() +
    ggtitle(geneID) + xlab("Transcripts")

    p = p + ylab("Proportions")
    p
}

# posthoc procedure to improve the false discovery rate (FDR) and overall false discovery rate (OFDR) control. It sets the p-values and adjusted p-values for transcripts with small per-sample proportion SD to 1
smallProportionSD <- function(d, filter = 0.1) {
        # Generate count table
        cts = as.matrix(subset(counts(d), select = -c(gene_id, feature_id)))
        # Summarise count total per gene
        gene.cts = rowsum(cts, counts(d)$gene_id)
        # Use total count per gene as count per transcript
        total.cts = gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
        # Calculate proportion of transcript in gene
        props = cts/total.cts
        rownames(props) = rownames(total.cts)
        
        # Calculate standard deviation
        propSD = sqrt(matrixStats::rowVars(props))
        # Check if standard deviation of per-sample proportions is < 0.1
        propSD < filter
}

run_DTU_noPCA = function(count_matrix, samples, tx_meta, subtype, alpha,
                         BBB=FALSE, ncores=NA, add_cov=NA) {
    cat('Starting with', nrow(tx_meta), 'variables\n')
    if (is.na(subtype)) {
        keep_me = rep(TRUE, nrow(count_matrix))
    } else {
        keep_me = grepl(tx_meta$transcript_biotype, pattern=sprintf('%s$', subtype))
    }
    cat('Keeping', sum(keep_me), subtype, 'variables\n')
    my_count_matrix = count_matrix[keep_me, ]
    my_tx_meta = tx_meta[keep_me, ]

    # removing variables where more than half of the subjects have zero counts
    keep_me = rowSums(my_count_matrix==0) < .25*ncol(my_count_matrix)
    my_count_matrix = my_count_matrix[keep_me, ]
    cat('Keeping', nrow(my_count_matrix), 'after zero removal\n')

    # removing variables with zero or near-zero variance
    library(caret)
    pp_order = c('zv', 'nzv')
    pp = preProcess(t(my_count_matrix), method = pp_order)
    X = t(predict(pp, t(my_count_matrix)))
    cat('Keeping', nrow(X), 'after NZ and NZV filtering\n')

    # keep only the remaining transcripts and their corresponding genes
    txdf.sub = my_tx_meta[match(rownames(X), my_tx_meta$transcript_id),]
    counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME)
    counts = cbind(counts, X)

    library(DRIMSeq)
    samples$group = samples$Diagnosis
    samples$sample_id = as.character(samples$submitted_name)
    d0 = dmDSdata(counts = counts, samples = samples)

    n = nrow(samples)
    n.small = min(table(samples$group))

    d = DRIMSeq::dmFilter(d0,
                        min_samps_feature_expr = n.small, min_feature_expr = 10,
                        min_samps_feature_prop = n.small, min_feature_prop = 0.1,
                        min_samps_gene_expr = n, min_gene_expr = 10)

    countData = round(as.matrix(counts(d)[,-c(1:2)]))
    cat('Keeping', nrow(countData), 'after DRIMSeq expression filtering\n')

    data.pm = samples
    # replace the one subject missing population PCs by the median of their
    # self-declared race and ethnicity
    idx = (data.pm$Race.x=='White' & data.pm$Ethnicity.x=='Non-Hispanic' &
           !is.na(data.pm$C1))
    pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
    med_pop = apply(data[idx, pop_pcs], 2, median)
    data.pm[which(is.na(data.pm$C1)), pop_pcs] = med_pop

    if (BBB) {
        data.pm$BBB = factor(sapply(1:nrow(data.pm),
                                    function(x) sprintf('%s_%s',
                                             as.character(data.pm[x,'brainbank']),
                                             as.character(data.pm[x, 'batch']))))
        use_pcs = c('BBB', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI')
    } else {
        use_pcs = c('batch', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI')
    }

    # add more covariates for robustness testing
    if (! is.na(add_cov)) {
        use_pcs = c(use_pcs, add_cov)
    }

    fm_str = sprintf('~ group + %s', paste0(use_pcs, collapse = ' + '))
    cat('Using formula:', fm_str, '\n')

    # scaling num_vars to assure convergence
    # removed pH because of too many NAs, RIN because we have RINe for everyone
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
                'C1', 'C2', 'C3', 'C4', 'C5')
    for (var in num_vars) {
        data.pm[, var] = scale(data.pm[, var])
    }

    design = model.matrix(as.formula(fm_str), data = data.pm)

    set.seed(42)
    system.time({
        if (is.na(ncores)) {
            d <- dmPrecision(d, design = design)
            d <- dmFit(d, design = design)
            d <- dmTest(d, coef = "groupCase")
        } else {
            multicoreParam <- BiocParallel::MulticoreParam(workers = ncores)
            d <- dmPrecision(d, design = design, BPPARAM=multicoreParam)
            d <- dmFit(d, design = design, BPPARAM=multicoreParam)
            d <- dmTest(d, coef = "groupCase", BPPARAM=multicoreParam)
        }
    })
    res.g = DRIMSeq::results(d)
    res.t = DRIMSeq::results(d, level = "feature")

    # make NA pvalues to be 1 so they don't screw up future steps
    no.na <- function(x) ifelse(is.na(x), 1, x)
    res.g$pvalue <- no.na(res.g$pvalue)
    res.t$pvalue <- no.na(res.t$pvalue)

    cat('Genes surviving FDR q <', alpha, ':',
        sum(res.g$adj_pvalue < alpha, na.rm=T), '\n')
    cat('Transcripts surviving FDR q<', alpha, ':',
        sum(res.t$adj_pvalue < alpha, na.rm=T), '\n')

    filt = smallProportionSD(d)

    res.t.filt = DRIMSeq::results(d, level = "feature")
    res.t.filt$pvalue[filt] = 1
    res.t.filt$adj_pvalue[filt] = 1
    res.t.filt$pvalue <- no.na(res.t.filt$pvalue)
    cat('Transcripts removed due to small SD:', sum(filt, na.rm=T), 'out of',
        length(filt), '\n')
    cat('Transcripts surviving SD filtering and FDR q<', alpha, ':',
        sum(res.t.filt$adj_pvalue < alpha, na.rm=T), '\n')

    # strp <- function(x) substr(x, 1, 15)
    # # Construct a vector of per-gene p-values for the screening stage
    # pScreen = res.g$pvalue
    # names(pScreen) = strp(res.g$gene_id)
    # # Construct a one column matrix of the per-transcript confirmation p-values
    # pConfirmation = matrix(res.t.filt$pvalue, ncol = 1)
    # dimnames(pConfirmation) = list(strp(res.t.filt$feature_id), "transcript")
    # # res.t is used twice to construct a 4-column data.frame that contain both original IDs and IDs without version numbers
    # tx2gene = data.frame(res.t[,c("feature_id", "gene_id")], 
    #                     res.t[,c("feature_id", "gene_id")])
    # # remove version from gene name
    # for (i in 1:2) tx2gene[,i] = strp(tx2gene[,i])

    # library(stageR)
    # stageRObj = stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, 
    #                     pScreenAdjusted = FALSE, tx2gene = tx2gene[,1:2])
    # stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = alpha)

    # drim.padj = getAdjustedPValues(stageRObj, order = FALSE,
    #                             onlySignificantGenes = TRUE)
    # # this summarizes the adjusted p-values from the two-stage analysis. Only genes that passed the filter are included in the table.
    # drim.padj = merge(tx2gene, drim.padj, by.x = c("gene_id","feature_id"),
    #                 by.y = c("geneID","txID"))

    # cat('Screened genes at FDR q<', alpha, ':',
    #     length(unique(drim.padj[drim.padj$gene < alpha,]$gene_id)), '\n')
    # cat('Transcripts passing OFDR:', sum(drim.padj$transcript < alpha), '\n')

    # cat('stageR q <', alpha, '\n')
    # gene_ids = getSignificantGenes(stageRObj)
    # tx_ids = getSignificantTx(stageRObj)
    # if (nrow(tx_ids) > 0) {
    #     print(gene_ids)
    #     print(tx_ids)
    # }
    # cat('Genes where expression switches among isoforms:',
    #     length(unique(drim.padj[drim.padj$transcript < alpha, 'gene_id'])), '\n')

    # # condensing the counts to be converted to proportions
    # drim.prop = reshape2::melt(counts[counts$feature_id %in% proportions(d)$feature_id,], id = c("gene_id", "feature_id"))
    # drim.prop = drim.prop[order(drim.prop$gene_id, drim.prop$variable,
    #                     drim.prop$feature_id),]
    # # Calculate proportions from counts
    # library(dplyr)
    # library(ggplot2)
    # library(ggbeeswarm)
    # drim.prop2 = drim.prop %>%
    #         group_by(gene_id, variable) %>%
    #         mutate(total = sum(value)) %>%
    #         group_by(variable, add=TRUE) %>%
    #         mutate(prop = value/total)
    # drim.prop3 = reshape2::dcast(drim.prop2[,c(1,2,3,6)],
    #                             gene_id + feature_id ~ variable)

    # my_res = list(res.g=res.g, res.t=res.t, res.t.filt = res.t.filt,
    #               dds=d, fm_str=fm_str, drim.padj = drim.padj,
    #               stageRObj=stageRObj, drim.prop=drim.prop3)
    my_res = list(res.g=res.g, res.t=res.t, res.t.filt = res.t.filt,
                  dds=d, fm_str=fm_str)
    return(my_res)
}
```

This should take care of the actual code. Let's load the data and see what we
get:

```r
myregion = 'ACC'

load('~/data/isoforms/tximport_rsem_DTU.RData')
txi = rsem
# I'll just use the metadata from here
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)
# samples has only the metadata now
samples = data[, !grepl(colnames(data), pattern='^ENS')]

# remove samples for the other brain region from the tx counts matrices
keep_me = colnames(txi$counts) %in% samples$submitted_name
for (i in c('abundance', 'counts', 'length')) {
    txi[[i]] = txi[[i]][, keep_me]
}
# sort samples to match order in tximport matrices
rownames(samples) = samples$submitted_name
samples = samples[colnames(txi$counts), ]

# cleaning up some metadata
samples$POP_CODE = as.character(samples$POP_CODE)
samples[samples$POP_CODE=='WNH', 'POP_CODE'] = 'W'
samples[samples$POP_CODE=='WH', 'POP_CODE'] = 'W'
samples$POP_CODE = factor(samples$POP_CODE)
samples$Individual = factor(samples$hbcc_brain_id)
samples[samples$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
samples[samples$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
samples$MoD = factor(samples$Manner.of.Death)
samples$batch = factor(as.numeric(samples$run_date))
samples$Diagnosis = factor(samples$Diagnosis, levels=c('Control', 'Case'))
samples$substance_group = factor(samples$substance_group)
samples$comorbid_group = factor(samples$comorbid_group_update)
samples$evidence_level = factor(samples$evidence_level)
samples$brainbank = factor(samples$bainbank)

# removing everything but autosomes
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "TXNAME"), columns=c('GENEID','TXCHROM'),
               "TXNAME")
# keep only the remaining transcripts and their corresponding genes
txdf.sub = txdf[match(substr(rownames(txi$counts), 1, 15), txdf$TXNAME),]
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('transcript_id', 'transcript_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
tx_meta = merge(txdf.sub, bt_slim, by.x='TXNAME', by.y='transcript_id')
tx_meta$transcript_id = rownames(txi$counts)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
count_matrix = txi$counts[imautosome, ]
tx_meta = tx_meta[imautosome, ]
```

At this point, no filtering has been done, except for keeping only autosomes.
And I added all annotations I wanted. Now, it's just a matter of running the
rest of the analysis.

I had to remove the stageR bit because I was getting some weird errors. I can
always recalculate it later.

Let's run both combinations. First, no BBB:

```r
dtu_acc = list()
for (st in c('lncRNA', 'protein_coding', 'all')) {
    st2 = ifelse(st == 'all', NA, st)
    dtu_acc[[st]] = run_DTU_noPCA(count_matrix, samples, tx_meta,
                                  st2, .05, BBB=T, ncores=7)
}
###
dtu_cau = list()
for (st in c('lncRNA', 'protein_coding', 'all')) {
    st2 = ifelse(st == 'all', NA, st)
    dtu_cau[[st]] = run_DTU_noPCA(count_matrix, samples, tx_meta,
                                  st2, .05, BBB=T, ncores=7)
}
save(dtu_acc, dtu_cau, file='~/data/post_mortem/DTU_03252021_noPCA.RData')
```

I didn't get any genes after filtering for pseudogenes.

# 2021-03-26 16:53:05

I ended up not running stageR because I was getting lots of errors, but I can
always run it posthoc. For now, let's check if the numbers of individual genes
or transcripts varies by our different iterations:

```r
for (s in c('_noPCA', '_BBB_noPCA')) {
    load(sprintf('~/data/post_mortem/DTU_03252021%s.RData', s))
    for (r in c('acc', 'cau')) {
        for (st in c('all', 'protein_coding')) {
            cat('====', s, r, st, '====\n')
            res_str = sprintf('res.g = as.data.frame(dtu_%s$%s$res.g)', r, st)
            eval(parse(text=res_str))
            res_str = sprintf('res.t = as.data.frame(dtu_%s$%s$res.t)', r, st)
            eval(parse(text=res_str))
            res_str = sprintf('res.t.filt = as.data.frame(dtu_%s$%s$res.t.filt)',
                              r, st)
            eval(parse(text=res_str))

            for (alpha in c(.05, .1)) {
                cat('Genes surviving FDR q <', alpha, ':',
                    sum(res.g$adj_pvalue < alpha, na.rm=T), '\n')
                cat('Transcripts surviving FDR q<', alpha, ':',
                    sum(res.t$adj_pvalue < alpha, na.rm=T), '\n')
                cat('Transcripts surviving SD filtering and FDR q<', alpha, ':',
                    sum(res.t.filt$adj_pvalue < alpha, na.rm=T), '\n')
            }
        }
    }
}
```

```
==== _noPCA acc all ====
Genes surviving FDR q < 0.05 : 0 
Transcripts surviving FDR q< 0.05 : 2 
Transcripts surviving SD filtering and FDR q< 0.05 : 2 
Genes surviving FDR q < 0.1 : 0 
Transcripts surviving FDR q< 0.1 : 2 
Transcripts surviving SD filtering and FDR q< 0.1 : 2 

==== _noPCA acc protein_coding ====
Genes surviving FDR q < 0.05 : 0 
Transcripts surviving FDR q< 0.05 : 2 
Transcripts surviving SD filtering and FDR q< 0.05 : 2 
Genes surviving FDR q < 0.1 : 1 
Transcripts surviving FDR q< 0.1 : 3 
Transcripts surviving SD filtering and FDR q< 0.1 : 2 

==== _noPCA cau all ====
Genes surviving FDR q < 0.05 : 27 
Transcripts surviving FDR q< 0.05 : 32 
Transcripts surviving SD filtering and FDR q< 0.05 : 27 
Genes surviving FDR q < 0.1 : 45 
Transcripts surviving FDR q< 0.1 : 65 
Transcripts surviving SD filtering and FDR q< 0.1 : 55 

==== _noPCA cau protein_coding ====
Genes surviving FDR q < 0.05 : 18 
Transcripts surviving FDR q< 0.05 : 20 
Transcripts surviving SD filtering and FDR q< 0.05 : 19 
Genes surviving FDR q < 0.1 : 23 
Transcripts surviving FDR q< 0.1 : 28 
Transcripts surviving SD filtering and FDR q< 0.1 : 27 

==== _BBB_noPCA acc all ====
Genes surviving FDR q < 0.05 : 0 
Transcripts surviving FDR q< 0.05 : 2 
Transcripts surviving SD filtering and FDR q< 0.05 : 2 
Genes surviving FDR q < 0.1 : 0 
Transcripts surviving FDR q< 0.1 : 2 
Transcripts surviving SD filtering and FDR q< 0.1 : 2 

==== _BBB_noPCA acc protein_coding ====
Genes surviving FDR q < 0.05 : 0 
Transcripts surviving FDR q< 0.05 : 3 
Transcripts surviving SD filtering and FDR q< 0.05 : 2 
Genes surviving FDR q < 0.1 : 4 
Transcripts surviving FDR q< 0.1 : 3 
Transcripts surviving SD filtering and FDR q< 0.1 : 2 

==== _BBB_noPCA cau all ====
Genes surviving FDR q < 0.05 : 7 
Transcripts surviving FDR q< 0.05 : 11 
Transcripts surviving SD filtering and FDR q< 0.05 : 8 
Genes surviving FDR q < 0.1 : 25 
Transcripts surviving FDR q< 0.1 : 29 
Transcripts surviving SD filtering and FDR q< 0.1 : 17 

==== _BBB_noPCA cau protein_coding ====
Genes surviving FDR q < 0.05 : 5 
Transcripts surviving FDR q< 0.05 : 6 
Transcripts surviving SD filtering and FDR q< 0.05 : 6 
Genes surviving FDR q < 0.1 : 7 
Transcripts surviving FDR q< 0.1 : 12 
Transcripts surviving SD filtering and FDR q< 0.1 : 7 
```

Potential advantage for protein_coding in ACC, but the other way around for
Caudate. Small advantage for BBB in ACC, big advantage for batch in Caudate! So,
not much here either... I'll just have to pick one and go with it.

## Overall

It seems that BBB was better for ACC in DTE and DTU. 

# 2021-03-29 10:09:19

Let's run the PRS version because it does take a while:

```r
run_DTU_PRS_noPCA = function(count_matrix, tx_meta, data, subtype, prs, alpha,
                             BBB=FALSE, ncores=NA, add_cov=NA) {
    cat('Starting with', nrow(tx_meta), 'variables\n')
    if (is.na(subtype)) {
        keep_me = rep(TRUE, nrow(count_matrix))
    } else {
        keep_me = grepl(tx_meta$transcript_biotype,
                        pattern=sprintf('%s$', subtype))
    }
    cat('Keeping', sum(keep_me), subtype, 'variables\n')
    my_count_matrix = count_matrix[keep_me, ]
    my_tx_meta = tx_meta[keep_me, ]

    # removing variables where more than half of the subjects have zero counts
    keep_me = rowSums(my_count_matrix==0) < .25*ncol(my_count_matrix)
    my_count_matrix = my_count_matrix[keep_me, ]
    cat('Keeping', nrow(my_count_matrix), 'after zero removal\n')

    # removing variables with zero or near-zero variance
    library(caret)
    pp_order = c('zv', 'nzv')
    pp = preProcess(t(my_count_matrix), method = pp_order)
    X = t(predict(pp, t(my_count_matrix)))
    cat('Keeping', nrow(X), 'after NZ and NZV filtering\n')

    # keep only the remaining transcripts and their corresponding genes
    txdf.sub = my_tx_meta[match(rownames(X), my_tx_meta$transcript_id),]
    counts = data.frame(gene_id = txdf.sub$GENEID, feature_id = txdf.sub$TXNAME)
    counts = cbind(counts, X)

    library(DRIMSeq)
    data$group = data$Diagnosis
    data$sample_id = as.character(data$submitted_name)
    d0 = dmDSdata(counts = counts, samples = data)

    n = nrow(data)
    n.small = min(table(data$group))

    d = DRIMSeq::dmFilter(d0,
                        min_samps_feature_expr = n.small, min_feature_expr = 10,
                        min_samps_feature_prop = n.small, min_feature_prop = 0.1,
                        min_samps_gene_expr = n, min_gene_expr = 10)

    countData = round(as.matrix(counts(d)[,-c(1:2)]))
    cat('Keeping', nrow(countData), 'after DRIMSeq expression filtering\n')

    data.pm = data
    # replace the one subject missing population PCs by the median of their
    # self-declared race and ethnicity
    idx = (data.pm$Race.x=='White' & data.pm$Ethnicity.x=='Non-Hispanic' &
           !is.na(data.pm$C1))
    pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
    med_pop = apply(data[idx, pop_pcs], 2, median)
    data.pm[which(is.na(data.pm$C1)), pop_pcs] = med_pop

    if (BBB) {
        data.pm$BBB = factor(sapply(1:nrow(data.pm),
                                    function(x) sprintf('%s_%s',
                                             as.character(data.pm[x,'brainbank']),
                                             as.character(data.pm[x, 'batch']))))
        use_pcs = c('BBB', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI')
    } else {
        use_pcs = c('batch', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI')
    }

    # add more covariates for robustness testing
    if (! is.na(add_cov)) {
        use_pcs = c(use_pcs, add_cov)
    }

    fm_str = sprintf('~ %s + %s', prs, paste0(use_pcs, collapse = ' + '))
    cat('Using formula:', fm_str, '\n')

    # scaling num_vars to assure convergence
    # removed pH because of too many NAs, RIN because we have RINe for everyone
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
                'C1', 'C2', 'C3', 'C4', 'C5', prs)
    for (var in num_vars) {
        data.pm[, var] = scale(data.pm[, var])
    }

    design = model.matrix(as.formula(fm_str), data = data.pm)

    set.seed(42)
    system.time({
        if (is.na(ncores)) {
            d <- dmPrecision(d, design = design)
            d <- dmFit(d, design = design)
            d <- dmTest(d, coef = prs)
        } else {
            multicoreParam <- BiocParallel::MulticoreParam(workers = ncores)
            d <- dmPrecision(d, design = design, BPPARAM=multicoreParam)
            d <- dmFit(d, design = design, BPPARAM=multicoreParam)
            d <- dmTest(d, coef = prs, BPPARAM=multicoreParam)
        }
    })
    res.g = DRIMSeq::results(d)

    return(res.g)
}
```

Then we load the data as usual:

```r
myregion = 'ACC'

load('~/data/isoforms/tximport_rsem_DTU.RData')
txi = rsem
# I'll just use the metadata from here
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
data = data[data$Region==myregion, ]
library(gdata)
more = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid_group_update',
                            'substance_group', 'evidence_level')],
             by='hbcc_brain_id', all.x=T, all.y=F)
# samples has only the metadata now
samples = data[, !grepl(colnames(data), pattern='^ENS')]

# remove samples for the other brain region from the tx counts matrices
keep_me = colnames(txi$counts) %in% samples$submitted_name
for (i in c('abundance', 'counts', 'length')) {
    txi[[i]] = txi[[i]][, keep_me]
}
# sort samples to match order in tximport matrices
rownames(samples) = samples$submitted_name
samples = samples[colnames(txi$counts), ]

# cleaning up some metadata
samples$POP_CODE = as.character(samples$POP_CODE)
samples[samples$POP_CODE=='WNH', 'POP_CODE'] = 'W'
samples[samples$POP_CODE=='WH', 'POP_CODE'] = 'W'
samples$POP_CODE = factor(samples$POP_CODE)
samples$Individual = factor(samples$hbcc_brain_id)
samples[samples$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
samples[samples$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
samples$MoD = factor(samples$Manner.of.Death)
samples$batch = factor(as.numeric(samples$run_date))
samples$Diagnosis = factor(samples$Diagnosis, levels=c('Control', 'Case'))
samples$substance_group = factor(samples$substance_group)
samples$comorbid_group = factor(samples$comorbid_group_update)
samples$evidence_level = factor(samples$evidence_level)
samples$brainbank = factor(samples$bainbank)

# removing everything but autosomes
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "TXNAME"), columns=c('GENEID','TXCHROM'),
               "TXNAME")
# keep only the remaining transcripts and their corresponding genes
txdf.sub = txdf[match(substr(rownames(txi$counts), 1, 15), txdf$TXNAME),]
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('transcript_id', 'transcript_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
tx_meta = merge(txdf.sub, bt_slim, by.x='TXNAME', by.y='transcript_id')
tx_meta$transcript_id = rownames(txi$counts)
imautosome = which(tx_meta$TXCHROM != 'X' &
                   tx_meta$TXCHROM != 'Y' &
                   tx_meta$TXCHROM != 'MT')
count_matrix = txi$counts[imautosome, ]
tx_meta = tx_meta[imautosome, ]

# Grabbing PRS
fname = '~/data/post_mortem/genotyping/1KG/merged_PM_1KG_PRS_12032020.csv'
prs = read.csv(fname)
prs$hbcc_brain_id = sapply(prs$IID,
                          function(x) {
                              br = strsplit(x, '_')[[1]][2];
                              as.numeric(gsub(br, pattern='BR',
                                              replacement=''))})
imWNH = samples$C1 > 0 & samples$C2 < -.075
wnh_brains = samples[which(imWNH),]$hbcc_brain_id
# using the most appropriate PRS, make sure we don't switch subject order
m = merge(samples, prs, by='hbcc_brain_id', sort=F)
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
m[, prs_names] = NA
keep_me = m$hbcc_brain_id %in% wnh_brains
m[keep_me, prs_names] = m[keep_me, 65:76]
m[!keep_me, prs_names] = m[!keep_me, 53:64]
data.prs = m[, c(1:51, 77:88)]
count_matrix = count_matrix[, samples$hbcc_brain_id %in% data.prs$hbcc_brain_id]
data = data.prs
```

And now we go ahead and run the new analysis using PRS as a predictor:

```r
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = list()
for (st in c('lncRNA', 'protein_coding', 'all')) {
    all_res[[st]] = list()
    st2 = ifelse(st == 'all', NA, st)
    for (prs in prs_names) {
        res = run_DTU_PRS_noPCA(count_matrix, tx_meta, data, st2, prs, .05,
                          BBB=T, ncores=30)
        all_res[[st]][[prs]] = res
    }
}
dtuPRS_acc = all_res

all_res = list()
for (st in c('lncRNA', 'protein_coding', 'all')) {
    all_res[[st]] = list()
    st2 = ifelse(st == 'all', NA, st)
    for (prs in prs_names) {
        res = run_DTU_PRS_noPCA(count_matrix, tx_meta, data, st2, prs, .05,
                          BBB=T, ncores=30)
        all_res[[st]][[prs]] = res
    }
}
dtuPRS_cau = all_res









save(dtuPRS_acc, dtuPRS_cau,
     file='~/data/post_mortem/DTU_PRS_03302021_BBB_noPCA.RData')
```

Let's start summarizing the results:

```r
library(GeneOverlap)
load('~/data/post_mortem/DTU_03252021_BBB_noPCA.RData')

all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc')) {
    res.acc = dtu_acc[[subtypes[[st]]]]$res.g
    res.cau = dtu_cau[[subtypes[[st]]]]$res.g
    
    both_res = merge(as.data.frame(res.acc), as.data.frame(res.cau), by='gene_id',
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t, 'gene_id']
        dx_genes = both_res[both_res$pvalue.dx < t, 'gene_id']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        pval2 = NA
        this_res = c(subtypes[[st]], t, 'abs', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/DTU_upDown_overlap_results_03292021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

And repeat it for the stageR p-values:

```r
library(GeneOverlap)
library(stageR)
load('~/data/post_mortem/DTU_03252021_BBB_noPCA.RData')

all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('all', 'pc', 'lnc')) {
    res.acc = dtu_acc[[subtypes[[st]]]]
    res.cau = dtu_cau[[subtypes[[st]]]]

    strp <- function(x) substr(x, 1, 15)
    pScreen = res.acc$res.g$pvalue
    names(pScreen) = strp(res.acc$res.g$gene_id)
    pConfirmation = matrix(res.acc$res.t.filt$pvalue, ncol = 1)
    dimnames(pConfirmation) = list(strp(res.acc$res.t.filt$feature_id),
                                   "transcript")
    tx2gene = data.frame(res.acc$res.t[,c("feature_id", "gene_id")], 
                         res.acc$res.t[,c("feature_id", "gene_id")])
    for (i in 1:2) tx2gene[,i] = strp(tx2gene[,i])
    stageRObj = stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, 
                        pScreenAdjusted = FALSE, tx2gene = tx2gene[,1:2])
    stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = .05)
    df = getAdjustedPValues(stageRObj, onlySignificantGenes=FALSE, order=TRUE)
    df2 = df[, c('geneID', 'gene')]
    df.acc = df2[!duplicated(df2$geneID), ]

    pScreen = res.cau$res.g$pvalue
    names(pScreen) = strp(res.cau$res.g$gene_id)
    pConfirmation = matrix(res.cau$res.t.filt$pvalue, ncol = 1)
    dimnames(pConfirmation) = list(strp(res.cau$res.t.filt$feature_id),
                                   "transcript")
    tx2gene = data.frame(res.cau$res.t[,c("feature_id", "gene_id")], 
                         res.cau$res.t[,c("feature_id", "gene_id")])
    for (i in 1:2) tx2gene[,i] = strp(tx2gene[,i])
    stageRObj = stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, 
                        pScreenAdjusted = FALSE, tx2gene = tx2gene[,1:2])
    stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = .05)
    df = getAdjustedPValues(stageRObj, onlySignificantGenes=FALSE, order=TRUE)
    df2 = df[, c('geneID', 'gene')]
    df.cau = df2[!duplicated(df2$geneID), ]
    
    both_res = merge(as.data.frame(df.acc), as.data.frame(df.cau), by='geneID',
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$gene.prs < t, 'geneID']
        dx_genes = both_res[both_res$gene.dx < t, 'geneID']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        pval2 = NA
        this_res = c(subtypes[[st]], t, 'abs', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('subtype', 'nomPvalThresh', 'direction',
                      'caudateGenes', 'accGenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/DTUstageR_upDown_overlap_results_03292021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```

I also have to run GSEA for DTU:

```r
library(WebGestaltR)

data_dir = '~/data/post_mortem/'
ncpu=32

load('~/data/post_mortem/DTU_03252021_BBB_noPCA.RData')

for (st in c('protein_coding', 'all')) {
    for (myregion in c('acc', 'caudate')) {
        res_str = ifelse(myregion == 'acc',
                         sprintf('res = dtu_acc$%s$res.g', st),
                         sprintf('res = dtu_cau$%s$res.g', st))
        eval(parse(text=res_str))

        tmp2 = data.frame(geneid=res$gene_id, rank=-log(res$pvalue))
        tmp2 = tmp2[order(tmp2$rank, decreasing=T),]

        res_str = sprintf('DTU_%s_%s', myregion, st)
        DBs = c(sprintf('my_%s_sets', myregion))
        for (db in DBs) {
            cat(res_str, db, '\n')
            db_file = sprintf('~/data/post_mortem/%s.gmt', db)
            project_name = sprintf('WG16_%s_%s_10K', res_str, db)
            enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabaseFile=db_file,
                                enrichDatabaseType="genesymbol",
                                interestGene=tmp2,
                                outputDirectory = data_dir,
                                interestGeneType="ensembl_gene_id",
                                sigMethod="top", topThr=20,
                                minNum=3, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800))
        }

        DBs = c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant')
        for (db in DBs) {
            cat(res_str, db, '\n')
            project_name = sprintf('WG16_%s_%s_10K', res_str, db)

            enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase=db,
                                        interestGene=tmp2,
                                        interestGeneType="ensembl_gene_id",
                                        sigMethod="top", topThr=20,
                                        outputDirectory = data_dir,
                                        minNum=5, projectName=project_name,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, perNum=10000))
        }

        DBs = c('KEGG', 'Panther', 'Reactome', 'Wikipathway')
        for (db in DBs) {
            cat(myregion, db, '\n')
            project_name = sprintf('WG16_%s_%s_10K', res_str, db)

            enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase=sprintf('pathway_%s', db),
                                        interestGene=tmp2,
                                        interestGeneType="ensembl_gene_id",
                                        sigMethod="top", minNum=3,
                                        outputDirectory = data_dir,
                                        projectName=project_name,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, topThr=20, perNum=10000))
        }
    }
}
```

I won't bother running the stageR GSEA because the adjusted pvalues screw up the
ranks.

Now, the single gene results:

```r
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
mydir = '~/data/post_mortem/'

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]

library(stageR)
load('~/data/post_mortem/DTU_03252021_BBB_noPCA.RData')
for (r in c('acc', 'cau')) {
    for (st in c('all', 'protein_coding', 'lncRNA')) {
        res_str = sprintf('res = dtu_%s[["%s"]]', r, st)
        eval(parse(text=res_str))
        fname = sprintf('%s/DTU_%s_%s_BBB_annot_03292021.csv', mydir, r, st)

        strp <- function(x) substr(x, 1, 15)
        pScreen = res$res.g$pvalue
        names(pScreen) = strp(res$res.g$gene_id)
        pConfirmation = matrix(res$res.t.filt$pvalue, ncol = 1)
        dimnames(pConfirmation) = list(strp(res$res.t.filt$feature_id),
                                    "transcript")
        tx2gene = data.frame(res$res.t[,c("feature_id", "gene_id")], 
                            res$res.t[,c("feature_id", "gene_id")])
        for (i in 1:2) tx2gene[,i] = strp(tx2gene[,i])
        stageRObj = stageRTx(pScreen = pScreen, pConfirmation = pConfirmation, 
                            pScreenAdjusted = FALSE, tx2gene = tx2gene[,1:2])
        stageRObj = stageWiseAdjustment(stageRObj, method = "dtu", alpha = .05)
        df = getAdjustedPValues(stageRObj, onlySignificantGenes=FALSE, order=TRUE)
        
        colnames(df)[3:4] = c('padj_gene', 'padj_transcript')
        df2 = merge(df, mart, sort=F,
                    by.x='geneID', by.y='ensembl_gene_id', all.x=T, all.y=F)
        df2 = merge(df2, bt_slim, sort=F,
                    by.x='geneID', by.y='gene_id', all.x=T, all.y=F)
        write.csv(df2, row.names=F, file=fname)
    }
}
```

Table of overlap between PRS and DTU:

```r
library(GeneOverlap)
load('~/data/post_mortem/DTU_PRS_03302021_BBB_noPCA.RData')
load('~/data/post_mortem/DTU_03252021_BBB_noPCA.RData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA', pg='pseudogene')
st = 'all'
res.dx = dtu_cau[[subtypes[[st]]]]$res.g
# res.dx = dtu_acc[[subtypes[[st]]]]$res.g
for (p in prs_names) {
    cat(st, p, '\n')
    res_str = sprintf('res.prs = dtuPRS_cau$%s$%s', subtypes[st], p)
    # res_str = sprintf('res.prs = dtuPRS_acc$%s$%s', subtypes[st], p)
    eval(parse(text=res_str))

    both_res = merge(res.dx, res.prs, by='gene_id',
                        all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = both_res[both_res$pvalue.prs < t, 'gene_id']
        dx_genes = both_res[both_res$pvalue.dx < t, 'gene_id']
        go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                    genome.size=nrow(both_res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        pval2 = NA
        this_res = c(subtypes[[st]], p, t, 'abs', length(prs_genes),
                        length(dx_genes), length(inter), pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/DTU_cauUpDown_prs_overlap_results_03302021.csv'
# out_fname = '~/data/post_mortem/DTU_accUpDown_prs_overlap_results_03302021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```
