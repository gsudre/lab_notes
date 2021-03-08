# 2021-02-03 19:52:27

Let's repeat the PRS analysis we did for DGE, but now for DTE:

```r
run_DTE_PRS = function(count_matrix, tx_meta, data, subtype, prs, alpha) {
    cat('Starting with', nrow(tx_meta), 'variables\n')
    if (is.na(subtype)) {
        keep_me = rep(TRUE, nrow(count_matrix))
    } else {
        keep_me = grepl(tx_meta$transcript_biotype,
                        pattern=sprintf('%s$',subtype))
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

    # checking which PCs are associated with our potential nuiscance variables
    set.seed(42)
    mypca <- prcomp(t(X), scale=TRUE)
    # how many PCs to keep... using Kaiser thredhold, close to eigenvalues < 1
    library(nFactors)
    eigs <- mypca$sdev^2
    nS = nScree(x=eigs)
    keep_me = seq(1, nS$Components$nkaiser)

    mydata = data.frame(mypca$x[, keep_me])
    # create main metadata data frame including metadata and PCs
    data.pm = cbind(data, mydata)
    rownames(data.pm) = data$hbcc_brain_id
    # scaling PRS to make computations easier to converge
    data.pm[, prs] = scale(data.pm[, prs])
    cat('Using', nS$Components$nkaiser, 'PCs from possible', ncol(X), '\n')

    # check which PCs are associated at nominal p<.01
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
                'C1', 'C2', 'C3', 'C4', 'C5', 'pH', 'RIN', prs)
    pc_vars = colnames(mydata)
    num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                       dimnames=list(num_vars, pc_vars))
    num_pvals = num_corrs
    for (x in num_vars) {
        for (y in pc_vars) {
            res = cor.test(data.pm[, x], data.pm[, y], method='spearman')
            num_corrs[x, y] = res$estimate
            num_pvals[x, y] = res$p.value
        }
    }

    categ_vars = c('batch', 'MoD', 'substance_group', 'brainbank',
                'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level')
    categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                         dimnames=list(categ_vars, pc_vars))
    categ_pvals = categ_corrs
    for (x in categ_vars) {
        for (y in pc_vars) {
            res = kruskal.test(data.pm[, y], data.pm[, x])
            categ_corrs[x, y] = res$statistic
            categ_pvals[x, y] = res$p.value
        }
    }
    use_pcs = unique(c(which(num_pvals < .01, arr.ind = T)[, 'col'],
                    which(categ_pvals < .01, arr.ind = T)[, 'col']))
    # only use the ones not related to Diagnosis
    keep_me = c()
    for (pc in use_pcs) {
        keep_me = c(keep_me, num_pvals[prs, pc] > .05)
    }
    use_pcs = use_pcs[keep_me]
    
    fm_str = sprintf('~ %s + %s', prs, paste0(pc_vars[use_pcs],
                                                collapse = ' + '))
    cat('Found', length(use_pcs), 'PCs p < .01\n')
    cat('Using formula:', fm_str, '\n')

    # scaling PCs to assure convergence
    for (var in pc_vars[use_pcs]) {
        data.pm[, var] = scale(data.pm[, var])
    }

    # removing variables with low expression
    library(edgeR)
    design=model.matrix(as.formula(fm_str), data=data.pm)
    isexpr <- filterByExpr(X, design=design)
    countsExpr = X[isexpr,]
    metaExpr = data.frame(GENEID = substr(rownames(countsExpr), 1, 15))
    metaExpr = merge(metaExpr, my_tx_meta, by='GENEID', sort=F)
    cat('Keeping', nrow(countsExpr), 'after expression filtering\n')

    # preparing DESeqData and running main analysis
    countdata = round(countsExpr)
    colnames(countdata) = rownames(data.pm)
    library(DESeq2)
    # because DESeq doesn't remove outliers if there are continuous variables
    # in the formula, we need to do this iteratively
    nOutliers = Inf
    myCounts = round(countsExpr)
    while (nOutliers > 0) {
        dds <- DESeqDataSetFromMatrix(countData = myCounts,
                                    colData = data.pm,
                                    design = as.formula(fm_str))
        cat('Processing', nrow(dds), 'variables.\n')
        dds <- DESeq(dds)
        maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
        # outlier cut-off uses the 99% quantile of the F(p,m-p) distribution (with 
        # p the number of parameters including the intercept and m number of
        # samples).
        m <- ncol(dds)
        # number or parameters (PCs + Diagnosis + intercept)
        p <- length(use_pcs) + 2
        co = qf(.99, p, m - p)
        keep_me = which(maxCooks < co)
        nOutliers = nrow(myCounts) - length(keep_me)
        cat('Found', nOutliers, 'outliers.\n')
        myCounts = round(myCounts)[keep_me, ]
    }
    res <- results(dds, name = prs, alpha = alpha)
    cat(sprintf('FDR q < %.2f\n', alpha))
    print(summary(res))

    library(IHW)
    resIHW <- results(dds, name = prs, alpha = alpha, filterFun=ihw)
    cat(sprintf('IHW q < %.2f\n', alpha))
    print(summary(resIHW))

    return(resIHW)
}
```

Then we load the data as usual:

```r
myregion = 'ACC'

load('~/data/isoforms/tximport_rsem_DTE.RData')
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
for (prs in prs_names) {
    res = run_DTE_PRS(count_matrix, tx_meta, myregion, 'lncRNA', prs, .05)
    all_res[[prs]] = res
}
dtePRS_acc_lnc = all_res
# protein_coding, lncRNA, pseudogene

all_res = list()
for (prs in prs_names) {
    res = run_DTE_PRS(count_matrix, tx_meta, myregion, 'lncRNA', prs, .05)
    all_res[[prs]] = res
}
dtePRS_cau_lnc = all_res

save(dtePRS_acc_lnc, dtePRS_acc_pc, dtePRS_acc_pg,
     dtePRS_cau_lnc, dtePRS_cau_pc, dtePRS_cau_pg,
     file='~/data/post_mortem/DTE_PRS_02032021.RData')
```

And we have to re-run the skinny DTE to get the actual res, from 183:

```r
dte_acc = list() 
for (st in c('lncRNA', 'protein_coding', 'pseudogene')) {
    dte_acc[[st]] = run_skinny_DTE(count_matrix, tx_meta, myregion, st)
}
dte_cau = list() 
for (st in c('lncRNA', 'protein_coding', 'pseudogene')) {
    dte_cau[[st]] = run_skinny_DTE(count_matrix, tx_meta, myregion, st)
}
save(dte_acc, dte_cau, file='~/data/post_mortem/DTE_skinny_02032021.RData')
```

Let's then compute the PRS and DX overlaps:

```r
library(GeneOverlap)
load('~/data/post_mortem/DTE_PRS_02032021.RData')
load('~/data/post_mortem/DTE_skinny_02032021.RData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
subtypes = list(pc='protein_coding', lnc='lncRNA', pg='pseudogene')
for (st in c('pc', 'lnc', 'pg')) {
    res.dx = as.data.frame(dte_cau[[subtypes[[st]]]])
    for (p in prs_names) {
        cat(st, p, '\n')
        res_str = sprintf('res.prs = as.data.frame(dtePRS_cau_%s[["%s"]])', st, p)
        eval(parse(text=res_str))

        both_res = merge(as.data.frame(res.dx), as.data.frame(res.prs), by=0,
                         all.x=F, all.y=F, suffixes = c('.dx', '.prs'))
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs > 0,
                                 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx > 0,
                                'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            allUp = union(both_res[both_res$stat.prs > 0, 'Row.names'],
                          both_res[both_res$stat.dx > 0, 'Row.names'])
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
            go.obj <- testGeneOverlap(go.obj)
            pval2 = getPval(go.obj)
            this_res = c(subtypes[[st]], p, t, 'up', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t & both_res$stat.prs < 0,
                                 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t & both_res$stat.dx < 0,
                                'Row.names']
            go.obj <- newGeneOverlap(prs_genes, dx_genes,
                                     genome.size=nrow(both_res))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval1 = getPval(go.obj)
            allDown = union(both_res[both_res$stat.prs < 0, 'Row.names'],
                            both_res[both_res$stat.dx < 0, 'Row.names'])
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
            go.obj <- testGeneOverlap(go.obj)
            pval2 = getPval(go.obj)
            this_res = c(subtypes[[st]], p, t, 'down', length(prs_genes),
                         length(dx_genes), length(inter), pval1, pval2)
            all_res = rbind(all_res, this_res)
        }
        for (t in c(.05, .01, .005, .001)) {
            prs_genes = both_res[both_res$pvalue.prs < t, 'Row.names']
            dx_genes = both_res[both_res$pvalue.dx < t, 'Row.names']
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
}
colnames(all_res) = c('subtype', 'PRS', 'nomPvalThresh', 'direction',
                      'PRSgenes', 'PMgenes', 'overlap', 'pvalWhole',
                      'pvalDirOnly')
out_fname = '~/data/post_mortem/allDTE_caudateUpDown_prs_overlap_results_02032021.csv'
write.csv(all_res, file=out_fname, row.names=F)
```