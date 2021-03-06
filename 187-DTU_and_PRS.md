# 2021-02-04 07:25:56

I'l repeat a similar analysis we ran for DTE and PRS, but this will take
muuuuuch longer because that's how the DTU analysis goes.

```r
run_DTU_PRS = function(count_matrix, tx_meta, data, subtype, prs, alpha, ncores=NA) {
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
                'C1', 'C2', 'C3', 'C4', 'C5', prs, 'pH', 'RIN')
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
for (st in c('lncRNA', 'protein_coding')) {
    all_res[[st]] = list()
    for (prs in prs_names) {
        res = run_DTU_PRS(count_matrix, tx_meta, myregion, st, prs, .05)
        all_res[[st]][[prs]] = res
    }
}
dtuPRS_acc = all_res

all_res = list()
for (st in c('lncRNA', 'protein_coding')) {
    all_res[[st]] = list()
    for (prs in prs_names) {
        res = run_DTU_PRS(count_matrix, tx_meta, myregion, st, prs, .05)
        all_res[[st]][[prs]] = res
    }
}
dtuPRS_cau = all_res

save(dtuPRS_acc_lnc, dtuPRS_acc_pc, dtuPRS_acc_pg,
     dtuPRS_cau_lnc, dtuPRS_cau_pc, dtuPRS_cau_pg,
     file='~/data/post_mortem/DTU_PRS_02032021.RData')
```

We just decided to go without splitting by subtypes... since this takes forever
to run, let's just run all:

```r
ncores=31
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
dtuPRS_acc = list()
for (prs in prs_names) {
    dtuPRS_acc[[prs]] = run_DTU_PRS(count_matrix, tx_meta, data, NA, prs, .05,
                                    ncores=ncores)
}

dtuPRS_cau = list()
for (prs in prs_names) {
    dtuPRS_cau[[prs]] = run_DTU_PRS(count_matrix, tx_meta, data, NA, prs, .05,
                                    ncores=ncores)
}

save(dtuPRS_acc, dtuPRS_cau, file='~/data/post_mortem/DTU_PRS_03082021.RData')
```

Let's also run it for pc and lnc:

```r
ncores=31
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
dtuPRS_acc_pc = list()
for (prs in prs_names) {
    dtuPRS_acc_pc[[prs]] = run_DTU_PRS(count_matrix, tx_meta, data,
                                       'protein_coding', prs, .05,
                                       ncores=ncores)
}
save(dtuPRS_acc_pc, file='~/data/tmp/acc_pc.RData')

dtuPRS_cau_pc = list()
for (prs in prs_names) {
    dtuPRS_cau_pc[[prs]] = run_DTU_PRS(count_matrix, tx_meta, data,
                                       'protein_coding', prs, .05,
                                       ncores=ncores)
}
save(dtuPRS_cau_pc, file='~/data/tmp/cau_pc.RData')
```

And I saved the results for all, lnc, and pc to
~/data/post_mortem/DTU_PRS_03112021.RData. Recall that pseudogenes never leave
anything for analysis after the filters.

Let's then compute the PRS and DX overlaps:

```r
library(GeneOverlap)
load('~/data/post_mortem/DTU_PRS_03112021.RData')
load('~/data/post_mortem/DTU_03082021.RData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
subtypes = list(all='all', pc='protein_coding', lnc='lncRNA')
for (st in c('all', 'pc', 'lnc')) {
    # res.dx = dtu_cau[[subtypes[[st]]]]$res.g
    res.dx = dtu_acc[[subtypes[[st]]]]$res.g
    for (p in prs_names) {
        cat(st, p, '\n')
        # res_str = sprintf('res.prs = dtuPRS_cau_%s[["%s"]]', st, p)
        res_str = sprintf('res.prs = dtuPRS_acc_%s[["%s"]]', st, p)
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
    # out_fname = '~/data/post_mortem/allDTU_caudateUpDown_prs_overlap_results_03112021.csv'
    out_fname = '~/data/post_mortem/allDTU_accUpDown_prs_overlap_results_03112021.csv'
    write.csv(all_res, file=out_fname, row.names=F)
}
```
