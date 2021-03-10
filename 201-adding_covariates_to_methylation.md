# 2021-03-10 09:57:12

Looking over Alex's slides again, I noticed he tested for some covariates I
didn't have in my gene expression data, so it might be worth checking out
whether they'd require removing different PCs. I'm using the code from note 196
as my base.

```r
run_methyl = function(mVals, samples, subtype, ann450kSub) {
    cat('Starting with', nrow(mVals), 'variables\n')
    if (subtype == 'all') {
        keep_me = rep(TRUE, nrow(mVals))
    } else {
      keep_me = grepl(ann450kSub$Relation_to_Island,
                      pattern=sprintf('%s$', subtype))
    }
    cat('Keeping', sum(keep_me), subtype, 'variables\n')
    mVals = mVals[keep_me, ]
    ann450kSub = ann450kSub[keep_me, ]

    # removing variables with zero or near-zero variance
    library(caret)
    pp_order = c('zv', 'nzv')
    pp = preProcess(t(mVals), method = pp_order)
    X = t(predict(pp, t(mVals)))
    cat('Keeping', nrow(X), 'after NZ and NZV filtering\n')

    # remove the 2 probes with infinity
    bad_probes = rownames(which(abs(mVals)==Inf, arr.ind = T))
    X = X[!(rownames(X) %in% bad_probes), ]

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
    data.pm = cbind(samples, mydata)
    rownames(data.pm) = samples$hbcc_brain_id
    cat('Using', nS$Components$nkaiser, 'PCs from possible', ncol(X), '\n')

    # check which PCs are associated at nominal p<.01
    qc_vars = c("Restoration", "Staining.Green", "Staining.Red",
                "Extension.Green", "Extension.Red", "Hybridization.High.Medium",
                "Hybridization.Medium.Low", "Target.Removal.1",
                "Target.Removal.2", "Bisulfite.Conversion.I.Green",
                "Bisulfite.Conversion.I.Red", "Bisulfite.Conversion.II",
                "Specificity.I.Green", "Specificity.I.Red", "Specificity.II",
                "Non.polymorphic.Green", "Non.polymorphic.Red")
    num_vars = c('Age', 'PMI', 'C1', 'C2', 'C3', 'C4', 'C5', 'pH', qc_vars)
    pc_vars = colnames(mydata)
    num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                    dimnames=list(num_vars, pc_vars))
    num_pvals = num_corrs
    for (x in num_vars) {
        for (y in pc_vars) {
            res = cor.test(samples[, x], mydata[, y], method='spearman')
            num_corrs[x, y] = res$estimate
            num_pvals[x, y] = res$p.value
        }
    }
    categ_vars = c('Diagnosis', 'MoD', 'substance_group', 'brain_bank',
                'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level', 'Kit',
                'Sample_Plate', 'Sample_Group', 'Sample_Well', 'Scanner',
                'Sentrix_ID', 'Sentrix_Position')
    categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                    dimnames=list(categ_vars, pc_vars))
    categ_pvals = categ_corrs
    for (x in categ_vars) {
        for (y in pc_vars) {
            res = kruskal.test(mydata[, y], samples[, x])
            categ_corrs[x, y] = res$statistic
            categ_pvals[x, y] = res$p.value
        }
    }
    use_pcs = unique(c(which(num_pvals < .01, arr.ind = T)[, 'col'],
                    which(categ_pvals < .01, arr.ind = T)[, 'col']))
    # only use the ones not related to Diagnosis
    keep_me = c()
    for (pc in use_pcs) {
        keep_me = c(keep_me, categ_pvals['Diagnosis', pc] > .05)
    }
    use_pcs = use_pcs[keep_me]
    fm_str = sprintf('~ Diagnosis + %s', paste0(pc_vars[use_pcs],
                                                collapse = ' + '))
    cat('Found', length(use_pcs), 'PCs p < .01\n')
    cat('Using formula:', fm_str, '\n')

    # scaling PCs to assure convergence
    for (var in pc_vars[use_pcs]) {
        data.pm[, var] = scale(data.pm[, var])
    }

    library(limma)
    design = model.matrix(as.formula(fm_str), data.pm)
    fit <- lmFit(X, design)
    fit2 <- eBayes( fit )

    DMPs <- topTable(fit2, num=Inf, coef='DiagnosisCase', genelist=ann450kSub)
    library(missMethyl)
    fitvar <- varFit(X, design = design, coef = c(1,2))
    topDV <- topVar(fitvar, coef=2, number=nrow(X))
    m = merge(topDV, ann450kSub, by=0, sort=F)

    res = list(DMPs = DMPs, topDV=m, fm_str=fm_str, design=design,
               pcs = rbind(categ_pvals, num_pvals))
    return(res)
}
```

The thing here is that we need to add a few more covariates. Specifically, they
come from control variables and some stuff that Alex dug out in his metadata
matrix. So, let's load those data:

```r
rootDir <- "~/data/methylation_post_mortem"
rawDataDir <- paste(rootDir, '/Shaw2019', sep="")
samplesFile <- paste(rootDir, '/samples.csv', sep="")

# disable scientific notation
samples <- read.csv(samplesFile)
samples$Sample_Group <- as.factor(as.character(samples$Diagnosis))
samples$Basename <- paste0(rawDataDir, '/',  samples$Sentrix_ID, '/',
                           samples$Sentrix_ID,'_',samples$Sentrix_Position)
samples$Sentrix_ID <- as.factor(samples$Sentrix_ID)
samples$Sample_Group <- as.factor(samples$Sample_Group)
samples$Row <- as.factor(substr(samples$Sentrix_Position, 1, 3))
samples$Scanner <- as.factor(samples$Scanner)
samples$Region <- as.factor(samples$Region)
samples$Sex <- as.factor(samples$Sex)
samples$Kit <- as.factor(samples$Kit)
samples$Sample_Plate <- as.factor(samples$Sample_Plate)

library(ewastools)
raw = read_idats(samples$Basename, quiet = FALSE)
ctrls = control_metrics(raw)
df_ctrls = as.data.frame(ctrls)
rownames(df_ctrls) <- samples$Sample_Name
more_covars = cbind(samples, df_ctrls)
saveRDS(more_covars, file='~/data/methylation_post_mortem/more_covars.rds')
```

Let's use some of the work we've done in the past:

```r
library(minfi)
load('~/data/methylation_post_mortem/filt_ACC_02182021.RData')
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
# and b values for display and interpretation
bVals <- getBeta(mSetSqFlt)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
```

Now it's just a mater of adding those variables to our sample:

```r
more_covars = readRDS('~/data/methylation_post_mortem/more_covars.rds')
more_covars = more_covars[more_covars$Region == 'ACC',]
samples2 = merge(samples, more_covars[, c('hbcc_brain_id', "Sentrix_ID",
                                          "Sentrix_Position", 'Kit',
                                          colnames(more_covars)[26:42])],
                by='hbcc_brain_id', all.x=T, all.y=F)

res = run_methyl(mVals, samples2, 'all', ann450kSub)
```

Hum... so, it does change things, because different PCs were picked... let's see
what changes in the results.

```r
res_acc = list()
for (st in c('all', 'Island', 'Shelf', 'Shore', 'Sea')) {
    res_acc[[st]] = run_methyl(mVals, samples2, st, ann450kSub)
}
idx = ann450kSub$Enhancer == "TRUE"
res_acc[['enhancer']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="Body") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon"))
res_acc[['body']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200"))
res_acc[['promoter1']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="5\'UTR") )
res_acc[['promoter2']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
save(res_acc, file='~/data/methylation_post_mortem/res_ACC_03102021.RData')
```

Now, let's repeat the same stuff for the Caudate:

```r
library(minfi)
load('~/data/methylation_post_mortem/filt_Caudate_02222021.RData')
# calculate M-values for statistical analysis
mVals <- getM(mSetSqFlt)
# and b values for display and interpretation
bVals <- getBeta(mSetSqFlt)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]
more_covars = readRDS('~/data/methylation_post_mortem/more_covars.rds')
more_covars = more_covars[more_covars$Region == 'Caudate',]
samples2 = merge(samples, more_covars[, c('hbcc_brain_id', "Sentrix_ID",
                                          "Sentrix_Position", 'Kit',
                                          colnames(more_covars)[26:42])],
                by='hbcc_brain_id', all.x=T, all.y=F)

res_cau = list()
for (st in c('all', 'Island', 'Shelf', 'Shore', 'Sea')) {
    res_cau[[st]] = run_methyl(mVals, samples2, st, ann450kSub)
}
idx = ann450kSub$Enhancer == "TRUE"
res_cau[['enhancer']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="Body") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon"))
res_cau[['body']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200"))
res_cau[['promoter1']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="5\'UTR") )
res_cau[['promoter2']] = run_methyl(mVals[idx,], samples2, 'all', ann450kSub[idx,])
save(res_cau, file='~/data/methylation_post_mortem/res_Caudate_03102021.RData')
```

## PRS

Because I added a few more covariates, we need to also add them to the PRS
analysis:

```r
run_methyl_PRS = function(mVals, samples, subtype, ann450kSub, prs) {
    cat('Starting with', nrow(mVals), 'variables\n')
    if (subtype == 'all') {
        keep_me = rep(TRUE, nrow(mVals))
    } else {
      keep_me = grepl(ann450kSub$Relation_to_Island,
                      pattern=sprintf('%s$', subtype))
    }
    cat('Keeping', sum(keep_me), subtype, 'variables\n')
    mVals = mVals[keep_me, ]
    ann450kSub = ann450kSub[keep_me, ]

    # removing variables with zero or near-zero variance
    library(caret)
    pp_order = c('zv', 'nzv')
    pp = preProcess(t(mVals), method = pp_order)
    X = t(predict(pp, t(mVals)))
    cat('Keeping', nrow(X), 'after NZ and NZV filtering\n')

    # remove the 2 probes with infinity
    bad_probes = rownames(which(abs(mVals)==Inf, arr.ind = T))
    X = X[!(rownames(X) %in% bad_probes), ]

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
    data.pm = cbind(samples, mydata)
    rownames(data.pm) = samples$hbcc_brain_id
    # scaling PRS to make computations easier to converge
    data.pm[, prs] = scale(data.pm[, prs])
    cat('Using', nS$Components$nkaiser, 'PCs from possible', ncol(X), '\n')

    # check which PCs are associated at nominal p<.01
    qc_vars = c("Restoration", "Staining.Green", "Staining.Red",
                "Extension.Green", "Extension.Red", "Hybridization.High.Medium",
                "Hybridization.Medium.Low", "Target.Removal.1",
                "Target.Removal.2", "Bisulfite.Conversion.I.Green",
                "Bisulfite.Conversion.I.Red", "Bisulfite.Conversion.II",
                "Specificity.I.Green", "Specificity.I.Red", "Specificity.II",
                "Non.polymorphic.Green", "Non.polymorphic.Red")
    num_vars = c('Age', 'PMI', 'C1', 'C2', 'C3', 'C4', 'C5', 'pH', qc_vars, prs)
    pc_vars = colnames(mydata)
    num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                    dimnames=list(num_vars, pc_vars))
    num_pvals = num_corrs
    for (x in num_vars) {
        for (y in pc_vars) {
            res = cor.test(samples[, x], mydata[, y], method='spearman')
            num_corrs[x, y] = res$estimate
            num_pvals[x, y] = res$p.value
        }
    }
    categ_vars = c('MoD', 'substance_group', 'brain_bank',
                'comorbid_group', 'POP_CODE', 'Sex', 'evidence_level', 'Kit',
                'Sample_Plate', 'Sample_Group', 'Sample_Well', 'Scanner',
                'Sentrix_ID', 'Sentrix_Position')
    categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                    dimnames=list(categ_vars, pc_vars))
    categ_pvals = categ_corrs
    for (x in categ_vars) {
        for (y in pc_vars) {
            res = kruskal.test(mydata[, y], samples[, x])
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

    library(limma)
    design = model.matrix(as.formula(fm_str), data.pm)
    fit <- lmFit(X, design)
    fit2 <- eBayes( fit )

    DMPs <- topTable(fit2, num=Inf, coef=prs, genelist=ann450kSub)
    library(missMethyl)
    fitvar <- varFit(X, design = design, coef = c(1,2))
    topDV <- topVar(fitvar, coef=2, number=nrow(X))
    m = merge(topDV, ann450kSub, by=0, sort=F)

    res = list(DMPs = DMPs, topDV=m, fm_str=fm_str, design=design,
               pcs = rbind(categ_pvals, num_pvals))
    return(res)
}
```

```r
library(minfi)
load('~/data/methylation_post_mortem/filt_ACC_02182021.RData')
mVals = getM(mSetSqFlt)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

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
m[keep_me, prs_names] = m[keep_me, 83:94]
m[!keep_me, prs_names] = m[!keep_me, 71:82]
data.prs = m[, c(1:69, 95:106)]

mVals.prs = mVals[, samples$hbcc_brain_id %in% data.prs$hbcc_brain_id]

more_covars = readRDS('~/data/methylation_post_mortem/more_covars.rds')
more_covars = more_covars[more_covars$Region == 'ACC',]
samples2 = merge(data.prs, more_covars[, c('hbcc_brain_id', "Sentrix_ID",
                                          "Sentrix_Position", 'Kit',
                                          colnames(more_covars)[26:42])],
                 by='hbcc_brain_id', all.x=T, all.y=F)

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = list()
for (prs in prs_names) {
    st_res = list()
    for (st in c('all', 'Island', 'Shelf', 'Shore', 'Sea')) {
        st_res[[st]] = run_methyl_PRS(mVals.prs, samples2, st, ann450kSub, prs)
    }
    idx = ann450kSub$Enhancer == "TRUE"
    st_res[['enhancer']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                          ann450kSub[idx,])
    idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="Body") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon"))
    st_res[['body']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                      ann450kSub[idx,])
    idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200"))
    st_res[['promoter1']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                          ann450kSub[idx,])
    idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="5\'UTR") )
    st_res[['promoter2']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                          ann450kSub[idx,])
    all_res[[prs]] = st_res
}
save(all_res, prs_names,
     file='~/data/methylation_post_mortem/res_ACC_PRS_03102021.RData')
```

As usual, repeat it for the caudate:

```r
library(minfi)
load('~/data/methylation_post_mortem/filt_Caudate_02222021.RData')
mVals = getM(mSetSqFlt)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

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
m[keep_me, prs_names] = m[keep_me, 83:94]
m[!keep_me, prs_names] = m[!keep_me, 71:82]
data.prs = m[, c(1:69, 95:106)]

mVals.prs = mVals[, samples$hbcc_brain_id %in% data.prs$hbcc_brain_id]

more_covars = readRDS('~/data/methylation_post_mortem/more_covars.rds')
more_covars = more_covars[more_covars$Region == 'Caudate',]
samples2 = merge(data.prs, more_covars[, c('hbcc_brain_id', "Sentrix_ID",
                                          "Sentrix_Position", 'Kit',
                                          colnames(more_covars)[26:42])],
                 by='hbcc_brain_id', all.x=T, all.y=F)

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = list()
for (prs in prs_names) {
    st_res = list()
    for (st in c('all', 'Island', 'Shelf', 'Shore', 'Sea')) {
        st_res[[st]] = run_methyl_PRS(mVals.prs, samples2, st, ann450kSub, prs)
    }
    idx = ann450kSub$Enhancer == "TRUE"
    st_res[['enhancer']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                          ann450kSub[idx,])
    idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="Body") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon"))
    st_res[['body']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                      ann450kSub[idx,])
    idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200"))
    st_res[['promoter1']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                          ann450kSub[idx,])
    idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon") |
        grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="5\'UTR") )
    st_res[['promoter2']] = run_methyl_PRS(mVals.prs[idx,], samples2, 'all',
                                          ann450kSub[idx,])
    all_res[[prs]] = st_res
}
save(all_res, prs_names,
     file='~/data/methylation_post_mortem/res_Caudate_PRS_03102021.RData')
```

## GSEA

As the results start to come in, let's set up GSEA:

```r
# myregion='ACC'
myregion='Caudate'
library(methylGSA)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
load(sprintf('~/data/methylation_post_mortem/res_%s_03102021.RData', myregion))
outdir = '~/data/methylation_post_mortem/more_covs/'
res_str = ifelse(myregion=='ACC', 'res = res_acc[["all"]]',
                 'res = res_cau[["all"]]')
eval(parse(text=res_str))

ranks = res$DMPs[, 'P.Value']
names(ranks) = rownames(res$DMPs)
for (gs in c('GO', 'KEGG', 'Reactome')) {
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(gs, g, '\n')
        res = methylglm(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.type = gs, parallel=F)
        fname = sprintf('%s/%s_%s_DMP_glm_%s.csv', outidr, myregion, g, gs)
        write.csv(res, row.names=F, file=fname)
    }
}

for (gs in c('GO', 'KEGG', 'Reactome')) {
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(gs, g, '\n')
        res = methylRRA(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.type = gs, method='GSEA')
        fname = sprintf('%s/%s_%s_DMP_RRA_%s.csv', outdir, myregion, g, gs)
        write.csv(res, row.names=F, file=fname)
    }
}

ranks = res$topDV[, 'P.Value']
names(ranks) = res$topDV$Row.names
for (gs in c('GO', 'KEGG', 'Reactome')) {
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(gs, g, '\n')
        res = methylglm(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.type = gs, parallel=F)
        fname = sprintf('%s/%s_%s_topVar_glm_%s.csv', outdir, myregion, g, gs)
        write.csv(res, row.names=F, file=fname)
    }
}
for (gs in c('GO', 'KEGG', 'Reactome')) {
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(gs, g, '\n')
        res = methylRRA(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.type = gs, method='GSEA')
        fname = sprintf('%s/%s_%s_topVar_RRA_%s.csv', outdir, myregion, g, gs)
        write.csv(res, row.names=F, file=fname)
    }
}
```

Let's try our own sets for GSEA too:

```r
library(WebGestaltR)
db_file = sprintf('~/data/post_mortem/my_%s_sets.gmt', myregion)
gmt = readGmt(db_file) # already in gene symbols
sets = list()
for (d in unique(gmt$description)) {
    genes = gmt[gmt$description==d, 'gene']
    sets[[d]] = genes
}

ranks = res$DMPs[, 'P.Value']
names(ranks) = rownames(res$DMPs)
for (g in c('all', 'body', 'promoter1', 'promoter2')) {
    cat(g, '\n')
    res2 = methylglm(cpg.pval = ranks, minsize = 3, group=g,
                maxsize = 500, GS.list = sets, parallel=F)
    fname = sprintf('%s/%s_%s_DMP_glm_myset.csv', outdir, myregion, g)
    write.csv(res2, row.names=F, file=fname)
}
for (g in c('all', 'body', 'promoter1', 'promoter2')) {
    cat(g, '\n')
    res2 = methylRRA(cpg.pval = ranks, minsize = 3, group=g,
                maxsize = 500, GS.list = sets, method='GSEA')
    fname = sprintf('%s/%s_%s_DMP_RRA_myset.csv', outdir, myregion, g)
    write.csv(res2, row.names=F, file=fname)
}

ranks = res$topDV[, 'P.Value']
names(ranks) = res$topDV$Row.names
for (g in c('all', 'body', 'promoter1', 'promoter2')) {
    cat(g, '\n')
    res2 = methylglm(cpg.pval = ranks, minsize = 3, group=g,
                maxsize = 500, GS.list = sets, parallel=F)
    fname = sprintf('%s/%s_%s_topVar_glm_myset.csv', outdir, myregion, g)
    write.csv(res2, row.names=F, file=fname)
}
for (g in c('all', 'body', 'promoter1', 'promoter2')) {
    cat(g, '\n')
    res2 = methylRRA(cpg.pval = ranks, minsize = 3, group=g,
                maxsize = 500, GS.list = sets, method='GSEA')
    fname = sprintf('%s/%s_%s_topVar_RRA_myset.csv', outdir, myregion, g)
    write.csv(res2, row.names=F, file=fname)
}
```

Let's see what we get if we run the same GO sets as we did for DGE and DTE:

```r
library(WebGestaltR)

DBs = c('Biological_Process', 'Cellular_Component', 'Molecular_Function')
for (db in DBs) {
    db_file = sprintf('~/data/post_mortem/hsapiens_geneontology_%s_noRedundant_entrezgene.gmt', db)
    gmt = readGmt(db_file) # already in gene symbols
    sets = list()
    for (d in unique(gmt$description)) {
        genes = gmt[gmt$description==d, 'gene']
        sets[[d]] = genes
    }

    ranks = res$DMPs[, 'P.Value']
    names(ranks) = rownames(res$DMPs)
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(g, db, '\n')
        res = methylglm(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.list = sets, parallel=F,
                    GS.idtype='ENTREZID')
        fname = sprintf('%s/%s_%s_DMP_glm_%s.csv', outdir, myregion, g, db)
        write.csv(res, row.names=F, file=fname)
    }
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(g, db, '\n')
        res = methylRRA(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.list = sets, method='GSEA',
                    GS.idtype='ENTREZID')
        fname = sprintf('%s/%s_%s_DMP_RRA_%s.csv', outdir, myregion, g, db)
        write.csv(res, row.names=F, file=fname)
    }
}

for (db in DBs) {
    db_file = sprintf('~/data/post_mortem/hsapiens_geneontology_%s_noRedundant_entrezgene.gmt', db)
    gmt = readGmt(db_file) # already in gene symbols
    sets = list()
    for (d in unique(gmt$description)) {
        genes = gmt[gmt$description==d, 'gene']
        sets[[d]] = genes
    }

    ranks = res$topDV[, 'P.Value']
    names(ranks) = res$topDV$Row.names
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(g, db, '\n')
        res = methylglm(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.list = sets, parallel=F,
                    GS.idtype='ENTREZID')
        fname = sprintf('%s/%s_%s_topVar_glm_%s.csv', outdir, myregion, g, db)
        write.csv(res, row.names=F, file=fname)
    }
    for (g in c('all', 'body', 'promoter1', 'promoter2')) {
        cat(g, db, '\n')
        res = methylRRA(cpg.pval = ranks, minsize = 3, group=g,
                    maxsize = 500, GS.list = sets, method='GSEA',
                    GS.idtype='ENTREZID')
        fname = sprintf('%s/%s_%s_topVar_RRA_%s.csv', myregion, g, db)
        write.csv(res, row.names=F, file=fname)
    }
}
```

# TODO
* PRS
* GSEA
* spit out results











Like the ACC analysis, let's run GSEA:


Let's see if there is anything interesting here. Starting with our own sets:

 * GWAS at 0.038961039 for promoter1_DMP_RRA and 0.016987381 for
   promoter1_DMP_glm, somewhat confirmatory. Both nominal though.
 * Not much interesting stuff happening though. adult dev at .009 for
   all_topVar_glm, carried by body at .016 all nominally, somewhat unlikely to
   survive correction without TWAS and GWAS. The body result is confirmed by
   topVar_RRA at 0.002016416 nominally, which could potentially survive. It's
   topVar though...

Biological Functions: (q < .05)
 * only hits at all_DMP_RRA: cell-cell adhesion mediated by cadherin.
   promoter2_DMP_RRA also had a couple hits, but not an intersection:
   columnar/cuboidal epithelial cell differentiation, and response to estradiol.

Cellular Components: (q<.05)
 * all_DMP_RRA: non-motile cilium
 * everything else was for topVar: all_topVar_RRA got axon, body_topVar_RRA got
   cell-substrate junction and receptor complex, while its glm counterpart got
   myelin sheath. It's interesting, but topVar and a single hit for FDR q < .05.

Molecular Processes: (q < .05)
 * all_DMP_glm got oxidoreductase activity, acting on single donors with incorporation of molecular oxygen
 * all_topVar_RRA got bHLH transcription factor binding
 * nothing that interesting...

KEGG: (q < .05)
 * all_DMP_RRA got Prostate cancer
 * body_DMP_RRA got Notch signaling pathway
 * body_topVar_glm got Leukocyte transendothelial migration and Fat digestion and absorption
 * again, nothing great

Reactome: (q < .05)
 * multiple hits here... let's see
 * all_DMP_glm:
Acyl chain remodelling of PI
Homo sapiens: Keratinization
Homo sapiens: Classical Kir channels
Homo sapiens: Acyl chain remodelling of PS
Homo sapiens: Formation of the cornified envelope
 * all_DMP_RRA:
Homo sapiens: Signaling by Type 1 Insulin-like Growth Factor 1 Receptor (IGF1R)
Homo sapiens: FRS-mediated FGFR4 signaling
Homo sapiens: Myogenesis
Homo sapiens: Parasite infection
Homo sapiens: Leishmania phagocytosis
Homo sapiens: FCGR3A-mediated phagocytosis
Homo sapiens: IGF1R signaling cascade
Homo sapiens: PI-3K cascade:FGFR3
Homo sapiens: Formation of Senescence-Associated Heterochromatin Foci (SAHF)
Homo sapiens: FRS-mediated FGFR3 signaling
Homo sapiens: PI-3K cascade:FGFR4
Homo sapiens: Pre-NOTCH Processing in Golgi
Homo sapiens: Adherens junctions interactions
Homo sapiens: SHC-mediated cascade:FGFR3
Homo sapiens: FGFR4 ligand binding and activation
Homo sapiens: IRS-related events triggered by IGF1R
Homo sapiens: Integration of provirus
Homo sapiens: SHC-mediated cascade:FGFR4
 * body_DMP_glm: Homo sapiens: Keratinization
 * body_DMP_RRA:
Homo sapiens: FRS-mediated FGFR4 signaling
Homo sapiens: Downstream signaling of activated FGFR3
Homo sapiens: Downstream signaling of activated FGFR4
Homo sapiens: Signaling by FGFR3 in disease
Homo sapiens: Signaling by FGFR3 point mutants in cancer
Homo sapiens: PI-3K cascade:FGFR4
Homo sapiens: FRS-mediated FGFR3 signaling
Homo sapiens: Uptake and function of anthrax toxins
Homo sapiens: Regulation of RUNX1 Expression and Activity
Homo sapiens: SHC-mediated cascade:FGFR4
Homo sapiens: Reduction of cytosolic Ca++ levels
Homo sapiens: PI-3K cascade:FGFR3
Homo sapiens: Signaling by FGFR3
Homo sapiens: FRS-mediated FGFR2 signaling
Homo sapiens: Netrin-1 signaling
Homo sapiens: Signaling by FGFR4
Homo sapiens: Phospholipase C-mediated cascade; FGFR4
Homo sapiens: Downstream signaling of activated FGFR2
 * promoter2_DMP_glm:
Homo sapiens: Classical Kir channels
Homo sapiens: Hydrolysis of LPC
Homo sapiens: Acyl chain remodelling of PI
 * all_topVar_glm:
Homo sapiens: GLI proteins bind promoters of Hh responsive genes to promote transcription
Homo sapiens: RUNX1 regulates transcription of genes involved in BCR signaling
 * body_topVar_glm: Homo sapiens: Proton-coupled monocarboxylate transport
 * body_topVar_RRA:
Homo sapiens: HDMs demethylate histones
Homo sapiens: Interleukin-23 signaling
Homo sapiens: Lysosomal oligosaccharide catabolism
Homo sapiens: Presynaptic depolarization and calcium channel opening
 * promoter1_topVar_glm: Homo sapiens: RUNX1 regulates transcription of genes involved in interleukin signaling
 * promoter2_topVar_glm: Homo sapiens: RUNX1 regulates transcription of genes involved in interleukin signaling

FGF3 is fibroblast growth factor receptor 3. These proteins play a role in
several important cellular processes, including regulation of cell growth and
division (proliferation), determination of cell type, formation of blood vessels
(angiogenesis), wound healing, and embryo development.
(https://medlineplus.gov/genetics/gene/fgfr3/#conditions). FGFR4 does something
similar. 

# Maybe overall story?

I think we can take the DMP_RRA results and make a story out of it. For example,
in ACC we get GWAS in all (p = 0.045954046), driven by body (p=0.02997003) and promoter1
(p=0.046953047). In cellular ocmponents, DMP_all_RRA: euchromatin, and DMP_promoter2_RRA I
get transcription repressor complex, exon-exon junction complex, and
transcription regulator complex. This is interesting as it relates to the
changes we're seeing in the transcriptome? Then, for molecular, we get
promoter2_DMP_RRA: transcription corepressor activity. And all_DMP_RRA: steroid hormone receptor activity, ligand-activated
transcription factor activity, SNAP receptor activity, repressing transcription
factor binding. Still in ACC, for KEGG all_DMP_RRA has SNARE interactions in
vesicular transport (q = .0507), but it does make a parallel to the all_DMP_RRA
finding in molecular function ontology. Not much in the reactome though. Only 
Homo sapiens: Nuclear Receptor transcription pathway and Homo sapiens:
Interferon Signaling for all_DMP_RRA.

Turning to the Caudate, there is GWAS at 0.038961039 for promoter1_DMP_RRA.
Cellular componets gets all_DMP_RRA: non-motile cilium. And the Reactome has all
that stuff about FGFR3 and 4.

There isn't anything in a single probe resolution, and my on-going efforts to
boost the FDr aren't doing anything. The per-region analysis also didn't come up
with anything.

If we go for topVar we get more stuff, and possible a different/stronger story.
If anything, we have a few single probe hits, and possibly even for the whole
region. GSEA results are very different too. But maybe we don't need to rely on
that.

# 2021-02-26 14:50:10

Let's do the same subdivisions we just did for ACC, but now for Caudate:

```r
library(minfi)
load('~/data/methylation_post_mortem/filt_Caudate_02222021.RData')
mVals <- getM(mSetSqFlt)

library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# get the table of results for the first contrast (naive - rTreg)
ann450kSub <- ann450k[match(rownames(mVals),ann450k$Name),
                      c(1:4,12:19,24:ncol(ann450k))]

load('~/data/methylation_post_mortem/res_Caudate_02222021.RData')
idx = ann450kSub$Enhancer == "TRUE"
res_cau[['enhancer']] = run_methyl(mVals[idx,], samples, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="Body") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon"))
res_cau[['body']] = run_methyl(mVals[idx,], samples, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200"))
res_cau[['promoter1']] = run_methyl(mVals[idx,], samples, 'all', ann450kSub[idx,])
idx = (grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS1500") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="TSS200") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="1stExon") |
       grepl(x=ann450kSub$UCSC_RefGene_Group, pattern="5\'UTR") )
res_cau[['promoter2']] = run_methyl(mVals[idx,], samples, 'all', ann450kSub[idx,])
save(res_cau, file='~/data/methylation_post_mortem/res_Caudate_02262021.RData')
```

Still, nothing of significance.

# 2021-03-05 10:10:18

Let's add some annotations to the GO sets results:

```r
library(WebGestaltR)
mydir = '~/data/methylation_post_mortem/'
DBs = c('Biological_Process', 'Cellular_Component', 'Molecular_Function')
for (db in DBs) {
    enrich_db = sprintf('geneontology_%s_noRedundant', db)
    a = loadGeneSet(enrichDatabase = enrich_db)
    for (r in c('ACC', 'Caudate')) {
        for (v in c('DMP', 'topVar')) {
            for (g in c('all', 'body', 'promoter1', 'promoter2')) {
                for (y in c('glm', 'RRA')) {
                    fname = sprintf('%s_%s_%s_%s_%s', r, g, v, y, db)
                    cat(fname, '\n')
                    df = read.csv(paste0(mydir, fname, '.csv'))
                    df$GO = sapply(df$ID,
                                   function(x) { a = strsplit(x, '/')[[1]];
                                                 a[length(a)] })
                    df2 = merge(df, a$geneSetDes, sort=F, by.x='GO',
                                by.y='geneSet')
                    first_cols = c('description', 'padj')
                    others = setdiff(colnames(df2), first_cols)
                    df2 = df2[, c(first_cols, others)]
                    write.csv(df2, row.names=F,
                              file=paste0(mydir, fname, '_annot.csv'))
                }
            }
        }
    }
}
```