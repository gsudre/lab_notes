# 2020-12-10 15:10:46

Let's see if our main DGE PM results change if we categorize them between
lncRNA, protein coding, and pseudogenes, like the Science paper:

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)
# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]

# grabing the annotations from the isoform data
df = read.delim('~/data/isoforms/shaw_adhd.rsem_output.tpm.tsv')
a = lapply(df[,1], function(x) strsplit(as.character(x), split="\\|"))
meta_iso = t(data.frame(a))
colnames(meta_iso) = c('id1', 'ensembleID', 'id2', 'id3', 'iso_name',
                        'hgnc_symbol','id4', 'read_type')
meta_iso = meta_iso[, c('ensembleID', 'hgnc_symbol', 'read_type')]
meta_genes = as.data.frame(meta_iso[!duplicated(meta_iso[, 'ensembleID']), ])
```

```
r$> all(sort(rownames(count_matrix)) == sort(meta_genes$ensembleID))                          
[1] TRUE
```

OK, so we can just sort both and keep going with the annotations.

```r
count_matrix = count_matrix[sort(rownames(count_matrix)), ]
meta_genes = meta_genes[sort(meta_genes$ensembleID, index.return=T)$ix, ]
rownames(meta_genes) = NULL

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

library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)
```

Note that I'm not restricting to autossome only, and I'm doing some general
cleaning before PCA. I'll also do PCA using everything, just so I don't have to
use different sets of PCs for different read types.

```r
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
G_list2 = meta_genes[meta_genes$ensembleID %in% rownames(geneCounts), ]
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])

num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data[, x], mydata[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}
categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(mydata[, y], data[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}
```

```
r$> which(num_pvals < .01, arr.ind = T)                                                                                       
                        row col
clusters                  2   1
RINe                      4   1
RINe                      4   2
PMI                       5   2
clusters                  2   4
pcnt_optical_duplicates   1   7
clusters                  2   7
Age                       3   8
Age                       3   9
C4                        9   9

r$> which(categ_pvals < .01, arr.ind = T)                                                                                     
                row col
batch             1   1
batch             1   2
batch             1   7
MoD               3   9
substance_group   4   9
```

```r
data2 = cbind(data, mydata)
form = ~ Diagnosis + PC1 + PC2 + PC4 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_acc = topTable(fit2, coef='DiagnosisCase', number=Inf)
```

Up until here the analysis hasn't changed much. We added a new PC, we are not
keeping to autossome only, and kept some duplicate hgnc there. That affected
some stuff, but we are only going to crop the results now. But before we do it,
let's run Caudate as well:

```r
myregion = 'Caudate'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)
# at this point we have 55 samples for ACC
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
data = data[, !grepl(colnames(data), pattern='^ENS')]

# grabing the annotations from the isoform data
df = read.delim('~/data/isoforms/shaw_adhd.rsem_output.tpm.tsv')
a = lapply(df[,1], function(x) strsplit(as.character(x), split="\\|"))
meta_iso = t(data.frame(a))
colnames(meta_iso) = c('id1', 'ensembleID', 'id2', 'id3', 'iso_name',
                        'hgnc_symbol','id4', 'read_type')
meta_iso = meta_iso[, c('ensembleID', 'hgnc_symbol', 'read_type')]
meta_genes = as.data.frame(meta_iso[!duplicated(meta_iso[, 'ensembleID']), ])

count_matrix = count_matrix[sort(rownames(count_matrix)), ]
meta_genes = meta_genes[sort(meta_genes$ensembleID, index.return=T)$ix, ]
rownames(meta_genes) = NULL

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

library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)

library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
G_list2 = meta_genes[meta_genes$ensembleID %in% rownames(geneCounts), ]
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])

num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(mydata)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(data[, x], mydata[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}
categ_vars = c('batch', 'Diagnosis', 'MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex')
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(mydata[, y], data[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}
```

```
r$> which(num_pvals < .01, arr.ind = T)                                                                                       
                        row col
RINe                      4   1
pcnt_optical_duplicates   1   3
clusters                  2   6
C1                        6   6
C2                        7   6
C7                       12   6
C1                        6   8
C2                        7   8
C3                        8   8
C4                        9   8
Age                       3   9

r$> which(categ_pvals < .01, arr.ind = T)                                                                                     
         row col
batch      1   1
batch      1   3
batch      1   5
batch      1   6
MoD        3   6
MoD        3   8
POP_CODE   6   8
```

```r
data2 = cbind(data, mydata)
form = ~ Diagnosis + PC1 + PC3 + PC5 + PC6 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
rnaseq_caudate = topTable(fit2, coef='DiagnosisCase', number=Inf)

save(rnaseq_acc, rnaseq_caudate,
     file='~/data/rnaseq_derek/results_withEverything_12102020.RData')
```

Now we can check if there is anything to FDR:

```r
for (region in c('acc', 'caudate')) {
    eval(parse(text=sprintf('tmp = rnaseq_%s', region)))
    for (dt in c('protein_coding$', 'lncRNA$', 'pseudogene$')) {
        idx = grepl(x=tmp$read_type, pattern=dt)
        res = tmp[idx, ]
        tmp$FDR = p.adjust(tmp$P.Value, method='fdr')
        print(sprintf('%s, %s = %d / %d', region, dt,
                        sum(tmp$FDR<.1), nrow(res)))
    }
}
```

```
[1] "acc, protein_coding$ = 0 / 14652"
[1] "acc, lncRNA$ = 0 / 5357"
[1] "acc, pseudogene$ = 0 / 1625"
[1] "caudate, protein_coding$ = 0 / 14654"
[1] "caudate, lncRNA$ = 0 / 5376"
[1] "caudate, pseudogene$ = 0 / 1668"
```

Nothing to see here... does it change our WG and PRS results?

```r
library(WebGestaltR)
library(tidyverse)
data_dir = '~/data/rnaseq_derek/'
load(sprintf('%s/results_withEverything_12102020.RData', data_dir))
ncpu=7

for (region in c('acc', 'caudate')) {
    eval(parse(text=sprintf('tmp = rnaseq_%s', region)))
    for (dt in c('protein_coding', 'pseudogene', 'lncRNA')) {
        idx = grepl(x=tmp$read_type, pattern=sprintf('%s$', dt))
        res = tmp[idx, ]
        # selecting the best probe hit to run enrichment with
        res = res %>% group_by(hgnc_symbol) %>%
                slice_min(n=1, P.Value, with_ties=F)
        res = as.data.frame(res)

        ranks = -log(res$P.Value) * sign(res$t)
        tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
        tmp2 = tmp2[order(ranks, decreasing=T),]

        # my own GMTs
        db = sprintf('my_%s_sets', region)
        cat(region, dt, db, '\n')
        project_name = sprintf('%s_%s_%s', region, dt, db)
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabaseFile=db_file,
                                    enrichDatabaseType="genesymbol",
                                    interestGene=tmp2,
                                    outputDirectory = data_dir,
                                    interestGeneType="genesymbol",
                                    sigMethod="top", topThr=150000,
                                    minNum=3, projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000, maxNum=800))
        if (class(enrichResult) != "try-error") {
            out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                                region, dt, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }

        DBs = c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant')
        for (db in DBs) {
            cat(region, dt, db, '\n')
            project_name = sprintf('%s_%s_%s', region, dt, db)
            enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase=db,
                                        interestGene=tmp2,
                                        interestGeneType="genesymbol",
                                        sigMethod="top", topThr=150000,
                                        outputDirectory = data_dir,
                                        minNum=5, projectName=project_name,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, perNum=10000, maxNum=800))
            if (class(enrichResult) != "try-error") {
                out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                                    region, dt, db)
                write.csv(enrichResult, file=out_fname, row.names=F)
            }
        }
    }
}
```

And for PRS we need to re-run the code above to get data2 and genes, and then:

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/results_withEverything_12102020.RData')

fname = '~/data/post_mortem/genotyping/1KG/merged_PM_1KG_PRS_12032020.csv'
prs = read.csv(fname)
prs$hbcc_brain_id = sapply(prs$IID,
                          function(x) {
                              br = strsplit(x, '_')[[1]][2];
                              as.numeric(gsub(br, pattern='BR',
                                              replacement=''))})
imWNH = data$C1 > 0 & data$C2 < -.075
wnh_brains = data[which(imWNH),]$hbcc_brain_id
# using the most appropriate PRS
m = merge(data, prs, by='hbcc_brain_id')
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
m[, prs_names] = NA
keep_me = m$hbcc_brain_id %in% wnh_brains
m[keep_me, prs_names] = m[keep_me, 64:75]
m[!keep_me, prs_names] = m[!keep_me, 52:63]
data_app = m[, c(1:51, 76:87)]

rownames(mydata) = data$hbcc_brain_id
data2 = merge(data_app, mydata, by.x='hbcc_brain_id', by.y=0)
genes2 = genes[, data$hbcc_brain_id %in% data2$hbcc_brain_id]

all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    if (myregion=='ACC') {
        form = as.formula(sprintf('~ %s + PC1 + PC2 + PC4 + PC7 + PC8 + PC9', p))
        res_dx = rnaseq_acc
        out_fname = '~/data/post_mortem/allSplit_acc_prs_overlap_results.csv'
    } else {
        form = as.formula(sprintf('~ %s + PC1 + PC3 + PC5 + PC6 + PC8 + PC9', p))
        res_dx = rnaseq_caudate
        out_fname = '~/data/post_mortem/allSplit_caudate_prs_overlap_results.csv'
    }
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res_prs = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        for (dt in c('protein_coding', 'lncRNA', 'pseudogene')) {
            idx = grepl(x=res_prs$read_type, pattern=sprintf('%s$', dt))
            res_prs2 = res_prs[idx, ]
            idx = grepl(x=res_dx$read_type, pattern=sprintf('%s$', dt))
            res_dx2 = res_dx[idx, ]

            prs_genes = res_prs2[res_prs2$P.Value < t, 'hgnc_symbol']
            dx_genes = res_dx2[res_dx2$P.Value < t, 'hgnc_symbol']
            go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res_prs))
            go.obj <- testGeneOverlap(go.obj)
            inter = intersect(prs_genes, dx_genes)
            pval = getPval(go.obj)
            this_res = c(dt, p, t, length(prs_genes), length(dx_genes), length(inter),
                        pval)
            all_res = rbind(all_res, this_res)
        }
    }
}
colnames(all_res) = c('read_type', 'PRS', 'nomPvalThresh', 'PRsgenes', 'PMgenes',
                      'overlap', 'pval')
write.csv(all_res, file=out_fname, row.names=F)
```



