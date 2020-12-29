# 2020-12-28 10:23:57

Let's repeat the WGCNA results but this time evaluating how robust those
clusters actually are:

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
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]
# we're down from 60K to 38K samples by only looking at the ones with hgnc symbol. We might be losing too much here, so it's a step to reconsider in the future

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
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]
library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
genes = calcNormFactors( genes)

library(WGCNA)
datExpr0 = t(cpm(genes, log=T))
rm_vars = data[,c('pcnt_optical_duplicates', 'clusters', 'RINe', 'batch')]
rm_vars$batch = as.numeric(rm_vars$batch)
datExpr1 = empiricalBayesLM(datExpr0, rm_vars, verbose=2) 
datExpr = datExpr1$adjustedData

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

enableWGCNAThreads()
net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       maxBlockSize=nGenes,
                       verbose = 3)
save(net, datExpr, data,
     file='~/data/WGCNA/pmACC_unsigned_pw6_mms30_mch25.RData')
```

And let's check if we're getting the same values as before:

```r
library(WGCNA)
load('~/data/WGCNA/pmACC_unsigned_pw6_mms30_mch25.RData')
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];

MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

keep_me = !is.na(data$C1)
mydata = cbind(MEs, data)[keep_me,]
res_MEs = MEs[keep_me,]
library(MASS)
myps = c()
for (v in colnames(MEs)) {
    fm_str = sprintf('%s ~ Diagnosis + Age + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 +
                     C10 + MoD + substance_group + comorbid_group + POP_CODE +
                     Sex', v)
    res.lm <- lm(as.formula(fm_str), data = mydata)
    step <- stepAIC(res.lm, direction = "both", trace = F,
                    scope = list(lower = as.formula('~ Diagnosis')))
    myp = summary(step)$coefficients['DiagnosisCase', 'Pr(>|t|)']
    cat(sprintf('%s', v), 'p =', myp, '\n')
    myps = c(myps, myp)
}
names(myps) = colnames(MEs)
```

So, there's some stuff there. I'm just not sure about a few things:

What modules should we look into?

```
r$> myps[myps<.05]                                 
  MEroyalblue MElightyellow       MEblack 
  0.020849906   0.023761104   0.002744299 
```

And according to this I can just run the target data and reference data using
the same input, to calculate robustness:

https://support.bioconductor.org/p/115857/

Peter is one of the authors, so we're safe here too.

```r
multiExpr = multiData(Set1 = datExpr, Set2 = datExpr)
colorList = list(Set1 = moduleColors)
enableWGCNAThreads()
mp = modulePreservation(multiExpr, colorList, referenceNetworks=1,
                          nPermutations = 10,
                          networkType = "unsigned",
                          randomSeed = 42,
                          quickCor=0, corFnc = "cor", parallelCalculation = T,
                          verbose = 4, 
                          indent = 0)
```

So, this code works. But as Philip pointed out, the current results don't add
much. Can we get a few modules that map to different GO sets? That would
contribute a lot more to the story. Otherwise, this network analysis will likely
not go in, even if it's robust.

So, let's see what other networks we can create here:

# 2020-12-29 14:39:06

This actually might be a good idea if we can find complementary results to the
main PM results, but also if these results correlate better to PRS?

```r
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft_bs = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                           corFnc = 'bicor', networkType = 'signed')
sft_bu = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                           corFnc = 'bicor', networkType = 'unsigned')
sft_bh = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                           corFnc = 'bicor', networkType = 'signed hybrid')
source('~/tmp/code1.R')                               
library(doParallel)
sft_bh2 = hpickSoftThreshold(datExpr, powerVector = powers, corFnc = bicor,
                             networkType = "hybrid2", verbose = 5)
```

Let's compare some of those results:

```r
quartz()
par(mfrow = c(2,2));
cex1 = 0.9;

plot(sft_bs$fitIndices[,1], -sign(sft_bs$fitIndices[,3])*sft_bs$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = 'signed');
text(sft_bs$fitIndices[,1], -sign(sft_bs$fitIndices[,3])*sft_bs$fitIndices[,2],
    labels=powers,cex=cex1,col="red");

plot(sft_bu$fitIndices[,1], -sign(sft_bu$fitIndices[,3])*sft_bu$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = 'unsigned');
text(sft_bu$fitIndices[,1], -sign(sft_bu$fitIndices[,3])*sft_bu$fitIndices[,2],
    labels=powers,cex=cex1,col="red");

plot(sft_bh$fitIndices[,1], -sign(sft_bh$fitIndices[,3])*sft_bh$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = 'hybrid');
text(sft_bh$fitIndices[,1], -sign(sft_bh$fitIndices[,3])*sft_bh$fitIndices[,2],
    labels=powers,cex=cex1,col="red");

plot(sft_bh2$fitIndices[,1], -sign(sft_bh2$fitIndices[,3])*sft_bh2$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = 'hybrid2');
text(sft_bh2$fitIndices[,1], -sign(sft_bh2$fitIndices[,3])*sft_bh2$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
```

![](images/2020-12-29-14-54-17.png)

We clearly have a lot of work to do here. Apparently different network
configurations might give different results. Let's go slowly then:

```r
nets = blockwiseModules(datExpr, power = 9, corType='bicor',
                        TOMType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = F, maxBlockSize=nGenes, verbose = 3)

netu = blockwiseModules(datExpr, power = 5, corType='bicor',
                        TOMType = "unsigned", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = F, maxBlockSize=nGenes, verbose = 3)

neth = blockwiseModules(datExpr, power = 4, corType='bicor',
                        TOMType = "signed hybrid", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = F, maxBlockSize=nGenes, verbose = 3)

adj = Hadjacency(datExpr=datExpr, type='hybrid2', power=12, corFnc='bicor')
neth2 = blockwiseModules(adj, power = 4, corType='bicor',
                        TOMType = "signed hybrid", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = F, maxBlockSize=nGenes, verbose = 3)
```

# TODO
 * how stable are these networks? (http://pages.stat.wisc.edu/~yandell/statgen/ucla/WGCNA/wgcna.html)
 * check FDR after stability analysis
 * use signed networks? (that's what Science paper did, and removed all      covariates first too)
 * use robustness (bicor)
 * try csuWGCNA (https://github.com/RujiaDai/csuWGCNA, like in Science paper)
 * is bicor the best thing to use here?