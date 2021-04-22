# 2021-04-22 11:16:58

Philip suggested I should re-run the robustness and check other metrics,
basically seeing if our actula results hold up. And also aplit the covariates as
follows:

```
[10:23 AM] Shaw, Philip (NIH/NHGRI) [E]
    I found the correlations for robustness a little hard to interpret----they will almost inevitably be very highly correlated---- but does adjustment changed the main results (which woudl be the N passing IHW; the results of GSEA; the correlations).   
​[10:24 AM] Shaw, Philip (NIH/NHGRI) [E]
    It wouldn't be feasible to run these adjusting for coviarates one at a time.  
​[10:24 AM] Shaw, Philip (NIH/NHGRI) [E]
    But perhaps for 'blocks' of covairates it might. 
​[10:24 AM] Shaw, Philip (NIH/NHGRI) [E]
    I bet that a reviewer will ask for the fully adjusted....lol
​[10:25 AM] Shaw, Philip (NIH/NHGRI) [E]
    SO the presentaiton of robustness would become a table-----along the top are heading for the main findings.....the ACC 'hits'; teh GSEA (brief summary); and the sig correlations between disorders.   
​[10:25 AM] Shaw, Philip (NIH/NHGRI) [E]
    the rows would be 
​[10:25 AM] Shaw, Philip (NIH/NHGRI) [E]
    (1) the unadjsuted
​[10:25 AM] Shaw, Philip (NIH/NHGRI) [E]
    (2) adjusted for demo
​[10:25 AM] Shaw, Philip (NIH/NHGRI) [E]
    (3) adjusted for clinical
​[10:25 AM] Shaw, Philip (NIH/NHGRI) [E]
    (4) adjusted for technical
(5) adjusted for all (or at least those passing 0,05)
```

Philip also sent a table in a Word document in Teams as a template of how the
reporting for that should be. I'm thinking we can add overlaps to that too.

```r
basic_DGE = function(myregion, add_cov=NA) {
    data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
    rownames(data) = data[,1]
    data[,1] = NULL
    data = round(data)
    sub_name = gsub(x=colnames(data), pattern='X', replacement='')
    colnames(data) = sub_name
    # this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
    # its other replicate
    data = data[, ! colnames(data) %in% c('66552')]
    # outliers based on PCA plots
    outliers = c('68080','68096', '68108', '68084', '68082')
    data = data[, ! colnames(data) %in% outliers]

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
    df$BBB = factor(sapply(1:nrow(df),
                            function(x) sprintf('%s_%s',
                                        as.character(df[x,'brainbank']),
                                        as.character(df[x, 'batch']))))
    df$BBB2 = NA                                                                        
    df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
    df[df$batch==3, 'BBB2'] = 2                                                         
    df[df$batch==4, 'BBB2'] = 3      
    df$BBB2 = factor(df$BBB2)                                                   
                        
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
    if (is.na(add_cov)) {
        fm_str = '~ BBB2 + Diagnosis'
    } else {
        cov_str = paste0(add_cov, collapse = ' + ')
        fm_str = sprintf('~ %s + BBB2 + Diagnosis', cov_str)
        # making sure any numeric covariates are scaled
        num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
                'C1', 'C2', 'C3', 'C4', 'C5')
        for (var in num_vars) {
            df[, var] = scale(df[, var])
        }
    }
    cat('Running', fm_str)
    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = df,
                                  design = as.formula(fm_str))

    keep_me = colData(dds)$Region == myregion
    dds = dds[, keep_me]
    design = model.matrix(as.formula(fm_str), data=colData(dds))
    dds = DESeq(dds)

    library(edgeR)
    isexpr <- filterByExpr(counts(dds), design=design)
    ddsExpr = dds[isexpr, ]
    ddsExpr = DESeq(ddsExpr)

    return(ddsExpr)
}
```

Now, let's run the usual covariates, then the grouped version Philip suggested:

```r
myregion = 'ACC'
covs = c('clusters', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI',
              'comorbid_group', 'pcnt_optical_duplicates',
              'C4', 'C5', 'MoD', 'substance_group', 'evidence_level')
adjusted = list()
for (mycov in covs) {
    cat(mycov, '\n')
    dds = basic_DGE(myregion, add_cov=mycov)
    adjusted[[mycov]] = dds
}

mycov = c('Age', 'Sex', 'C1', 'C2', 'C3')
dds = basic_DGE(myregion, add_cov=mycov)
adjusted[['demo3']] = dds

mycov = c('Age', 'Sex', 'C1', 'C2', 'C3', 'C4', 'C5')
dds = basic_DGE(myregion, add_cov=mycov)
adjusted[['demo5']] = dds

mycov = c('PMI', 'comorbid_group', 'substance_group', 'MoD', 'evidence_level')
dds = basic_DGE(myregion, add_cov=mycov)
adjusted[['clin']] = dds

mycov = c('clusters', 'RINe', 'pcnt_optical_duplicates')
dds = basic_DGE(myregion, add_cov=mycov)
adjusted[['tech']] = dds

# splitting it into smaller files
covs = c('clusters', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI',
         'comorbid_group', 'pcnt_optical_duplicates',
         'C4', 'C5', 'MoD', 'substance_group', 'evidence_level',
         'demo3', 'demo5', 'clin', 'tech')
for (mycov in covs) {
    cat(mycov, '\n')
    out_fname = sprintf('~/data/post_mortem/cov_dds_BBB2_04222021_%s_%s.RData',
                    myregion, mycov)
    dds = adjusted[[mycov]]
    save(dds, file=out_fname)
}
```

We will also need to code running GSEA for each of them. I'll just do
developmental first, then I'll do the other DBs:

```r
library(WebGestaltR)
library(DESeq2)

data_dir = '~/data/post_mortem/'
ncpu=31
region = 'ACC'

covs = c('clusters', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI',
         'comorbid_group', 'pcnt_optical_duplicates',
         'C4', 'C5', 'MoD', 'substance_group', 'evidence_level',
         'demo3', 'demo5', 'clin', 'tech')

for (mycov in covs) {
    load(sprintf('~/data/post_mortem/cov_dds_BBB2_04222021_%s_%s.RData',
                 region, mycov))
    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))
    
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)
    
    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    res_str = sprintf('ROB1_DGE_%s_%s_BBB2', region, mycov)
    DBs = c(sprintf('%s_developmental', tolower(region)))
    for (db in DBs) {
        cat(res_str, db, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        project_name = sprintf('%s_%s_10K', res_str, db)
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
}
```

And the same thing for the other DBs:

```r
for (mycov in covs) {
    load(sprintf('~/data/post_mortem/cov_dds_BBB2_04222021_%s_%s.RData',
                 region, mycov))
    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))
    
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)
    
    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    res_str = sprintf('ROB1_DGE_%s_%s_BBB2', region, mycov)
    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(res_str, db, '\n')
        project_name = sprintf('%s_%s_10K', res_str, db)

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
}
```

OK, so let's generate something akin to that table Philip suggested

```r
covs = c('clusters', 'Age', 'Sex', 'C1', 'C2', 'C3', 'RINe', 'PMI',
         'comorbid_group', 'pcnt_optical_duplicates',
         'C4', 'C5', 'MoD', 'substance_group', 'evidence_level',
         'demo3', 'demo5', 'clin', 'tech')
covs = c('clusters', 'Age')

library(DESeq2)
load('~/data/post_mortem/basic_DGE_04202021.RData')
res.acc = as.data.frame(results(dds.ACC, name = "Diagnosis_Case_vs_Control",
                                alpha=.05))
library(IHW)
acc.IHW = adj_pvalues(ihw(pvalue ~ baseMean,  data=res.acc, alpha=0.05))

library(GeneOverlap)
resmat = data.frame()
# store all the metrics in a list
metrics = list()
ACC_IHW_p05 = rownames(res.acc)[which(acc.IHW < .05)]
metrics$ACC_IHW_p05 = length(ACC_IHW_p05)
ACC_FDR_p05 = rownames(res.acc)[which(res.acc$padj < .05)]
metrics$ACC_FDR_p05 = length(ACC_FDR_p05)
ACC_p00001 = rownames(res.acc)[which(res.acc$pvalue < .00001)]
metrics$ACC_p00001 = length(ACC_p00001)
metrics$baseIHWOverlap = NA
metrics$IHWOverlapPval = NA
metrics$baseFDROverlap = NA
metrics$FDROverlapPval = NA
metrics$baseP00001Overlap = NA
metrics$P00001OverlapPval = NA
r = 'ACC'
df = read.table(sprintf('~/data/post_mortem/Project_WG26_DGE_%s_BBB2_%s_developmental_10K/enrichment_results_WG26_DGE_%s_BBB2_%s_developmental_10K.txt',
                          r, tolower(r), r, tolower(r)),
                  header=1, sep='\t')[, 1:6]
metrics$ACC_panDev_adjPval = df[df$link=='overlap', 'FDR']
metrics$ACC_adult_adjPval = df[df$link=='adult (>19 yrs)', 'FDR']
r = 'Caudate'
df = read.table(sprintf('~/data/post_mortem/Project_WG26_DGE_%s_BBB2_%s_developmental_10K/enrichment_results_WG26_DGE_%s_BBB2_%s_developmental_10K.txt',
                          r, tolower(r), r, tolower(r)),
                  header=1, sep='\t')[, 1:6]
metrics$Caudate_infant_adjPval = df[df$link=='infant (0-2 yrs)', 'FDR']

# grab all significant gene sets for ACC, and just check later the adjusted
# pvalues for other conditions bind each covariate iteration as a row


resmat = rbind(resmat, metrics)
rownames(resmat)[nrow(resmat)] = 'base'
for (mycov in covs) {
    cat(mycov, '\n')
    load(sprintf('~/data/post_mortem/cov_dds_BBB2_04222021_ACC_%s.RData',
                 mycov))
    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))

    # store all the metrics in a list
    metrics = list()
    res.IHW = adj_pvalues(ihw(pvalue ~ baseMean,  data=res, alpha=0.05))
    res_IHW_p05 = rownames(res)[which(res.IHW < .05)]
    metrics$ACC_IHW_p05 = length(res_IHW_p05)
    res_FDR_p05 = rownames(res)[which(res$padj < .05)]
    metrics$ACC_FDR_p05 = length(res_FDR_p05)
    res_p00001 = rownames(res)[which(res$pvalue < .00001)]
    metrics$ACC_p00001 = length(ACC_p00001)

    gsize = length(intersect(rownames(res), rownames(res.acc)))
    metrics$baseIHWOverlap = length(intersect(ACC_IHW_p05, res_IHW_p05))
    go.obj <- newGeneOverlap(ACC_IHW_p05, res_IHW_p05, genome.size=gsize)
    go.obj <- testGeneOverlap(go.obj)
    metrics$IHWOverlapPval = getPval(go.obj)
    metrics$baseFDROverlap = length(intersect(ACC_FDR_p05, res_FDR_p05))
    go.obj <- newGeneOverlap(ACC_FDR_p05, res_FDR_p05, genome.size=gsize)
    go.obj <- testGeneOverlap(go.obj)
    metrics$FDROverlapPval = getPval(go.obj)
    metrics$baseP00001Overlap = length(intersect(ACC_p00001, res_p00001))
    go.obj <- newGeneOverlap(ACC_p00001, res_p00001, genome.size=gsize)
    go.obj <- testGeneOverlap(go.obj)
    metrics$P00001OverlapPval = getPval(go.obj)

    r = 'ACC'
    df = read.table(sprintf('~/data/post_mortem/Project_ROB1_DGE_%s_%s_BBB2_%s_developmental_10K/enrichment_results_ROB1_DGE_%s_%s_BBB2_%s_developmental_10K.txt',
                            r, mycov, tolower(r), r, mycov, tolower(r)),
                    header=1, sep='\t')[, 1:6]
    metrics$ACC_panDev_adjPval = df[df$link=='overlap', 'FDR']
    metrics$ACC_adult_adjPval = df[df$link=='adult (>19 yrs)', 'FDR']
    r = 'Caudate'
    df = read.table(sprintf('~/data/post_mortem/Project_ROB1_DGE_%s_%s_BBB2_%s_developmental_10K/enrichment_results_ROB1_DGE_%s_%s_BBB2_%s_developmental_10K.txt',
                            r, mycov, tolower(r), r, mycov, tolower(r)),
                    header=1, sep='\t')[, 1:6]
    metrics$Caudate_infant_adjPval = df[df$link=='infant (0-2 yrs)', 'FDR']
  
    # bind each covariate iteration as a row
    resmat = rbind(resmat, metrics)
    rownames(resmat)[nrow(resmat)] = mycov
}
```

# TODO
 * maybe always remove the same outliers as detected by the main model, instead
   of leaving it up to DESeq2 every time?
 * combine covariates as Philip mentioned
 * run SV as well  