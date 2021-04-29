# 2021-04-29 11:04:29

WE saw that DX was dependent on comorbidity and substance abuse. But, for WNH,
we removed everyone else but WNH. Why don't we try the same approach for
comorbidity and substance, and see how it goes?

```r
pca_DGE_clean = function(myregion, fm_str, varvec) {
    # varvec is a list of variable name and the group to keep
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
    if (!is.na(varvec)) {
        for (v in 1:length(varvec)) {
            keep_me = which(df[, names(varvec)[v]] == varvec[v])
            df = df[keep_me, ]
        }
    }

    data = data[, colnames(data) %in% df$submitted_name]
    df = df[df$submitted_name %in% colnames(data), ]
    df = df[order(df$submitted_name), ]
    data = data[, order(df$submitted_name)]

    keep_me = df$Region == myregion
    data = data[, keep_me]
    df = df[keep_me, ]

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
    imWNH = which(df$C1 > 0 & df$C2 < -.075)
    df$POP_BIN = 'other'
    df[imWNH, 'POP_BIN'] = 'WNH'
    df$POP_BIN = factor(df$POP_BIN)        
    # bining so DESeq2 can do its own filyering automatically
    # breaks = quantile(df$RINe, probs = seq(0, 1, by = 0.25))
    # df$RINc = cut(df$RINe, breaks=breaks, labels=c('q1', 'q2', 'q3', 'q4'),
    #             include.lowest=T)
    df$RINc = cut(df$RINe, breaks = 5, include.lowest=T)  

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

    # remove constant genes (including zeros) as it breaks PCA
    const_genes = apply(data, 1, sd) == 0
    data = data[!const_genes, ]

    library("DESeq2")
    # making sure any numeric covariates are scaled
    num_vars = c('pcnt_optical_duplicates', 'clusters', 'Age', 'RINe', 'PMI',
            'C1', 'C2', 'C3', 'C4', 'C5')
    for (var in num_vars) {
        df[, var] = scale(df[, var])
    }

    cat('Running', fm_str, '\n')
    dds <- DESeqDataSetFromMatrix(countData = data,
                                  colData = df,
                                  design = as.formula(fm_str))

    min_subjs = min(table(df$Diagnosis))
    keep <- rowSums(counts(dds) == 0) <= min_subjs
    dds <- dds[keep,]
    dds = DESeq(dds)

    library(edgeR)
    design = model.matrix(as.formula(fm_str), data=colData(dds))
    isexpr <- filterByExpr(counts(dds), design=design)
    ddsExpr = dds[isexpr, ]
    ddsExpr = DESeq(ddsExpr)

    return(ddsExpr)
}
```

OK, let's do the analysis then. I'm also keeping the old dates to make the
scripts run easier:

```r
varvec = vector(mode='numeric')
varvec[1] = 0
names(varvec) = 'substance_group'
dds.ACC = pca_DGE_clean('ACC', '~ RINe + BBB2 + Diagnosis', varvec)
dds.Caudate = pca_DGE_clean('Caudate', '~ RINe + BBB2 + Diagnosis', varvec)
save(dds.ACC, dds.Caudate,
     file='~/data/post_mortem/pca_DGE_RINe_substanceClean_04262021.RData')

varvec = vector(mode='character')
varvec[1] = 'no'
names(varvec) = 'comorbid_group_update'
dds.ACC = pca_DGE_clean('ACC', '~ RINe + BBB2 + Diagnosis', varvec)
dds.Caudate = pca_DGE_clean('Caudate', '~ RINe + BBB2 + Diagnosis', varvec)
save(dds.ACC, dds.Caudate,
     file='~/data/post_mortem/pca_DGE_RINe_comorbidClean_04262021.RData')
```

And we need to run the rest of the posthoc analysis for these two new sets:

```r
library(WebGestaltR)
library(DESeq2)

data_dir = '~/data/post_mortem/'
ncpu=31
region = 'ACC'

covs = c('comorbidClean', 'substanceClean')

for (mycov in covs) {
    load(sprintf('~/data/post_mortem/pca_DGE_RINe_%s_04262021.RData', mycov))

    res_str = sprintf('dds = dds.%s', region)
    eval(parse(text=res_str))
    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))
    
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)
    
    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    res_str = sprintf('ROB2_DGE_%s_%s_RINe_BBB2', region, mycov)
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
                                    sigMethod="top", topThr=150000,
                                    outputDirectory = data_dir,
                                    minNum=3, projectName=project_name,
                                    isOutput=F, isParallel=T,
                                    nThreads=ncpu, perNum=10000))
        out_fname = sprintf('%s/%s.csv', data_dir, project_name)
        write.csv(enrichResult, file=out_fname, row.names=F, quote=T)
    }
    for (db in c('KEGG', 'Panther', 'Reactome', 'Wikipathway')) {
        cat(res_str, db, '\n')
        project_name = sprintf('%s_%s_10K', res_str, db)

        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=sprintf('pathway_%s', db),
                                    interestGene=tmp2,
                                    interestGeneType="ensembl_gene_id",
                                    sigMethod="top", minNum=3, 
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    isOutput=F, isParallel=T,
                                    nThreads=ncpu, topThr=150000, perNum=10000))
        out_fname = sprintf('%s/%s.csv', data_dir, project_name)
        write.csv(enrichResult, file=out_fname, row.names=F, quote=T)
    }
}
```

Now we do MAGMA:

```r
library(org.Hs.eg.db)
library(GenomicFeatures)
G_list0 = select(org.Hs.eg.db, keys(org.Hs.eg.db, 'ENSEMBL'),
                 columns=c('ENSEMBL', 'ENTREZID', 'SYMBOL'), 'ENSEMBL')
library(dplyr)
library(DESeq2)
for (g in c('comorbidClean', 'substanceClean')) {
    cat(g, '\n')
    load(sprintf('~/data/post_mortem/pca_DGE_RINe_%s_04262021.RData', g))
    for (r in c('ACC', 'Caudate')) {
        res_str = sprintf('dds = dds.%s', r)
        eval(parse(text=res_str))

        res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))

        res$GENEID = substr(rownames(res), 1, 15)
        G_list <- G_list0[!is.na(G_list0$ENSEMBL),]
        G_list = G_list[G_list$ENSEMBL!='',]
        G_list <- G_list[!duplicated(G_list$ENSEMBL),]
        imnamed = res$GENEID %in% G_list$ENSEMBL
        res = res[imnamed, ]
        res2 = merge(res, G_list, sort=F, all.x=F, all.y=F, by.x='GENEID',
                    by.y='ENSEMBL')
        ranks = res2 %>% group_by(ENTREZID) %>% slice_min(n=1, pvalue, with_ties=F)
        myres = data.frame(gene=ranks$ENTREZID,
                        signed_rank=sign(ranks$log2FoldChange)*-log(ranks$pvalue),
                        unsigned_rank=-log(ranks$pvalue))
        out_fname = sprintf('~/data/post_mortem/MAGMA_RINe_BBB2_%s_dge_%s.tab',
                            g, r)
        write.table(myres, row.names=F, sep='\t', file=out_fname, quote=F)
    }
}
```

Then, for MAGMA we only need to run the last command:

```bash
module load MAGMA
cd ~/data/tmp
for g in 'substanceClean' 'comorbidClean'; do
    for r in 'ACC' 'Caudate'; do
        magma --gene-results genes_BW.genes.raw \
            --gene-covar ~/data/post_mortem/MAGMA_RINe_BBB2_${g}_dge_${r}.tab \
            --out ~/data/post_mortem/MAGMA_RINe_BBB2_gc_${g}_dge_${r};
    done;
done
```

```
(base) [sudregp@cn3158 post_mortem]$ for f in `ls MAGMA_RINe_BBB2_gc_*Clean_dge_*.gsa.out`; do echo $f; cat $f; done
MAGMA_RINe_BBB2_gc_comorbidClean_dge_ACC.gsa.out
# MEAN_SAMPLE_SIZE = 55374
# TOTAL_GENES = 14821
# TEST_DIRECTION = one-sided, positive (set), two-sided (covar)
# CONDITIONED_INTERNAL = gene size, gene density, inverse mac, log(gene size), log(gene density), log(inverse mac)
VARIABLE           TYPE  NGENES         BETA     BETA_STD           SE            P
signed_rank       COVAR   14809     -0.01204     -0.01874    0.0051165     0.018633
unsigned_rank     COVAR   14809    -0.004466   -0.0049087    0.0071283      0.53099
MAGMA_RINe_BBB2_gc_comorbidClean_dge_Caudate.gsa.out
# MEAN_SAMPLE_SIZE = 55374
# TOTAL_GENES = 14803
# TEST_DIRECTION = one-sided, positive (set), two-sided (covar)
# CONDITIONED_INTERNAL = gene size, gene density, inverse mac, log(gene size), log(gene density), log(inverse mac)
VARIABLE           TYPE  NGENES         BETA     BETA_STD           SE            P
signed_rank       COVAR   14797    0.0014783    0.0026013    0.0046974        0.753
unsigned_rank     COVAR   14797    0.0027191     0.003304    0.0062765      0.66487
MAGMA_RINe_BBB2_gc_substanceClean_dge_ACC.gsa.out
# MEAN_SAMPLE_SIZE = 55374
# TOTAL_GENES = 14832
# TEST_DIRECTION = one-sided, positive (set), two-sided (covar)
# CONDITIONED_INTERNAL = gene size, gene density, inverse mac, log(gene size), log(gene density), log(inverse mac)
VARIABLE           TYPE  NGENES         BETA     BETA_STD           SE            P
signed_rank       COVAR   14818    -0.011074    -0.022048    0.0040788    0.0066361
unsigned_rank     COVAR   14818    -0.010543     -0.01605    0.0052946     0.046475
MAGMA_RINe_BBB2_gc_substanceClean_dge_Caudate.gsa.out
# MEAN_SAMPLE_SIZE = 55374
# TOTAL_GENES = 14827
# TEST_DIRECTION = one-sided, positive (set), two-sided (covar)
# CONDITIONED_INTERNAL = gene size, gene density, inverse mac, log(gene size), log(gene density), log(inverse mac)
VARIABLE           TYPE  NGENES         BETA     BETA_STD           SE            P
signed_rank       COVAR   14820     0.001826    0.0028977    0.0051442      0.72263
unsigned_rank     COVAR   14820   0.00047722   0.00053239    0.0070515      0.94604
```

We might as well use the same code to export the single gene results:

```r
library(DESeq2)
library(IHW)
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
mydir = '~/data/post_mortem/'

library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
               "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]

for (g in c('_substanceClean', '_comorbidClean')) {
    load(sprintf('~/data/post_mortem/pca_DGE_RINe%s_04262021.RData', g))
    for (r in c('ACC', 'Caudate')) {
        res_str = sprintf('res = results(dds.%s, name = "Diagnosis_Case_vs_Control", alpha=.05)',
                        r)
        eval(parse(text=res_str))
        fname = sprintf('%s/DGE_%s_RINe_BBB2%s_annot_04262021.csv', mydir, r, g)

        df = as.data.frame(res)
        colnames(df)[ncol(df)] = 'padj.FDR'
        df$padj.IHW = adj_pvalues(ihw(pvalue ~ baseMean,  data=df, alpha=0.05))
        df$GENEID = substr(rownames(df), 1, 15)
        df2 = merge(df, mart, sort=F,
                    by.x='GENEID', by.y='ensembl_gene_id', all.x=T, all.y=F)
        df2 = merge(df2, bt_slim, sort=F,
                    by.x='GENEID', by.y='gene_id', all.x=T, all.y=F)
        df2 = df2[order(df2$pvalue), ]
        
        write.csv(df2, row.names=F, file=fname)
    }
}
```

Looking at all dev sets:

```r
library(WebGestaltR)
library(DESeq2)

data_dir = '~/data/post_mortem/'
ncpu=7

region = 'ACC'

covs = c('comorbidClean', 'substanceClean')

for (mycov in covs) {
    load(sprintf('~/data/post_mortem/pca_DGE_RINe_%s_04262021.RData', mycov))

    res_str = sprintf('dds = dds.%s', region)
    eval(parse(text=res_str))

    res = as.data.frame(results(dds, name = "Diagnosis_Case_vs_Control"))
    
    ranks = -log(res$pvalue) * sign(res$log2FoldChange)
    geneid = substring(rownames(res), 1, 15)
    
    tmp2 = data.frame(geneid=geneid, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    res_str = sprintf('WG30_DGE_%s_%s_RINe_BBB2', region, mycov)

    DBs = c('%s_manySets_co0.900', '%s_manySets_co0.950')
    for (db in DBs) {
        db2 = sprintf(db, tolower(region))
        cat(res_str, db2, '\n')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db2)
        project_name = sprintf('%s_%s_10K', res_str, db2)
        enrichResult <- try(WebGestaltR(enrichMethod="GSEA",
                            organism="hsapiens",
                            enrichDatabaseFile=db_file,
                            enrichDatabaseType="genesymbol",
                            interestGene=tmp2,
                            outputDirectory = data_dir,
                            interestGeneType="ensembl_gene_id",
                            sigMethod="top", topThr=50,
                            minNum=3, projectName=project_name,
                            isOutput=T, isParallel=T,
                            nThreads=ncpu, perNum=10000, maxNum=2000))
    }
}
```

For correlation to other disorders, just use the code from note 215, changing
the data that's loaded and then the filename to save.


# TODO
