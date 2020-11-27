# 2020-11-25 08:57:32

Let me see how the WG results change if using -logP*sign(logFC) for
rank.

```r
# bw
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
load(sprintf('%s/rnaseq_results_11122020.rData', data_dir))
ncpu=8

# region='acc'
region='caudate'
eval(parse(text=sprintf('res = rnaseq_%s', region)))

ranks = -log(res$P.Value) * sign(res$logFC)
tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
tmp2 = tmp2[order(ranks, decreasing=T),]

# my own GMTs
for (db in c('adhd_genes', sprintf('%s_developmental', region))) {
    cat(region, db, '\n')
    project_name = sprintf('%s_%s', region, db)
    db_file = sprintf('~/data/post_mortem/%s.gmt', db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabaseFile=db_file,
                                enrichDatabaseType="genesymbol",
                                interestGene=tmp2,
                                outputDirectory = data_dir,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                minNum=3, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800)
    out_fname = sprintf('%s/WG3_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
DBs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
for (db in DBs) {
    cat(region, db, '\n')
    project_name = sprintf('%s_%s', region, db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabase=db,
                                interestGene=tmp2,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                outputDirectory = data_dir,
                                minNum=5, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800)
    out_fname = sprintf('%s/WG3_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

And do the same thing for the imputation results:

```r
# bw
library(WebGestaltR)

data_dir = '~/data/expression_impute/'
phenotypes = read.table(sprintf('%s/phenos.txt', data_dir))[,1]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
ncpu=8

for (phen in phenotypes) {
    if (grepl(x=phen, pattern='Caudate') || grepl(x=phen, pattern='ATR')) {
        region = 'caudate'
    } else {
        region = 'ACC'
    }
    for (md in c('EN', 'MASHR')) {
        res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen),
                        header=1)
        id_num = sapply(res$gene, function(x) strsplit(x=x, split='\\.')[[1]][1])
        dups = duplicated(id_num)
        id_num = id_num[!dups]
        res$id_num = id_num

        imnamed = res$id_num %in% G_list$ensembl_gene_id
        res = res[imnamed, ]
        G_list2 = merge(G_list, res, by.x='ensembl_gene_id', by.y='id_num')
        imautosome = which(G_list2$chromosome_name != 'X' &
                        G_list2$chromosome_name != 'Y' &
                        G_list2$chromosome_name != 'MT')
        G_list2 = G_list2[imautosome, ]
        tmp = G_list2[, c('hgnc_symbol', 'pvalue', 'zscore')]
        ranks = -log(tmp$pvalue) * sign(tmp$zscore)
        tmp2 = data.frame(hgnc_symbol=tmp$hgnc_symbol, rank=ranks)
        tmp2 = tmp2[order(ranks, decreasing=T),]
        
        for (db in c('adhd_genes', sprintf('%s_developmental', region))) {
            cat(md, phen, db, '\n')
            project_name = sprintf('logP_%s_%s_%s', md, phen, db)
            # make sure no dots in the name
            project_name = gsub(x=project_name, pattern='\\.',
                                replacement='_')
            db_file = sprintf('~/data/post_mortem/%s.gmt', db)
            enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabaseFile=db_file,
                                        enrichDatabaseType="genesymbol",
                                        interestGene=tmp2,
                                        interestGeneType="genesymbol",
                                        sigMethod="top", topThr=150000,
                                        minNum=3,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, perNum=10000,
                                        outputDirectory = data_dir,
                                        projectName=project_name,
                                        maxNum=800)
            out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                                md, phen, db)
            out_fname = gsub(x=out_fname, pattern='\\.',
                                replacement='_')
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
        DBs = c('geneontology_Biological_Process_noRedundant',
                'geneontology_Cellular_Component_noRedundant',
                'geneontology_Molecular_Function_noRedundant')
        for (db in DBs) {
            cat(region, db, '\n')
            project_name = sprintf('logP_%s_%s_%s', md, phen, db)
            enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                        organism="hsapiens",
                                        enrichDatabase=db,
                                        interestGene=tmp2,
                                        interestGeneType="genesymbol",
                                        sigMethod="top", topThr=150000,
                                        outputDirectory = data_dir,
                                        minNum=5, projectName=project_name,
                                        isOutput=T, isParallel=T,
                                        nThreads=ncpu, perNum=10000, maxNum=800)
            out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                                md, phen, db)
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }
}
```

And then for splices:

```r
# bw
library(WebGestaltR)
ncpu=8
data_dir = '~/data/isoforms/'

r = 'acc'
df = read.csv(sprintf('~/data/post_mortem/david_pca_%s.csv', r))
for (p in colnames(df)[2:4]) {
    tmp = df[!is.na(df[, p]), c('geneName', p)]
    colnames(tmp) = c('hgnc_symbol', 'P.Value')
    dup_genes = tmp$hgnc_symbol[duplicated(tmp$hgnc_symbol)]
    res = tmp[!tmp$hgnc_symbol %in% dup_genes, ]
    for (g in dup_genes) {
        gene_data = tmp[tmp$hgnc_symbol==g, ]
        best_res = which.min(gene_data$P.Value)
        res = rbind(res, gene_data[best_res, ])
    }
    ranks = -log(res$P.Value)
    names(ranks) = res$hgnc_symbol
    ranks = sort(ranks, decreasing=T)
    tmp2 = data.frame(hgnc_symbol=names(ranks), rank=ranks)
    for (db in c('adhd_genes', sprintf('%s_developmental', r))) {
        cat(r, p, db, '\n')
        project_name = sprintf('logPiso_%s_%s_%s', r, p, db)
        # make sure no dots in the name
        project_name = gsub(x=project_name, pattern='\\.',
                            replacement='_')
        db_file = sprintf('~/data/post_mortem/%s.gmt', db)
        enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabaseFile=db_file,
                                    enrichDatabaseType="genesymbol",
                                    interestGene=tmp2,
                                    interestGeneType="genesymbol",
                                    sigMethod="top", topThr=150000,
                                    minNum=3,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000,
                                    outputDirectory = data_dir,
                                    projectName=project_name,
                                    maxNum=800)
        out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                            r, p, db)
        out_fname = gsub(x=out_fname, pattern='\\.',
                            replacement='_')
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
    DBs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant')
    for (db in DBs) {
        cat(r, db, '\n')
        project_name = sprintf('logPiso_%s_%s_%s', r, p, db)
        enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                    organism="hsapiens",
                                    enrichDatabase=db,
                                    interestGene=tmp2,
                                    interestGeneType="genesymbol",
                                    sigMethod="top", topThr=150000,
                                    outputDirectory = data_dir,
                                    minNum=5, projectName=project_name,
                                    isOutput=T, isParallel=T,
                                    nThreads=ncpu, perNum=10000, maxNum=800)
        out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                            r, p, db)
        write.csv(enrichResult, file=out_fname, row.names=F)
    }
}
```

Now we need to compile the WG results again. But let's grab the enrichment
results first, something like:

```bash
cd ~/data/isoforms/WG3
scp helix:~/data/isoforms/Project*/enrichment_results*txt .
```

```r
gs = c('geneontology_Biological_Process_noRedundant',
        'geneontology_Cellular_Component_noRedundant',
        'geneontology_Molecular_Function_noRedundant')
thresh = .05
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        fname = sprintf('%s/WG3/enrichment_results_%s_%s.txt', mydir, r, db)
        res = read.delim(fname)
        # PPI has no descriptions
        if (db == 'network_PPI_BIOGRID') {
            res$description = res$geneSet
        }
        good_res = res[which(res$FDR < thresh), ]
        nres = nrow(good_res)
        good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                    'Upregulated', 'Downregulated')
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$description,
                              direction = good_res$direction,
                              pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
            for (r in phenotypes) {
                r = gsub(x=r, pattern='\\.', replacement='_')
                fname = sprintf('%s/WG3/enrichment_results_logP_%s_%s_%s.txt',
                                mydir, md, r, db)
                res = read.delim(fname)
                if (db == 'network_PPI_BIOGRID') {
                    res$description = res$geneSet
                }
                good_res = res[which(res$pValue < thresh), ]
                good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                            'Upregulated', 'Downregulated')
                nres = nrow(good_res)
                this_res = data.frame(origin = rep(sprintf('imp_%s_%s', md, r), nres),
                                    setFamily = rep(db, nres),
                                    set = good_res$description,
                                    direction = good_res$direction,
                                    pval = good_res$pValue, FDR = good_res$FDR)
                all_res = rbind(all_res, this_res)
        }
    }
    mydir = '~/data/isoforms'
    for (r in c('acc', 'caudate')) {
        MDs = ifelse(r=='acc', c('p_val_Kallisto', 'p_val_RSEM', 'GM_of_both'),
                     c('Kallisto_p_val', 'RSEM_p_val', 'GM_of_both'))
        for (md in MDs) {
            fname = sprintf('%s/WG3/enrichment_results_logPiso_%s_%s_%s.txt',
                            mydir, r, md, db)
            res = read.delim(fname)
            if (db == 'network_PPI_BIOGRID') {
                res$description = res$geneSet
            }
            good_res = res[which(res$pValue < thresh), ]
            good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                        'Upregulated', 'Downregulated')
            nres = nrow(good_res)
            this_res = data.frame(origin = rep(sprintf('iso_%s_%s', r, md), nres),
                                setFamily = rep(db, nres),
                                set = good_res$description,
                                direction = good_res$direction,
                                pval = good_res$pValue, FDR = good_res$FDR)
            all_res = rbind(all_res, this_res)
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/WG3_goodSets_FDRp05ImpNominal.csv',
          row.names=F)
```

And now for our own GMTs:

```r
# all disorder results (no FDR cut off)
gs = c('adhd_genes')
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        fname = sprintf('%s/WG3/enrichment_results_%s_%s.txt', mydir, r, db)
        good_res = read.delim(fname)
        good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                    'Upregulated', 'Downregulated')
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$link,
                              direction = good_res$direction,
                              pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (r in phenotypes) {
            r = gsub(x=r, pattern='\\.', replacement='_')
            fname = sprintf('%s/WG3/enrichment_results_logP_%s_%s_%s.txt',
                            mydir, md, r, db)
            good_res = read.delim(fname)
            good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                        'Upregulated', 'Downregulated')
            nres = nrow(good_res)
            this_res = data.frame(origin = rep(sprintf('imp_%s_%s', md, r), nres),
                                setFamily = rep(db, nres),
                                set = good_res$link,
                                direction = good_res$direction,
                                pval = good_res$pValue, FDR = good_res$FDR)
            all_res = rbind(all_res, this_res)
        }
    }
    mydir = '~/data/isoforms'
    for (r in c('acc', 'caudate')) {
        MDs = ifelse(r=='acc', c('p_val_Kallisto', 'p_val_RSEM', 'GM_of_both'),
                     c('Kallisto_p_val', 'RSEM_p_val', 'GM_of_both'))
        for (md in MDs) {
            fname = sprintf('%s/WG3/enrichment_results_logPiso_%s_%s_%s.txt',
                            mydir, r, md, db)
            good_res = read.delim(fname)
            good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                        'Upregulated', 'Downregulated')
            nres = nrow(good_res)
            this_res = data.frame(origin = rep(sprintf('iso_%s_%s', r, md), nres),
                                setFamily = rep(db, nres),
                                set = good_res$link,
                                direction = good_res$direction,
                                pval = good_res$pValue, FDR = good_res$FDR)
            all_res = rbind(all_res, this_res)
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/WG3_disorders_goodSets.csv',
          row.names=F)


# all developmental results (no FDR cutoff)
all_res = data.frame(set=c(), origin=c(), pval=c(), FDR=c()) 
mydir = '~/data/rnaseq_derek'
for (r in c('acc', 'caudate')) {
    fname = sprintf('%s/WG3/enrichment_results_%s_%s_developmental.txt', mydir,
                    r, r)
    good_res = read.delim(fname)
    good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                'Upregulated', 'Downregulated')
    nres = nrow(good_res)
    this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                            set = good_res$link,
                            direction = good_res$direction,
                            pval = good_res$pValue, FDR = good_res$FDR)
    all_res = rbind(all_res, this_res)
}
mydir = '~/data/expression_impute'
db = 'developmental'
for (md in c('MASHR', 'EN')) {
    for (r in phenotypes) {
        if (grepl(x=r, pattern='Caudate') || grepl(x=r, pattern='ATR')) {
            region = 'caudate'
        } else {
            region = 'ACC'
        }
        r = gsub(x=r, pattern='\\.', replacement='_')
        fname = sprintf('%s/WG3/enrichment_results_logP_%s_%s_%s_%s.txt',
                        mydir, md, r, region, db)
        good_res = read.delim(fname)
        good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                    'Upregulated', 'Downregulated')
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('imp_%s_%s', md, r), nres),
                            set = good_res$link,
                            direction = good_res$direction,
                            pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
}
mydir = '~/data/isoforms'
for (r in c('acc', 'caudate')) {
    MDs = ifelse(r=='acc', c('p_val_Kallisto', 'p_val_RSEM', 'GM_of_both'),
                     c('Kallisto_p_val', 'RSEM_p_val', 'GM_of_both'))
    for (md in MDs) {
        fname = sprintf('%s/WG3/enrichment_results_logPiso_%s_%s_%s_%s.txt',
                        mydir, r, md, r, db)
        good_res = read.delim(fname)
        good_res$direction = ifelse(good_res$normalizedEnrichmentScore > 0,
                                    'Upregulated', 'Downregulated')
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('iso_%s_%s', r, md), nres),
                            set = good_res$link,
                            direction = good_res$direction,
                            pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
}
write.csv(all_res, file='~/data/post_mortem/WG3_dev_goodSets.csv',
          row.names=F, quote=F)
```

# 2020-11-27 16:16:56

I'm going to run just the sets we're interested in to get FDR:

```r
library(ActivePathways)
gmt = read.GMT('~/data/post_mortem/acc_developmental.gmt')
gmt2 = read.GMT('~/data/post_mortem/adhd_genes.gmt')
gmt[['GWAS1']] = gmt2[['GWAS1']]
gmt[['TWAS']] = gmt2[['TWAS']]
write.GMT(gmt, '~/data/post_mortem/my_acc_sets.gmt')

gmt = read.GMT('~/data/post_mortem/caudate_developmental.gmt')
gmt2 = read.GMT('~/data/post_mortem/adhd_genes.gmt')
gmt[['GWAS1']] = gmt2[['GWAS1']]
gmt[['TWAS']] = gmt2[['TWAS']]
write.GMT(gmt, '~/data/post_mortem/my_caudate_sets.gmt')
```

```r
library(WebGestaltR)

data_dir = '~/data/rnaseq_derek/'
load(sprintf('%s/rnaseq_results_11122020.rData', data_dir))
ncpu=2

for (region in c('acc', 'caudate')) {
    eval(parse(text=sprintf('res = rnaseq_%s', region)))

    ranks = -log(res$P.Value) * sign(res$logFC)
    tmp2 = data.frame(hgnc_symbol=res$hgnc_symbol, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]

    # my own GMTs
    db = sprintf('my_%s_sets', region)
    cat(region, db, '\n')
    project_name = sprintf('%s_%s', region, db)
    db_file = sprintf('~/data/post_mortem/%s.gmt', db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabaseFile=db_file,
                                enrichDatabaseType="genesymbol",
                                interestGene=tmp2,
                                outputDirectory = data_dir,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                minNum=3, projectName=project_name,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000, maxNum=800)
    out_fname = sprintf('%s/WG3_%s_%s_10K.csv', data_dir, region, db)
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

Just a note that I tried fixing the seed before calling WebGestaltR, and even
using rg_rnd within the call, but the results kept changing. So, let's just go
with what we have. I can't think of another way to do it other than changing the
functions within the package.

And repeat it for imputations:

```r
library(WebGestaltR)

data_dir = '~/data/expression_impute/'
phenotypes = c('res_rh_caudalanteriorcingulate_volume',
               'res_norm_ACC_vol', 'res_ACC_volume',
               'res_norm_Caudate_vol', 'res_Caudate_volume')

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
ncpu=2

for (phen in phenotypes) {
    if (grepl(x=phen, pattern='Caudate') || grepl(x=phen, pattern='ATR')) {
        region = 'caudate'
    } else {
        region = 'ACC'
    }
    res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen),
                    header=1)
    id_num = sapply(res$gene, function(x) strsplit(x=x, split='\\.')[[1]][1])
    dups = duplicated(id_num)
    id_num = id_num[!dups]
    res$id_num = id_num

    imnamed = res$id_num %in% G_list$ensembl_gene_id
    res = res[imnamed, ]
    G_list2 = merge(G_list, res, by.x='ensembl_gene_id', by.y='id_num')
    imautosome = which(G_list2$chromosome_name != 'X' &
                    G_list2$chromosome_name != 'Y' &
                    G_list2$chromosome_name != 'MT')
    G_list2 = G_list2[imautosome, ]
    tmp = G_list2[, c('hgnc_symbol', 'pvalue', 'zscore')]
    ranks = -log(tmp$pvalue) * sign(tmp$zscore)
    tmp2 = data.frame(hgnc_symbol=tmp$hgnc_symbol, rank=ranks)
    tmp2 = tmp2[order(ranks, decreasing=T),]
    
    db = sprintf('my_%s_sets', region)
    cat(md, phen, db, '\n')
    project_name = sprintf('logP_%s_%s_%s', md, phen, db)
    # make sure no dots in the name
    project_name = gsub(x=project_name, pattern='\\.',
                        replacement='_')
    db_file = sprintf('~/data/post_mortem/%s.gmt', db)
    enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                organism="hsapiens",
                                enrichDatabaseFile=db_file,
                                enrichDatabaseType="genesymbol",
                                interestGene=tmp2,
                                interestGeneType="genesymbol",
                                sigMethod="top", topThr=150000,
                                minNum=3,
                                isOutput=T, isParallel=T,
                                nThreads=ncpu, perNum=10000,
                                outputDirectory = data_dir,
                                projectName=project_name,
                                maxNum=800)
    out_fname = sprintf('%s/WG3_%s_%s_%s_10K.csv', data_dir,
                        md, phen, db)
    out_fname = gsub(x=out_fname, pattern='\\.',
                        replacement='_')
    write.csv(enrichResult, file=out_fname, row.names=F)
}
```

