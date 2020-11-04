# 2020-11-03 19:38:45

Now that I have some imputation results, and new methyl results, let's go back
to the overrepresentation analysis. First, do the imputation results overlap
with the RNAseq results?

```r
data_dir = '~/data/expression_impute/'
md = 'MASHR'

clean_imp_res = function(res) {
    G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
    G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
    G_list = G_list[G_list$hgnc_symbol!='',]
    G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]

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
    return(G_list2[, c('hgnc_symbol', 'pvalue')])
}

phen = 'res_ACC_thickness'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res_acc_imp = clean_imp_res(res)
phen = 'res_Caudate_volume'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res_caudate_imp = clean_imp_res(res)
phen = 'res_fa_cin_cin'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res_cincin_imp = clean_imp_res(res)
phen = 'res_FA_cc'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res_cc_imp = clean_imp_res(res)
phen = 'res_fa_ATR'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res_atr_imp = clean_imp_res(res)
```

Now that we have hgnc for imputation and p-values, let's grab the same thing for
the RNAseq results.

```r
load('~/data/rnaseq_derek/xmodal_results_10152020.RData')
library(GeneOverlap)
t = .05
ra = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
rc = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
ia = res_acc_imp[res_acc_imp$pvalue < t, 'hgnc_symbol']
ic = res_caudate_imp[res_caudate_imp$pvalue < t, 'hgnc_symbol']
icin = res_cincin_imp[res_cincin_imp$pvalue < t, 'hgnc_symbol']
icc = res_cc_imp[res_cc_imp$pvalue < t, 'hgnc_symbol']
iatr = res_atr_imp[res_atr_imp$pvalue < t, 'hgnc_symbol']
gom.obj <- newGOM(list(rnaseq_acc=ra, rnaseq_caudate=rc, imp_acc=ia,
                       imp_caudate=ic, imp_cincin=icin, imp_cc=icc,
                       imp_atr=iatr), , spec='hg19.gene')
getMatrix(gom.obj, name='pval')
```

```
                                  rnaseq_caudate     imp_acc imp_caudate     imp_cincin
rnaseq_acc     0.0000000000000000000000006440123 0.031106850  0.04388688 0.147468000176
rnaseq_caudate 1.0000000000000000000000000000000 0.001076264  0.78507449 0.182642209380
imp_acc        1.0000000000000000000000000000000 1.000000000  0.08241346 0.000004337019
imp_caudate    1.0000000000000000000000000000000 1.000000000  1.00000000 0.079242910482
imp_cincin     1.0000000000000000000000000000000 1.000000000  1.00000000 1.000000000000
imp_cc         1.0000000000000000000000000000000 1.000000000  1.00000000 1.000000000000
                                                                                        imp_cc
rnaseq_acc     0.33190707425639187633237270347308367490768432617187500000000000000000000000000
rnaseq_caudate 0.00071905751712048396329102217805484542623162269592285156250000000000000000000
imp_acc        0.00000442340390956567524650706185740389742022671271115541458129882812500000000
imp_caudate    0.00034323415039390559521867429815245031932136043906211853027343750000000000000
imp_cincin     0.00000000000000000000000000000000000000000000000000000000000000000000004498816
imp_cc         1.00000000000000000000000000000000000000000000000000000000000000000000000000000
                                                      imp_atr
rnaseq_acc     0.00116127558973589924476432777566969889448956
rnaseq_caudate 0.00136001566857726213788382096936402376741171
imp_acc        0.00276460879936775514448266122258246468845755
imp_caudate    0.00036225369932692475438348878036265432456275
imp_cincin     0.00000000000000000000005099693153116958824231
imp_cc         0.00000000000000000000000000000000000000913403
```

There seems to be a significant overlap when I'm very generous at nominal p of
.05. If I start going down, and I can go down to .005, I get within modality
results, but not across (i.e not between imputation and rnaseq).

Let me see if the elastic net results are any better. Also, I'll include the
methylation results and their slices:

```r
md = 'EN'
phen = 'res_ACC_thickness'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res2_acc_imp = clean_imp_res(res)
phen = 'res_Caudate_volume'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res2_caudate_imp = clean_imp_res(res)
phen = 'res_fa_cin_cin'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res2_cincin_imp = clean_imp_res(res)
phen = 'res_FA_cc'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res2_cc_imp = clean_imp_res(res)
phen = 'res_fa_ATR'
res = read.table(sprintf('%s/assoc_%s_%s.txt', data_dir, md, phen), header=1)
res2_atr_imp = clean_imp_res(res)
```

```r
data_dir = '~/data/methylation_post_mortem/'
for (m in c("acc", "caudate")){
    fname = sprintf('%s/%s_methyl_results_11032020.rds', data_dir, m)
    df <- readRDS(fname)
    #code from Gustavo's rank list example:
    idx = data$gene != ''
    tmp = data[idx, c('gene', 'P.Value')]
    tmp2 = c()
    for (g in unique(tmp$gene)) {
        gene_data = tmp[tmp$gene==g, ]
        best_res = which.min(abs(gene_data$P.Value))
        tmp2 = rbind(tmp2, gene_data[best_res, ])
    }
    # just to conform with imp and rna structures
    colnames(tmp2) = c('hgnc_symbol', 'pvalue')
    eval(parse(text=sprintf('res_%s_methyl = tmp2', m)))

    for (r in c("island", "opensea", "shelf", "shore")) {
        cat(m, r, '\n')
        #subset df based on cgi levels
        data <- df[df$cgi==r, ]
    
        #code from Gustavo's rank list example:
        idx = data$gene != ''
        tmp = data[idx, c('gene', 'P.Value')]
        tmp2 = c()
        for (g in unique(tmp$gene)) {
            gene_data = tmp[tmp$gene==g, ]
            best_res = which.min(abs(gene_data$P.Value))
            tmp2 = rbind(tmp2, gene_data[best_res, ])
        }
        # just to conform with imp and rna structures
        colnames(tmp2) = c('hgnc_symbol', 'pvalue')
        eval(parse(text=sprintf('res_%s_%s_methyl = tmp2', m, r)))
    }
}
```

# TODO
 * try EN results as well
 * try our own universe