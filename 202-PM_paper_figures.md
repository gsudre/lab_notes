# 2021-03-16 11:29:31

Starting with the plot for developmental results in DGE:

```r
keep_me = c('overlap_c0.9', 'dev5_c0.9', 'dev1_c0.9', 'dev4_c0.9', 'dev2_c0.9',
            'dev3_c0.9')
res = read.table('/Users/sudregp/data/post_mortem/final_results/DGE_enrichment/gene_ontologies/enrichment_results_WG7_dge_acc___protein_coding____my_acc_sets_10K.txt',
                sep='\t', header=1)
res = res[res$geneSet %in% keep_me, c('link', 'normalizedEnrichmentScore',
                                      'pValue')]
res$padj.BF = p.adjust(res$pValue, method='bonferroni')
res$padj.FDR = p.adjust(res$pValue, method='fdr')

dev = res
dev$Region = 'ACC'
res = read.table('/Users/sudregp/data/post_mortem/final_results/DGE_enrichment/gene_ontologies/enrichment_results_WG7_dge_cau___protein_coding____my_caudate_sets_10K.txt',
                sep='\t', header=1)
res = res[res$geneSet %in% keep_me, c('link', 'normalizedEnrichmentScore',
                                      'pValue')]
res$padj.BF = p.adjust(res$pValue, method='bonferroni')
res$padj.FDR = p.adjust(res$pValue, method='fdr')
res$Region = 'Caudate'
dev = rbind(dev, res)

df = matrix(nrow = 2, ncol = 6, dimnames=list(c('ACC', 'Caudate'),
                                              unique(dev$link)))
df = df[, c(1, 3, 5, 6, 4, 2)]
for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
        idx = dev$Region == rownames(df)[i] & dev$link == colnames(df)[j]
        df[i, j] = dev[idx, 'normalizedEnrichmentScore']
    }
}

df = df[, c(1, 3, 5, 6, 4, 2)]
for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
        idx = dev$Region == rownames(df)[i] & dev$link == colnames(df)[j]
        if (dev[idx, 'pValue'] == 0) {
            dev[idx, 'pValue'] = 1e-5
        }
        df[i, j] = (sign(dev[idx, 'normalizedEnrichmentScore']) *
                    -log(dev[idx, 'pValue']))
    }
}

library(corrplot)
mylim = max(abs(df))
quartz()
corrplot(df, is.corr=F, cl.lim=c(-mylim, mylim))
