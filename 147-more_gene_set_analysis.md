# 2020-11-16 06:31:02

I now have a bit more results with different imputation targets, and the RNAseq
with control as the reference, so let's use those results in our summaries. It's
also a good check to see if everything ran to completion.

I won't include the methylation results for now though, even though they should
be ready. This way our summaries are a bit easier to read.

## GSEA

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'pathway_KEGG', 'disease_Disgenet',
            'phenotype_Human_Phenotype_Ontology',
            'network_PPI_BIOGRID', 'disorders')
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]

thresh = .05
all_res = c()
for (db in gs) {
    db_res = c()
    mydir = '~/data/rnaseq_derek'
    res_var = c()
    for (r in c('acc', 'caudate')) {
        res = read.csv(sprintf('%s/WG2_%s_%s_10K.csv', mydir, r, db))
        db_res = c(db_res, sum(res$FDR < thresh, na.rm=T))
        eval(parse(text=sprintf('res_var = c(res_var, "res_rnaseq_%s")', r)))
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                res = read.csv(sprintf('%s/WG2_%s_%s_%s_%s_10K.csv', mydir, md, sc, r, db))
                db_res = c(db_res, sum(res$FDR < thresh, na.rm=T))
                eval(parse(text=sprintf('res_var = c(res_var, "res_%s_%s_%s")',
                                        md, sc, r)))
            }
        }
    }
    all_res = rbind(all_res, db_res)
}
colnames(all_res) = res_var
rownames(all_res) = gs
write.csv(t(all_res), file='~/data/post_mortem/webgestalt2_geneset_summary_FDRp05.csv')
```

Not everything is done yet... figuring out what's missing.

```bash
for f in `cut -d"," -f 1 phenos.txt`; do nfiles=`ls WG2*${f}*10K.csv 2>/dev/null | wc -l`; echo $f $nfiles; done

/bin/bash
mydir=~/data/expression_impute;
for r in `cat ~/data/expression_impute/phenos.txt`; do
    for db in 'geneontology_Biological_Process_noRedundant' \
                'geneontology_Cellular_Component_noRedundant' \
                'geneontology_Molecular_Function_noRedundant' \
                'pathway_KEGG' 'disease_Disgenet' \
                'phenotype_Human_Phenotype_Ontology' \
                'network_PPI_BIOGRID'; do
        for md in 'MASHR' 'EN'; do
            for sc in 'effect' 'zscore'; do
                if [ ! -e $mydir/WG2_${md}_${sc}_${r}_${db}_10K.csv ]; then
                    echo $r $md $sc $r $db;
                fi;
            done;
        done;
    done;
done
```

## camera

Doing the same for camera, which always finish much faster:

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'geneontology_Biological_Process',
            'geneontology_Cellular_Component',
            'geneontology_Molecular_Function',
            'pathway_KEGG', 'disease_Disgenet',
            'phenotype_Human_Phenotype_Ontology',
            'network_PPI_BIOGRID')
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]
thresh = .05
all_res = c()
for (db in gs) {
    db_res = c()
    mydir = '~/data/rnaseq_derek'
    res_var = c()
    for (r in c('acc', 'caudate')) {
        res = read.csv(sprintf('%s/camera2_%s_%s.csv', mydir, r, db))
        db_res = c(db_res, sum(res$FDR < thresh, na.rm=T))
        eval(parse(text=sprintf('res_var = c(res_var, "res_rnaseq_%s")', r)))
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                res = read.csv(sprintf('%s/camera2_%s_%s_%s_%s.csv', mydir, md, sc, r, db))
                db_res = c(db_res, sum(res$FDR < thresh, na.rm=T))
                eval(parse(text=sprintf('res_var = c(res_var, "res_%s_%s_%s")',
                                        md, sc, r)))
            }
        }
    }
    all_res = rbind(all_res, db_res)
}
colnames(all_res) = res_var
rownames(all_res) = gs
write.csv(t(all_res), file='~/data/post_mortem/camera2_geneset_summary_FDRp05.csv')
```

The disorders, developmental, and adhd_genes sets will be done separately
because I don't think they'd survive FDR:

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'geneontology_Biological_Process',
            'geneontology_Cellular_Component',
            'geneontology_Molecular_Function',
            'pathway_KEGG', 'disease_Disgenet',
            'phenotype_Human_Phenotype_Ontology',
            'network_PPI_BIOGRID', 'disorders')
thresh = .1
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        res = read.csv(sprintf('%s/camera2_%s_%s.csv', mydir, r, db))
        good_res = res[which(res$FDR < thresh), ]
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$description,
                              pval = good_res$PValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                res = read.csv(sprintf('%s/camera2_%s_%s_%s_%s.csv', mydir, md, sc, r, db))
                good_res = res[which(res$FDR < thresh), ]
                nres = nrow(good_res)
                this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                    setFamily = rep(db, nres),
                                    set = good_res$description,
                                    pval = good_res$PValue, FDR = good_res$FDR)
                all_res = rbind(all_res, this_res)
            }
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/camera2_geneset_goodSets_FDRp1.csv',
          row.names=F, quote=F)

# all disorder results (no FDR cut off)
gs = c('disorders', 'adhd_genes')
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        good_res = read.csv(sprintf('%s/camera2_%s_%s.csv', mydir, r, db))
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$description,
                              pval = good_res$PValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                good_res = read.csv(sprintf('%s/camera2_%s_%s_%s_%s.csv', mydir, md, sc, r, db))
                nres = nrow(good_res)
                this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                    setFamily = rep(db, nres),
                                    set = good_res$description,
                                    pval = good_res$PValue, FDR = good_res$FDR)
                all_res = rbind(all_res, this_res)
            }
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/camera2_disorders_goodSets.csv',
          row.names=F, quote=F)


# all developmental results (no FDR cutoff)
all_res = data.frame(set=c(), origin=c(), pval=c(), FDR=c()) 
mydir = '~/data/rnaseq_derek'
for (r in c('acc', 'caudate')) {
    good_res = read.csv(sprintf('%s/camera2_%s_%s_developmental.csv', mydir, r, r))
    nres = nrow(good_res)
    this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                            set = good_res$description,
                            pval = good_res$PValue, FDR = good_res$FDR)
    all_res = rbind(all_res, this_res)
}
mydir = '~/data/expression_impute'
for (md in c('MASHR', 'EN')) {
    for (sc in c('effect', 'zscore')) {
        for (r in phenotypes) {
            if (grepl(x=r, pattern='Caudate') || grepl(x=r, pattern='ATR')) {
                region = 'Caudate'
            } else {
                region = 'ACC'
            }
            good_res = read.csv(sprintf('%s/camera2_%s_%s_%s_%s_developmental.csv',
                                        mydir, md, sc, r, region))
            nres = nrow(good_res)
            this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                set = good_res$description,
                                pval = good_res$PValue, FDR = good_res$FDR)
            all_res = rbind(all_res, this_res)
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/camera2_dev_goodSets.csv',
          row.names=F, quote=F)
```

# 2020-11-17 06:21:23

Once again, the effect results were much better than the zscore results using
camera, if we look at the big gene sets. Difference between MASHR and EN models
is a bit more subtle. Might have to do it based on gene set overlap. On the
other hand, the zscore results had more hits related to ADHD gene sets than the
effect results, so it's not out of the question yet. EN had more than MASHR by
about 90-50 ratio. If we take the Demonstis genes as a benchmark, then MASHR and
EN do equally well, but they are all zscore. TWAS1 leans towards EN, and all
zscore too. In the dev results, all hits were for EN, heavily leaning towards
effect, although zscore was the only one related to rd_cin_cin and infancy set.

MASHR effect was also the only intersection between rnaseq_acc and the
imputation results, for "serotonin receptor signaling pathway". 
So, these results are not very conclusive. For example, if I choose effect I
lose one set of results, and choosing zscore I lose the other. Maybe the GSEA
results will shed some light?


```r
cc = cor(data)
M = nrow(cc)
cnt = 0
for (j in 1:M) {
    print(j)
    for (k in 1:M) {
        cnt = cnt + (1 - cc[j, k]**2)
    }
}
meff = 1 + cnt / M
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

But that assumes samples as rows, which is not our case. And we cannot put a
square matrix of 650K in memory. So, maybe we can change this to compute on the
fly?


```r
# calculates Meff without computing the costly big cc matrix, but paying in run
# time to calculate each correlation in the loop.
# mydata is vars by samples
slow_meff = function(mydata) {
    M = nrow(mydata)
    cnt = 0
    for (j in 1:M) {
        # print(j)
        for (k in 1:M) {
            cnt = cnt + (1 - cor(mydata[j, ], mydata[k, ])**2)
        }
    }
    meff = 1 + cnt / M
    cat(sprintf('Galwey Meff = %.2f\n', meff))
    return(meff)
}
```

# TODO
 * Meff?
 * top X results, and then plot the worst p-value for each top and dataset
