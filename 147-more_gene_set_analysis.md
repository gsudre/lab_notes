# 2020-11-16 06:31:02

I now have a bit more results with different imputation targets, and the RNAseq
with control as the reference, so let's use those results in our summaries. It's
also a good check to see if everything ran to completion.

I won't include the methylation results for now though, even though they should
be ready. This way our summaries are a bit easier to read.

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

Doing something similar, but keeping the nominally significant hits for the
imputed data:

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
thresh = .05
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
                good_res = res[which(res$PValue < thresh), ]
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
write.csv(all_res, file='~/data/post_mortem/camera2_goodSets_FDRp05ImpNominal.csv',
          row.names=F)
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

## GSEA

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
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

for f in `cut -d"," -f 1 phenos.txt`; do nfiles=`ls WG/*${f}*v2.txt 2>/dev/null | wc -l`; echo $f $nfiles; done

/bin/bash
mydir=~/data/expression_impute;
for r in `cat ~/data/expression_impute/phenos.txt`; do
    for db in 'geneontology_Biological_Process_noRedundant' \
                'geneontology_Cellular_Component_noRedundant' \
                'geneontology_Molecular_Function_noRedundant' \
                'pathway_KEGG' 'disease_Disgenet' \
                'phenotype_Human_Phenotype_Ontology' \
                'network_PPI_BIOGRID' 'disorders' 'adhd_genes'; do
        for md in 'MASHR' 'EN'; do
            for sc in 'effect' 'zscore'; do
                if [ ! -e $mydir/WG2_${md}_${sc}_${r}_${db}_10K.csv ]; then
                    echo $r $md $sc $r $db;
                fi;
            done;
        done;
    done;
done

/bin/bash
mydir=~/data/expression_impute;
db='phenotype_Human_Phenotype_Ontology';
for r in `cat ~/data/expression_impute/phenos.txt`; do
    for md in 'MASHR' 'EN'; do
        for sc in 'effect' 'zscore'; do
            if [ ! -e $mydir/WG2_${md}_${sc}_${r}_${db}_10K.csv ]; then
                echo $r $md $sc $r $db;
            fi;
        done;
    done;
done
```

Still re-running some of the disGenet runs... but we can go ahead without it for
now just so we can have an overall idea of the results.

I realized the issue was that DisGenetd had a description with commans, and I
didn't switch that to ;. It takes too long to re-run it, so I'll change future
code (or just let quotes in!), and replace it by quotes in the files now. That
wasn't enough though... I'll likely have to re-run it. For now, I'll just do
this analysis without it.

Actually, I can just use the enrichment_results file for everything.

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'pathway_KEGG', 'disease_Disgenet',
            'network_PPI_BIOGRID',
            'phenotype_Human_Phenotype_Ontology')
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]

thresh = .05
all_res = c()
for (db in gs) {
    cat(db, '\n')
    db_res = c()
    mydir = '~/data/rnaseq_derek'
    res_var = c()
    for (r in c('acc', 'caudate')) {
        fname = sprintf('%s/WG/enrichment_results_%s_%s.txt', mydir, r, db)
        res = read.delim(fname)
        db_res = c(db_res, sum(res$FDR < thresh, na.rm=T))
        eval(parse(text=sprintf('res_var = c(res_var, "res_rnaseq_%s")', r)))
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                r = gsub(x=r, pattern='\\.', replacement='_')
                fname = sprintf('%s/WG/enrichment_results_%s_%s_%s_%s_v2.txt',
                                mydir, md, sc, r, db)
                res = read.delim(fname)
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
write.csv(t(all_res), file='~/data/post_mortem/WG_geneset_summary_FDRp05.csv')
```

The disorders, developmental, and adhd_genes sets will be done separately
because I don't think they'd survive FDR:

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'pathway_KEGG',
            'phenotype_Human_Phenotype_Ontology',
            'network_PPI_BIOGRID')
thresh = .05
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        fname = sprintf('%s/WG/enrichment_results_%s_%s.txt', mydir, r, db)
        res = read.delim(fname)
        # PPI has no descriptions
        if (db == 'network_PPI_BIOGRID') {
            res$description = res$geneSet
        }
        good_res = res[which(res$FDR < thresh), ]
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$description,
                              pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                r = gsub(x=r, pattern='\\.', replacement='_')
                fname = sprintf('%s/WG/enrichment_results_%s_%s_%s_%s_v2.txt',
                                mydir, md, sc, r, db)
                res = read.delim(fname)
                if (db == 'network_PPI_BIOGRID') {
                    res$description = res$geneSet
                }
                good_res = res[which(res$FDR < thresh), ]
                nres = nrow(good_res)
                this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                    setFamily = rep(db, nres),
                                    set = good_res$description,
                                    pval = good_res$pValue, FDR = good_res$FDR)
                all_res = rbind(all_res, this_res)
            }
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/WG_goodSets_FDRp05.csv',
          row.names=F)

# all disorder results (no FDR cut off)
gs = c('disorders', 'adhd_genes')
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        fname = sprintf('%s/WG/enrichment_results_%s_%s.txt', mydir, r, db)
        good_res = read.delim(fname)
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$link,
                              pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                r = gsub(x=r, pattern='\\.', replacement='_')
                fname = sprintf('%s/WG/enrichment_results_%s_%s_%s_%s_v2.txt',
                                mydir, md, sc, r, db)
                good_res = read.delim(fname)
                nres = nrow(good_res)
                this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                    setFamily = rep(db, nres),
                                    set = good_res$link,
                                    pval = good_res$pValue, FDR = good_res$FDR)
                all_res = rbind(all_res, this_res)
            }
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/WG_disorders_goodSets.csv',
          row.names=F)


# all developmental results (no FDR cutoff)
all_res = data.frame(set=c(), origin=c(), pval=c(), FDR=c()) 
mydir = '~/data/rnaseq_derek'
for (r in c('acc', 'caudate')) {
    fname = sprintf('%s/WG/enrichment_results_%s_%s_developmental.txt', mydir,
                    r, r)
    good_res = read.delim(fname)
    nres = nrow(good_res)
    this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                            set = good_res$link,
                            pval = good_res$pValue, FDR = good_res$FDR)
    all_res = rbind(all_res, this_res)
}
mydir = '~/data/expression_impute'
db = 'developmental'
for (md in c('MASHR', 'EN')) {
    for (sc in c('effect', 'zscore')) {
        for (r in phenotypes) {
            if (grepl(x=r, pattern='Caudate') || grepl(x=r, pattern='ATR')) {
                region = 'caudate'
            } else {
                region = 'ACC'
            }
            r = gsub(x=r, pattern='\\.', replacement='_')
            fname = sprintf('%s/WG/enrichment_results_%s_%s_%s_%s_%s_v2.txt',
                            mydir, md, sc, r, region, db)
            good_res = read.delim(fname)
            nres = nrow(good_res)
            this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                set = good_res$link,
                                pval = good_res$pValue, FDR = good_res$FDR)
            all_res = rbind(all_res, this_res)
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/WG_dev_goodSets.csv',
          row.names=F, quote=F)
```

Looking at the ADHD gene sets, we have 72 results for EN and 81 for MASHR, but
the Demonstis results are all MASHR. A bit more leaning towards zscore though.
In the dev results, I have 33 hits for EN, and only 3 for MASHR (all zscore). 
If we take th GWAS union, then 5 EN (2 effect), and 6 MASHR (all zscore). It
feels like it's best to use the EN zscore results, because then I get some
results for TWAS and CNV as well.

Doing something similar, but keeping the nominally significant hits for the
imputed data:

```r
gs = c('geneontology_Biological_Process_noRedundant',
            'geneontology_Cellular_Component_noRedundant',
            'geneontology_Molecular_Function_noRedundant',
            'pathway_KEGG',
            'phenotype_Human_Phenotype_Ontology',
            'network_PPI_BIOGRID')
thresh = .05
phenotypes = read.table('~/data/expression_impute/phenos.txt')[,1]
all_res = data.frame(set=c(), origin=c(), setFamily=c(), pval=c(), FDR=c()) 
for (db in gs) {
    mydir = '~/data/rnaseq_derek'
    for (r in c('acc', 'caudate')) {
        fname = sprintf('%s/WG/enrichment_results_%s_%s.txt', mydir, r, db)
        res = read.delim(fname)
        # PPI has no descriptions
        if (db == 'network_PPI_BIOGRID') {
            res$description = res$geneSet
        }
        good_res = res[which(res$FDR < thresh), ]
        nres = nrow(good_res)
        this_res = data.frame(origin = rep(sprintf('rnaseq_%s', r), nres),
                              setFamily = rep(db, nres),
                              set = good_res$description,
                              pval = good_res$pValue, FDR = good_res$FDR)
        all_res = rbind(all_res, this_res)
    }
    mydir = '~/data/expression_impute'
    for (md in c('MASHR', 'EN')) {
        for (sc in c('effect', 'zscore')) {
            for (r in phenotypes) {
                r = gsub(x=r, pattern='\\.', replacement='_')
                fname = sprintf('%s/WG/enrichment_results_%s_%s_%s_%s_v2.txt',
                                mydir, md, sc, r, db)
                res = read.delim(fname)
                if (db == 'network_PPI_BIOGRID') {
                    res$description = res$geneSet
                }
                good_res = res[which(res$pValue < thresh), ]
                nres = nrow(good_res)
                this_res = data.frame(origin = rep(sprintf('imp_%s_%s_%s', md, sc, r), nres),
                                    setFamily = rep(db, nres),
                                    set = good_res$description,
                                    pval = good_res$pValue, FDR = good_res$FDR)
                all_res = rbind(all_res, this_res)
            }
        }
    }
}
write.csv(all_res, file='~/data/post_mortem/WG_goodSets_FDRp05ImpNominal.csv',
          row.names=F)
```

# 2020-11-20 14:16:18

Let's see what Meff would be for our results:
How does camera (and GSEA) deal with signs in ranking? Let's run some quick tests:

```r
myregion = 'Caudate'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)
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

lcpm <- cpm(genes, log=TRUE)


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

# 2020-11-19 20:01:58

So, here's my current rationale for picking the phenotypes: we start with
rnaseq_acc and rnaseq_caudate. They have good results in both camera and WG in
developmental sets, adhd_genes, and ontologies. However, I cannot come up with
an imputed phenotype that works in the developmental sets and also in the
adhd_genes set (using GWAS1: Desmontis as the anchor). That was true for camera
and WG. So, if there will be no overlap between developmental and adhd_genes
results, we can just check within ADHD genes. Using camera first, rh_acc_volume
is a good one, with GWAS1, TWAS1, and a very recent CNV. The non-normalized
result is stronger. And that's using MASHR zscore. 

Looking at the WG results,
imp_MASHR_zscore_res_rh_caudalanteriorcingulate_volume is also the strongest
result in adhd_genes, with similar hits. I think this one is chosen then. And it
locks MASHR and zscore.

We still need a Caudate phenotype, as well as DTI phenotypes, although those are
not necessarily crucial. We could go with imp_MASHR_zscore_res_rd_ATR, which has
the same CNV hit, as well as Elia's, in camera. The ad counterpart also has the
2020 CNV hit. imp_MASHR_zscore_res_Caudate_volume has a hit with the other TWAS,
so we'd need to include it too, which is not so bad as they were both in 2019?
imp_MASHR_zscore_res_ad_cin_cin_r is related to one EWAS, and
imp_MASHR_zscore_res_rd_cin_cin_r to another. How about the WG results?

imp_MASHR_zscore_res_rd_ATR is related to both Elia's and 2020 CNV, and
imp_MASHR_zscore_res_ad_ATR with the 2020 CNV. imp_MASHR_zscore_res_rd_cin_cin
is related to GWAS1, and imp_MASHR_zscore_res_Caudate_volume with one of the
TWAS. Nothing for imp_MASHR_zscore_res_ad_cin_cin though.

At least the results seem to be quite consistent between camera and WG.

Finally, let's look at the ontology sets. I saw that not much survives FDR on
them, but maybe I can hit rnaseq results with FDR and imputation nominally?

I don't think there is a clear winner whether we need to go with camera or WG.
If anything, they somewhat agree, which is nice. The FDR to nominal ontology
results are a bit more interesting for WG, but they could just be spurious. The
developmental results make a bit more sense on camera. 

A few things I should still try:

 * Sign in ranks?
 * Overrepresentation again, maybe using Meff
 * dx*gene*region interaction
 * Up/down segmented analysis
 * Check back with David
 * Splicing factors?
 * lncRNA: “gene” sets? 
 * do pathway analysis similar to Derek's on a specific set of genes
 * maybe run camera estimating gene correlation?
 * try this approach: https://bioinformaticsbreakdown.com/how-to-gsea/ and also
   this one: https://bioinformatics-core-shared-training.github.io/cruk-summer-school-2018/RNASeq2018/html/06_Gene_set_testing.nb.html#gsea-analysis
 * try overlap with gene set lists
 * did david check diagnostic group with the PCs? email him directly

# TODO
 * Meff?
 * top X results, and then plot the worst p-value for each top and dataset
