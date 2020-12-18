# 2020-11-12 07:07:55

Afyer talking to Philip I'll re-run more options for imputations phenotypes, and
then check if they agree with ADHD clinical phenotypes. Worst case scenario we
just rely on literature. I'll start with DTI because those are easier, and I can
prepare the anatomical later as the DTI runs.

```r
gf = read.csv('~/data/expression_impute/gfWithDTIandClin_602_10292020.csv')
data = gf[gf$bestInFamily=='TRUE',]

keep_idx = data$goodVolumes >= 56 & data$norm.trans <= 1.5 & data$norm.rot <= .035
data = data[keep_idx,]

a = read.table('~/data/expression_impute/results/NCR_v3_ACC_predict_1KG_en.txt', header=1)
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
imp_data = merge(imp_data, data, by.x='IID', by.y='Subject.Code...Subjects',
                 all.x=F, all.y=F)
imwnh = imp_data$PC01<0 & imp_data$PC02>-.02
brain_vars = colnames(data)[c(29:88, 104:136)]
nogene_idx = !grepl(x=colnames(imp_data), pattern='^ENS')

data2 = imp_data[imwnh, nogene_idx]
for (m in c('ad', 'fa', 'rd')) {
    data2[, sprintf('%s_ATR', m)] = data2[, sprintf('%s_ATR_l', m)] + data2[, sprintf('%s_ATR_r', m)]
    data2[, sprintf('%s_cin_cin', m)] = data2[, sprintf('%s_cin_cin_l', m)] + data2[, sprintf('%s_cin_cin_r', m)]
    brain_vars = c(brain_vars, sprintf('%s_ATR', m), sprintf('%s_cin_cin', m))
}
for (v in brain_vars) {
    m = mean(data2[, v], na.rm=T)
    s = sd(data2[, v], na.rm=T)
    data2[which(data2[, v] > m + 3*s), v] = NA
    data2[which(data2[, v] < m - 3*s), v] = NA
}
```

Down to 233 subjects now. Residualizing:

```r
library(MASS)
for (v in brain_vars) {
    fm_str = sprintf('%s ~ Sex...Subjects + age_acq + scanner_update + sequence_type + norm.trans + norm.rot + goodVolumes', v)
    res.lm <- lm(as.formula(fm_str), data = data2, na.action=na.exclude)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    data2[, sprintf('res_%s', v)] = scale(residuals(step))
}
```

Let's now create the files to be used for association:

```r
a = read.table('~/data/expression_impute/results/NCR_v3_ACC_predict_1KG_en.txt', header=1)
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
imp_data = merge(imp_data, data2, by='IID', all.x=F, all.y=F)

imwnh = imp_data[imp_data$PC01.x<0 & imp_data$PC02.x>-.02,]$IID
data_dir = '~/data/expression_impute/'
phenotypes = list(ACC=c('res_fa_cin_cin', 'res_ad_cin_cin', 'res_rd_cin_cin',
                        'res_fa_cin_cin_l', 'res_ad_cin_cin_l',
                        'res_rd_cin_cin_l',
                        'res_fa_cin_cin_r', 'res_ad_cin_cin_r',
                        'res_rd_cin_cin_r'),
                  Caudate=c('res_fa_ATR', 'res_ad_ATR', 'res_rd_ATR',
                            'res_fa_ATR_l', 'res_ad_ATR_l', 'res_rd_ATR_l',
                            'res_fa_ATR_r', 'res_ad_ATR_r', 'res_rd_ATR_r'))
for (region in c('ACC', 'Caudate')) {
   for (my_phen in phenotypes[[region]]) {
       print(my_phen)
      data3 = data2[data2$IID %in% imwnh, ]
      data3 = data3[data3$bestInFamily==T, ]
      data3 = data3[, c('subject.id', my_phen)]
      colnames(data3)[1] = 'IID'
      colnames(data3)[2] = 'phen'
      data3 = data3[order(data3$IID), ]
      # it expects no more than the number of people we have in the phenotypes
    #   a = read.table(sprintf('%s/results/NCR_v3_%s_predict_1KG_en.txt',
    #                          data_dir, region), header=1)
       a = readRDS(sprintf('%s/results/NCR_v3_%s_1KG_mashr.rds', data_dir,
                           region))
      # remove FAMID from IID
      iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
      iid3 = gsub(x=iid2, pattern='SID.', replacement='')
      a$IID = as.numeric(iid3)
      b = a[a$IID %in% data3$IID, ]
      b = b[order(b$IID), ]
      data3$FID = b$FID # they're both sorted on IID
    #   write.table(b, file=sprintf('%s/DTI_cropped_imp_EN_%s.tab', data_dir,
    #                               region), row.names=F, quote=F, sep='\t')
      write.table(b, file=sprintf('%s/DTI_cropped_imp_MASHR_%s.tab', data_dir,
                                  region), row.names=F, quote=F, sep='\t')
      write.table(data3, file=sprintf('%s/phen_%s.tab', data_dir, my_phen),
                  row.names=F, quote=F, sep='\t')
   }
}
```

And we run the associations:

```bash
# laptop
source /Users/sudregp/opt/miniconda3/etc/profile.d/conda.sh
conda activate imlabtools
DATA=~/data/expression_impute;
METAXCAN=~/data/expression_impute/MetaXcan/software;
for m in fa ad rd; do
    for h in '' '_l' '_r'; do
        phen=res_${m}_cin_cin${h};
        python3 $METAXCAN/PrediXcanAssociation.py \
                --expression_file $DATA/DTI_cropped_imp_MASHR_ACC.tab \
            --input_phenos_file $DATA/phen_${phen}.tab \
            --input_phenos_column phen \
            --output $DATA/assoc_MASHR_${phen}.txt \
            --verbosity 9;
    done;
done
for m in fa ad rd; do
    for h in '' '_l' '_r'; do
        phen=res_${m}_ATR${h};
        python3 $METAXCAN/PrediXcanAssociation.py \
                --expression_file $DATA/DTI_cropped_imp_MASHR_Caudate.tab \
            --input_phenos_file $DATA/phen_${phen}.tab \
            --input_phenos_column phen \
            --output $DATA/assoc_MASHR_${phen}.txt \
            --verbosity 9;
    done;
done
```

Then switch the above to run the MASHR datasets, and it's just a matter of
running the gene set analysis in BW. Before I do that, let's create the anatomy
datasets too:

```r
gf = read.csv('~/data/expression_impute/gfWithMPRAGE_632_11102020.csv')
data = gf[gf$bestInFamily=='TRUE',]
keep_idx = data$QC...Scan <= 2.5 & data$external_score <=2.5 & data$internal_score <= 2.5
data = data[keep_idx,]

a = read.table('~/data/expression_impute/results/NCR_v3_ACC_predict_1KG_en.txt', header=1)
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
imp_data = merge(imp_data, data, by.x='IID', by.y='Subject.Code...Subjects',
                 all.x=F, all.y=F)
imwnh = imp_data$PC01<0 & imp_data$PC02>-.02
nogene_idx = !grepl(x=colnames(imp_data), pattern='^ENS')
data2 = imp_data[imwnh, nogene_idx]
brain_vars = colnames(data)[c(24:35)]

for (v in brain_vars) {
    m = mean(data2[, v], na.rm=T)
    s = sd(data2[, v], na.rm=T)
    data2[which(data2[, v] > m + 3*s), v] = NA
    data2[which(data2[, v] < m - 3*s), v] = NA
}
data2$norm_lh_ACC_vol = data2$lh_caudalanteriorcingulate_volume / data2$EstimatedTotalIntraCranialVol
data2$norm_rh_ACC_vol = data2$rh_caudalanteriorcingulate_volume / data2$EstimatedTotalIntraCranialVol
data2$norm_ACC_vol = data2$ACC_volume / data2$EstimatedTotalIntraCranialVol
data2$norm_lh_Caudate_vol = data2$Left.Caudate / data2$EstimatedTotalIntraCranialVol
data2$norm_rh_Caudate_vol = data2$Right.Caudate / data2$EstimatedTotalIntraCranialVol
data2$norm_Caudate_vol = data2$Caudate_volume / data2$EstimatedTotalIntraCranialVol
brain_vars = c(brain_vars, c('norm_lh_ACC_vol', 'norm_rh_ACC_vol',
                             'norm_ACC_vol', 'norm_lh_Caudate_vol','norm_rh_Caudate_vol', 'norm_Caudate_vol'))
data2$MeanThickness = data2$lh_MeanThickness_thickness + data2$rh_MeanThickness_thickness
data2$norm_lh_ACC_thi = data2$lh_caudalanteriorcingulate_thickness / data2$MeanThickness
data2$norm_rh_ACC_thi = data2$rh_caudalanteriorcingulate_thickness / data2$MeanThickness
data2$norm_ACC_thi = data2$ACC_thickness / data2$MeanThickness
brain_vars = c(brain_vars, c('norm_lh_ACC_thi', 'norm_rh_ACC_thi',
                             'norm_ACC_thi'))

library(nlme)
library(MASS)
for (v in brain_vars) {
    fm_str = sprintf('%s ~ Sex...Subjects + scanner_update + age_scan + QC...Scan + external_score + internal_score', v)
    fit <- lme(as.formula(fm_str), random=~1|FAMID, data = data2, na.action=na.exclude, method='ML')
    step <- stepAIC(fit, direction = "both", trace = F)
    data2[, sprintf('res_%s', v)] = scale(residuals(step))
}
```

And we run those phenotypes through TWAS. 

```r
a = read.table('~/data/expression_impute/results/NCR_v3_ACC_predict_1KG_en.txt', header=1)
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
imp_data = merge(imp_data, data2, by='IID', all.x=F, all.y=F)

imwnh = imp_data[imp_data$PC01.x<0 & imp_data$PC02.x>-.02,]$IID
data_dir = '~/data/expression_impute/'
for (bv in brain_vars) {
    my_phen = sprintf('res_%s', bv)
    if (grepl(x=bv, pattern='Caudate')) {
        region = 'Caudate'
    } else {
        region = 'ACC'
    }
    cat(my_phen, region, '\n')
    data3 = data2[data2$IID %in% imwnh, ]
    data3 = data3[data3$bestInFamily==T, ]
    data3 = data3[, c('IID', my_phen)]
    colnames(data3)[2] = 'phen'
    data3 = data3[order(data3$IID), ]
    # it expects no more than the number of people we have in the phenotypes
    a = read.table(sprintf('%s/results/NCR_v3_%s_predict_1KG_en.txt',
                            data_dir, region), header=1)
    # a = readRDS(sprintf('%s/results/NCR_v3_%s_1KG_mashr.rds', data_dir,
    #                     region))
    # remove FAMID from IID
    iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
    a$IID = iid2
    b = a[a$IID %in% data3$IID, ]
    b = b[order(b$IID), ]
    data3$FID = b$FID # they're both sorted on IID
    write.table(b, file=sprintf('%s/ANAT_cropped_imp_EN_%s.tab', data_dir,
                                region), row.names=F, quote=F, sep='\t')
    # write.table(b, file=sprintf('%s/ANAT_cropped_imp_MASHR_%s.tab', data_dir,
    #                             region), row.names=F, quote=F, sep='\t')
    write.table(data3, file=sprintf('%s/phen_%s.tab', data_dir, my_phen),
                row.names=F, quote=F, sep='\t')
}
```

Now we run the associations for these phenotypes as well:

```bash
# laptop
source /Users/sudregp/opt/miniconda3/etc/profile.d/conda.sh
conda activate imlabtools
DATA=~/data/expression_impute;
METAXCAN=~/data/expression_impute/MetaXcan/software;
for phen in res_lh_caudalanteriorcingulate_area \
    res_rh_caudalanteriorcingulate_area \
    res_lh_caudalanteriorcingulate_volume \
    res_rh_caudalanteriorcingulate_volume \
    res_lh_caudalanteriorcingulate_thickness \
    res_rh_caudalanteriorcingulate_thickness \
    res_ACC_area res_ACC_volume res_ACC_thickness \
    res_norm_lh_ACC_vol res_norm_rh_ACC_vol res_norm_ACC_vol \
    res_norm_lh_ACC_thi res_norm_rh_ACC_thi res_norm_ACC_thi; do
        python3 $METAXCAN/PrediXcanAssociation.py \
                --expression_file $DATA/ANAT_cropped_imp_EN_ACC.tab \
            --input_phenos_file $DATA/phen_${phen}.tab \
            --input_phenos_column phen \
            --output $DATA/assoc_EN_${phen}.txt \
            --verbosity 9;
done
for phen in res_Left.Caudate res_Right.Caudate res_Caudate_volume \
    res_norm_lh_Caudate_vol res_norm_rh_Caudate_vol res_norm_Caudate_vol; do
        python3 $METAXCAN/PrediXcanAssociation.py \
                --expression_file $DATA/ANAT_cropped_imp_EN_Caudate.tab \
            --input_phenos_file $DATA/phen_${phen}.tab \
            --input_phenos_column phen \
            --output $DATA/assoc_EN_${phen}.txt \
            --verbosity 9;
done
```

Now we get ready to run camera and GSEA in both sets of phenotypes together.

Let's create the new GMT with Gauri's list of ADHD genes:

```r
library(ActivePathways)
library(gdata)
df = read.xls('~/data/adhd_gene_list_gwas_twas_ewas_cnv_11112020_GS.xlsx',
             'to_import')
df = df[!is.na(df$geneList_code2), ]
nsets = length(unique(df$geneList_code2)) + length(unique(df$geneList_code1))
gmt = read.GMT('~/data/post_mortem/hsapiens_disease_Disgenet_entrezgene.gmt')
junk = gmt[seq(1, nsets)]
cnt = 1
for (i in unique(df$geneList_code2)) {
    myset = df[df$geneList_code2==i, ]
    desc = sprintf('%s; %s', i, myset[1, 'Author.year'])
    genes = unique(gsub(x=myset$Gene, pattern=' ', replacement=''))
    mygs = list(id=i, name=desc, genes=genes)
    junk[[cnt]] = mygs
    cnt = cnt + 1
}
for (i in unique(df$geneList_code1)) {
    myset = df[df$geneList_code1==i, ]
    genes = unique(gsub(x=myset$Gene, pattern=' ', replacement=''))
    desc = sprintf('%s; union', i)
    mygs = list(id=i, name=desc, genes=genes)
    junk[[cnt]] = mygs
    cnt = cnt + 1
}
write.GMT(junk, '~/data/post_mortem/adhd_genes.gmt')
```

Let's start with camera because it runs faster:

```r
library(limma)
library(WebGestaltR)  # for readGmt()

get_enrich_order2 = function( res, gene_sets ){
  if( !is.null(res$z.std) ){
    stat = res$z.std
  }else if( !is.null(res$F.std) ){
    stat = res$F.std
  }else if( !is.null(res$t) ){
    stat = res$t
  }else{
    stat = res$F
  }
  names(stat) = res$hgnc_symbol
  stat = stat[!is.na(names(stat))]
  index = ids2indices(gene_sets, names(stat))
  cameraPR( stat, index )
}

data_dir = '~/data/expression_impute/'
phenotypes = read.table(sprintf('%s/phenos.txt', data_dir))[,1]
G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]

for (phen in phenotypes) {
    if (grepl(x=phen, pattern='Caudate') || grepl(x=phen, pattern='ATR')) {
        region = 'Caudate'
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

        for (score in c('zscore', 'effect')) {
            tmp2 = G_list2[, c('hgnc_symbol', score)]
            tmp2 = tmp2[!is.na(tmp2[, score]), ]
            colnames(tmp2)[2] = 't' # quick hack so I don't have to change the function
            for (db in c('geneontology_Biological_Process_noRedundant',
                            'geneontology_Cellular_Component_noRedundant',
                            'geneontology_Molecular_Function_noRedundant',
                            'pathway_KEGG', 'disease_Disgenet',
                            'phenotype_Human_Phenotype_Ontology',
                            'network_PPI_BIOGRID',
                            'geneontology_Biological_Process',
                        'geneontology_Cellular_Component',
                        'geneontology_Molecular_Function')) {
                cat(md, score, phen, db, '\n')
                # assigning hgnc to our gene set lists
                gs = loadGeneSet(enrichDatabase=db)
                gmt = gs$geneSet
                a = idMapping(inputGene=gmt$gene, sourceIdType='entrezgene',
                            targetIdType='genesymbol')
                gmt2 = merge(gmt, a$mapped[, c('userId', 'geneSymbol')], by.x = 'gene',
                            by.y='userId', all.x=F, all.y=F)
                # and convert it to lists
                mylist = list()
                for (s in unique(gmt2$geneSet)) {
                    mylist[[s]] = unique(gmt2$geneSymbol[gmt2$geneSet==s])
                }
                res_camera = get_enrich_order2( tmp2, mylist )
                if (is.null(gs$geneSetDes)) {
                    # PPI doesn't have descriptions
                    m = cbind(rownames(res_camera), res_camera, NA)
                    colnames(m)[1] = 'Row.names'
                    colnames(m)[ncol(m)] = 'description'
                } else {
                    # attach gene set description
                    m = merge(res_camera, gs$geneSetDes, by.x=0, by.y=1)
                    m = m[order(m$PValue), ]
                    # make sure our CSV is not corrupted by extra commas
                    m$description = gsub(x=m$description, pattern=',', replacement=';')
                }
                out_fname = sprintf('%s/camera2_%s_%s_%s_%s.csv', data_dir,
                                    md, score, phen, db)
                write.csv(m, file=out_fname, quote=F, row.names=F)
            }
            # my own GMTs
            for (db in c('disorders', 'adhd_genes',
                         sprintf('%s_developmental', region))) {
                cat(md, score, phen, db, '\n')
                db_file = sprintf('~/data/post_mortem/%s.gmt', db)
                gmt = readGmt(db_file) # already in gene symbols
                # and convert it to lists
                mylist = list()
                for (s in unique(gmt$geneSet)) {
                    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
                }
                res_camera = get_enrich_order2( tmp2, mylist )
                # some massaging to look like the other results
                desc = gmt[, c(1,2)]
                desc = desc[!duplicated(desc$geneSet), ]
                m = merge(res_camera, desc, by.x=0, by.y=1)
                m = m[order(m$PValue), ]
                # make sure our CSV is not corrupted by extra commas
                m$description = gsub(x=m$description, pattern=',', replacement=';')
                out_fname = sprintf('%s/camera2_%s_%s_%s_%s.csv', data_dir,
                                    md, score, phen, db)
                write.csv(m, file=out_fname, quote=F, row.names=F)
            }
        }
    }
}
```

And then we do WebGestaltR:

```bash
# bw
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate radian
~/.local/bin/radian
```

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

        for (score in c('zscore', 'effect')) {
            tmp2 = G_list2[, c('hgnc_symbol', score)]
            # my own GMTs
            for (db in c('disorders', 'adhd_genes',
                         sprintf('%s_developmental', region))) {
                cat(md, score, phen, db, '\n')
                project_name = sprintf('%s_%s_%s_%s_v2', md, score, phen, db)
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
                                            projectName=project_name)
                out_fname = sprintf('%s/WG2_%s_%s_%s_%s_10K.csv', data_dir,
                                    md, score, phen, db)
                out_fname = gsub(x=out_fname, pattern='\\.',
                                 replacement='_')
                write.csv(enrichResult, file=out_fname, row.names=F)
            }
            for (db in c('geneontology_Biological_Process_noRedundant',
                            'geneontology_Cellular_Component_noRedundant',
                            'geneontology_Molecular_Function_noRedundant',
                            'pathway_KEGG', 'disease_Disgenet',
                            'phenotype_Human_Phenotype_Ontology',
                            'network_PPI_BIOGRID')) {
                cat(md, score, phen, db, '\n')
                project_name = sprintf('%s_%s_%s_%s_v2', md, score, phen, db)
                # make sure no dots in the name
                project_name = gsub(x=project_name, pattern='\\.',
                                    replacement='_')
                enrichResult <- WebGestaltR(enrichMethod="GSEA",
                                            organism="hsapiens",
                                            enrichDatabase=db,
                                            interestGene=tmp2,
                                            interestGeneType="genesymbol",
                                            sigMethod="top", topThr=150000,
                                            outputDirectory = data_dir,
                                            minNum=5, projectName=project_name,
                                            isOutput=T, isParallel=T,
                                            nThreads=ncpu, perNum=10000)
                out_fname = sprintf('%s/WG2_%s_%s_%s_%s_10K.csv', data_dir,
                                    md, score, phen, db)
                out_fname = gsub(x=out_fname, pattern='\\.',
                                 replacement='_')
                write.csv(enrichResult, file=out_fname, row.names=F)
            }
        }
    }
}
```

I'm actually thinking this will go much faster if I swarm it... let me get some
results with camera and once that's done I'll re-evaluate. I can also read those
two papers Philip sent in while these things run.

# 2020-11-19 18:00:15

It turns out that the overlap developmenta set has more than 500 genes, so we
need to increase that. I don't want to re-run everything, so let's just do the
developmental stuff:

```r
# bw
library(WebGestaltR)

data_dir = '~/data/expression_impute/'
phenotypes = read.table(sprintf('%s/phenos.txt', data_dir))[,1]

G_list0 = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
ncpu=31

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

        for (score in c('zscore', 'effect')) {
            tmp2 = G_list2[, c('hgnc_symbol', score)]
            db = sprintf('%s_developmental', region)
            cat(md, score, phen, db, '\n')
            project_name = sprintf('%s_%s_%s_%s_v2', md, score, phen, db)
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
                                        maxNum=1000)
            out_fname = sprintf('%s/WG2_%s_%s_%s_%s_10K.csv', data_dir,
                                md, score, phen, db)
            out_fname = gsub(x=out_fname, pattern='\\.',
                                replacement='_')
            write.csv(enrichResult, file=out_fname, row.names=F)
        }
    }
}
```

# 2020-12-03 20:00:04

Let's see if our currently chosen phenotypes show any sort of relationship to
our clinical metrics. Based on our results, we'll go with: Caudate_volume and
ACC_volume.

Let's try to do it with the WNH cohort we have been playing with. We can add
more people later if needed.

```r
phen = 'ACC_volume'
my_dir = '~/data/expression_impute/'
brain_data = read.delim(sprintf('%s/phen_res_%s.tab', my_dir, phen))
brain_data$SID = as.numeric(gsub(x=brain_data$IID, replacement = '',
                                 pattern = 'SID.'))

clin = read.csv('~/data/expression_impute//augmented_anon_clinical_10242020.csv')
clin_slim = clin[, c('SID', 'everADHD_dsm', 'maxOverTimeSX_inatt',
                     'maxOverTimeSX_hi')]
clin_slim = clin_slim[!duplicated(clin_slim$SID),]
m = merge(brain_data, clin_slim, by='SID')

num_vars = c('maxOverTimeSX_inatt', 'maxOverTimeSX_hi')
pvals = c()
for (x in num_vars) {
    res = cor.test(m[, x], m$phen)
    pvals = c(pvals, res$p.value)
}
categ_vars = c('everADHD_dsm')
for (x in categ_vars) {
    idx = m[,x] == levels(factor(m[,x]))[1]
    res = t.test(m[idx, 'phen'], m[!idx, phen])
    pvals = c(pvals, res$p.value)
}
names(pvals) = c(num_vars, categ_vars)
print(pvals)
```

```
(Caudate)
maxOverTimeSX_inatt    maxOverTimeSX_hi        everADHD_dsm 
          0.3590426           0.6421900           0.7521822 
(ACC)
maxOverTimeSX_inatt    maxOverTimeSX_hi        everADHD_dsm 
          0.3970843           0.7550794           0.4280699 
```

Nothing here.

# 2020-12-04 11:44:34

How about adding everyone we have? And we can do an LME:

```r
fname = '/Volumes/NCR/brain_data_compiled/all_freesurfer_12042020.csv'
df = read.csv(fname)
# restrict it to within our age range
df = df[df$age_scan >= 6.69 & df$age_scan <= 38.83,]
# keep_idx = (df$mprage_score <= 2.5 & df$ext_avg_score <=2.5 &
#             df$int_avg_score <= 2.5)
keep_idx = (df$mprage_score <= 2 & df$ext_avg_score <=2 &
            df$int_avg_score <= 2)
df = df[keep_idx,]
# keep only the best scan for each subject
brain_data = c()
for (sid in unique(df$SID)) {
    sdata = df[which(df$SID == sid), ]
    if (nrow(sdata) == 1) {
        brain_data = rbind(brain_data, sdata)
    } else {
        scores = rowMeans(sdata[, c('mprage_score', 'ext_avg_score',
                                    'int_avg_score')])
        brain_data = rbind(brain_data, sdata[which.min(scores), ])
    }
}
# bring in the clinical variables
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/expression_impute//augmented_anon_clinical_10242020.csv')
clin$SX_total = clin$SX_inatt + clin$SX_hi
clin$maxOverTimeSX_total = clin$maxOverTimeSX_inatt + clin$maxOverTimeSX_hi
clin_slim = clin[clin$age_clin!='child',]
clin_slim$age_clin = as.numeric(clin_slim$age_clin)
data2 = mergeOnClosestAge(brain_data, clin_slim, brain_data$SID, x.id='SID',
                          y.id='SID', x.age='age_scan', y.age='age_clin')
data2$ACC_vol = (data2$lh_caudalanteriorcingulate_volume +
                 data2$rh_caudalanteriorcingulate_volume)
data2$Caudate_vol = (data2$Left.Caudate + data2$Right.Caudate)
```

Now it's just a matter of setting up the LMEs. First, let's try it with stepAIC:

```r
library(nlme)
library(MASS)
data3 = data2[which(data2$DX_dsm != 'other'),]
brain_vars = c('ACC_vol', 'Caudate_vol')
clin_vars = colnames(clin)[c(3:4, 6:8, 12:13, 15:16)]
all_res = c()
for (cl in clin_vars) {
    for (v in brain_vars) {
        fm_str = sprintf('%s ~ %s + sex + scanner_update + age_scan + mprage_score + ext_avg_score + int_avg_score', v, cl)
        fit <- lme(as.formula(fm_str), random=~1|famID, data = data3,
                   na.action=na.omit, method='ML')
        myres = summary(fit)$tTable[2, ]
        myres$brain = v
        myres$clin = cl
        all_res = rbind(all_res, myres)
        # step <- stepAIC(fit, direction = "both", trace = F, na.action=na.exclude)
    }
}
```

Nothing here either...

# 2020-12-16 20:54:59

I'm finding it hard to believe that our brain phenotypes don't have any association
with the ADHD clinical phenotypes. It might be easier to go backwards after we
find a good phenotype, following the rationale I outlined in note 163. Of course
it could also be done in ABCD, and since I'm waiting on that PrediXcan
imputation, and I'm also waiting on the download of the S-PrediXcan models, I
might as well get the code ready to run our data and then it's just a matter of
running ABCD as well.

If that still doesn't work, even for ABCD, then I could come up with a
meta-phenotype. Basically, use PLS to predict ADHD using all our brain
phenotypes, and use that linear combination of brain phenotypes as our
meta-phenotype to be regressed against the predicted gene expression. That might
actually be a better approach in the sense that we don't know for sure which
brain phenotype is supposed to be most linked to gene expression in the brain.
It's just harder to understand.

# 2020-12-17 17:12:28

```r
fname = '/Volumes/NCR/brain_data_compiled/all_freesurfer_12042020.csv'
df = read.csv(fname)
# restrict it to within our age range
df = df[df$age_scan >= 6.69 & df$age_scan <= 38.83,]
keep_idx = (df$mprage_score <= 2.5 & df$ext_avg_score <=2.5 &
            df$int_avg_score <= 2.5)
# keep_idx = (df$mprage_score <= 2 & df$ext_avg_score <=2 &
#             df$int_avg_score <= 2)
df = df[keep_idx,]
# keep only the best scan for each subject
brain_data = c()
for (sid in unique(df$SID)) {
    sdata = df[which(df$SID == sid), ]
    if (nrow(sdata) == 1) {
        brain_data = rbind(brain_data, sdata)
    } else {
        scores = rowMeans(sdata[, c('mprage_score', 'ext_avg_score',
                                    'int_avg_score')])
        brain_data = rbind(brain_data, sdata[which.min(scores), ])
    }
}
# whole rhACC
brain_data$rh_ACC_volume = brain_data$rh_caudalanteriorcingulate_volume + brain_data$rh_rostralanteriorcingulate_volume
brain_data$rh_ACC_thickness = brain_data$rh_caudalanteriorcingulate_thickness + brain_data$rh_rostralanteriorcingulate_thickness
brain_data$rh_ACC_area = brain_data$rh_caudalanteriorcingulate_area + brain_data$rh_rostralanteriorcingulate_area
# whole caACC
brain_data$caACC_volume = brain_data$lh_caudalanteriorcingulate_volume + brain_data$rh_caudalanteriorcingulate_volume
brain_data$caACC_area = brain_data$lh_caudalanteriorcingulate_area + brain_data$rh_caudalanteriorcingulate_area
brain_data$caACC_thickness = brain_data$lh_caudalanteriorcingulate_thickness + brain_data$rh_caudalanteriorcingulate_thickness
# whole roACC
brain_data$roACC_volume = brain_data$lh_rostralanteriorcingulate_volume + brain_data$rh_rostralanteriorcingulate_volume
brain_data$roACC_area = brain_data$lh_rostralanteriorcingulate_area + brain_data$rh_rostralanteriorcingulate_area
brain_data$roACC_thickness = brain_data$lh_rostralanteriorcingulate_thickness + brain_data$rh_rostralanteriorcingulate_thickness
# whole ACC
brain_data$ACC_volume = brain_data$roACC_volume + brain_data$caACC_volume
brain_data$ACC_area = brain_data$roACC_area + brain_data$caACC_area
brain_data$ACC_thickness = brain_data$roACC_thickness + brain_data$caACC_thickness

brain_vars = c()
for (p in c('area', 'thickness', 'volume')) {
    for (bv in c('rh_ACC', 'caACC', 'roACC', 'ACC',
                 'rh_rostralanteriorcingulate', 'rh_caudalanteriorcingulate')) {
        brain_vars = c(brain_vars, sprintf('%s_%s', bv, p))
    }
}

# bring in the clinical variables
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('~/data/expression_impute//augmented_anon_clinical_10242020.csv')
clin$SX_total = clin$SX_inatt + clin$SX_hi
clin$maxOverTimeSX_total = clin$maxOverTimeSX_inatt + clin$maxOverTimeSX_hi
clin_slim = clin[clin$age_clin!='child',]
clin_slim$age_clin = as.numeric(clin_slim$age_clin)
data2 = mergeOnClosestAge(brain_data, clin_slim, brain_data$SID, x.id='SID',
                          y.id='SID', x.age='age_scan', y.age='age_clin')

# remove outliers
for (v in brain_vars) {
    m = mean(data2[, v], na.rm=T)
    s = sd(data2[, v], na.rm=T)
    data2[which(data2[, v] > m + 3*s), v] = NA
    data2[which(data2[, v] < m - 3*s), v] = NA
}
```

And we run the LMEs to use as much data as we can:

```r
library(nlme)
library(MASS)
data3 = data2[which(data2$DX_dsm != 'other'),]
data3$DX_dsm = factor(data3$DX_dsm)
data3$everADHD_dsm = factor(data3$everADHD_dsm)
data3$sex = factor(data3$sex)
data3$scanner_update = factor(data3$scanner_update)
clin_vars = colnames(clin)[c(3:4, 6:8, 12:13, 15:16)]
all_res = c()
for (cl in clin_vars) {
    cat(cl, '\n')
    for (v in brain_vars) {
        fm_str = sprintf('%s ~ %s + sex + scanner_update + age_scan + mprage_score + ext_avg_score + int_avg_score', v, cl)
        fit <- lme(as.formula(fm_str), random=~1|famID, data = data3,
                   na.action=na.omit, method='ML')
        step <- stepAIC(fit, direction = "both", trace = F,
                        scope = list(lower = as.formula(sprintf('~ %s', cl))))
        myres = summary(step)$tTable[2, ]
        myres$brain = v
        myres$clin = cl
        all_res = rbind(all_res, myres)
    }
}
```

That gives me 700 subjects (scans), but there was nothing there... I could try
the normalized version of the brain metrics, but let's first try other
populations:

```r
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
pcs$SID = as.numeric(gsub(x=pcs$IID, pattern='SID.', replacement=''))
imp_data = merge(brain_data, pcs, by='SID', all.x=F, all.y=F)
imwnh = imp_data$PC01<0 & imp_data$PC02>-.02

data2 = mergeOnClosestAge(imp_data[imwnh, ], clin_slim, imp_data$SID, x.id='SID',
                          y.id='SID', x.age='age_scan', y.age='age_clin')

# remove outliers
for (v in brain_vars) {
    m = mean(data2[, v], na.rm=T)
    s = sd(data2[, v], na.rm=T)
    data2[which(data2[, v] > m + 3*s), v] = NA
    data2[which(data2[, v] < m - 3*s), v] = NA
}

data3 = data2[which(data2$DX_dsm != 'other'),]
data3$DX_dsm = factor(data3$DX_dsm)
data3$everADHD_dsm = factor(data3$everADHD_dsm)
data3$sex = factor(data3$sex)
data3$scanner_update = factor(data3$scanner_update)
clin_vars = colnames(clin)[c(3:4, 6:8, 12:13, 15:16)]
all_res = c()
for (cl in clin_vars) {
    cat(cl, '\n')
    for (v in brain_vars) {
        fm_str = sprintf('%s ~ %s + sex + scanner_update + age_scan + mprage_score + ext_avg_score + int_avg_score', v, cl)
        fit <- try(lme(as.formula(fm_str), random=~1|famID, data = data3,
                   na.action=na.omit, method='ML'))
        if (length(fit) < 2) {
            myres = NA
        } else {
            step <- try(stepAIC(fit, direction = "both", trace = F,
                        scope = list(lower = as.formula(sprintf('~ %s', cl)))))
            if (length(step) < 2) {
                myres = summary(fit)$tTable[2, ]
            } else {
                myres = summary(step)$tTable[2, ]
            }
            myres$brain = v
            myres$clin = cl
        }
        all_res = rbind(all_res, myres)
    }
}
```

Now we're down to 403 scans, and still no significant results. Let's do just one
more cut, and include a single person per family:

```r
data3 = c()
for (f in unique(data2$famID)) {
    fdata = data2[data2$famID == f, ]
    if (nrow(fdata) == 1) {
        data3 = rbind(data3, fdata)
    } else {
        qc_score = rowMeans(fdata[, c('mprage_score', 'ext_avg_score',
                                      'int_avg_score')])
        data3 = rbind(data3, fdata[which.min(qc_score), ])
    }
}
for (v in brain_vars) {
    m = mean(data3[, v], na.rm=T)
    s = sd(data3[, v], na.rm=T)
    data3[which(data3[, v] > m + 3*s), v] = NA
    data3[which(data3[, v] < m - 3*s), v] = NA
}

data4 = data3[which(data3$DX_dsm != 'other'),]
data4$DX_dsm = factor(data4$DX_dsm)
data4$everADHD_dsm = factor(data4$everADHD_dsm)
data4$sex = factor(data4$sex)
data4$scanner_update = factor(data4$scanner_update)
clin_vars = colnames(clin)[c(3:4, 6:8, 12:13, 15:16)]
all_res = c()
for (cl in clin_vars) {
    cat(cl, '\n')
    for (v in brain_vars) {
        fm_str = sprintf('%s ~ %s + sex + scanner_update + age_scan + mprage_score + ext_avg_score + int_avg_score', v, cl)
        fit <- lm(as.formula(fm_str), data = data4, na.action=na.omit)
        step <- stepAIC(fit, direction = "both", trace = F,
                        scope = list(lower = as.formula(sprintf('~ %s', cl))))
        myres = summary(step)$coefficients[2, ]
        myres$brain = v
        myres$clin = cl
        myres$fm1 = fm_str
        myres$fm2 = as.character(step$terms)[3]
        all_res = rbind(all_res, myres)
    }
}
```

```
r$> all_res[all_res[, 'Pr(>|t|)'] < .05, c('t value', 'Pr(>|t|)', 'brain', 'clin', 'fm2')]                          
      t value   Pr(>|t|)    brain                                  clin                 
myres 3.098029  0.002192765 "rh_caudalanteriorcingulate_thickness" "SX_hi"              
myres 2.30815   0.02185271  "caACC_area"                           "DX_dsm"             
myres 2.007541  0.04583219  "caACC_volume"                         "DX_dsm"             
myres 2.350956  0.01956695  "rh_caudalanteriorcingulate_thickness" "everADHD_dsm"       
myres -2.057895 0.0408471   "caACC_area"                           "outcome_dsm"        
myres 1.976865  0.04941067  "caACC_thickness"                      "outcome_dsm"        
myres 2.7882    0.005801412 "rh_caudalanteriorcingulate_thickness" "outcome_dsm"        
myres -2.02599  0.04388622  "ACC_area"                             "maxOverTimeSX_inatt"
myres -2.0484   0.04162741  "rh_rostralanteriorcingulate_area"     "maxOverTimeSX_inatt"
myres -2.092556 0.03745286  "caACC_area"                           "maxOverTimeSX_hi"   
myres -2.112829 0.03566145  "ACC_area"                             "maxOverTimeSX_hi"   
myres 1.993297  0.04740029  "caACC_thickness"                      "maxOverTimeSX_hi"   
myres 2.791436  0.005686526 "rh_caudalanteriorcingulate_thickness" "maxOverTimeSX_hi"   
myres -2.052371 0.04124636  "caACC_area"                           "SX_total"           
myres 2.623474  0.009291295 "rh_caudalanteriorcingulate_thickness" "SX_total"           
myres -2.090015 0.03768274  "caACC_area"                           "maxOverTimeSX_total"
myres -2.205685 0.02836651  "ACC_area"                             "maxOverTimeSX_total"
myres 2.173595  0.03075196  "rh_caudalanteriorcingulate_thickness" "maxOverTimeSX_total"
      fm2                                                                    
myres "SX_hi + sex + scanner_update + age_scan + ext_avg_score"              
myres "DX_dsm + sex + int_avg_score"                                         
myres "DX_dsm + sex + age_scan + int_avg_score"                              
myres "everADHD_dsm + sex + scanner_update + age_scan + ext_avg_score"       
myres "outcome_dsm + sex + int_avg_score"                                    
myres "outcome_dsm + sex + scanner_update + age_scan + ext_avg_score"        
myres "outcome_dsm + scanner_update + age_scan + ext_avg_score"              
myres "maxOverTimeSX_inatt + sex + int_avg_score"                            
myres "maxOverTimeSX_inatt + sex + int_avg_score"                            
myres "maxOverTimeSX_hi + sex + int_avg_score"                               
myres "maxOverTimeSX_hi + sex + int_avg_score"                               
myres "maxOverTimeSX_hi + sex + scanner_update + age_scan + ext_avg_score"   
myres "maxOverTimeSX_hi + sex + scanner_update + age_scan + ext_avg_score"   
myres "SX_total + sex + int_avg_score"                                       
myres "SX_total + sex + scanner_update + age_scan + ext_avg_score"           
myres "maxOverTimeSX_total + sex + int_avg_score"                            
myres "maxOverTimeSX_total + sex + int_avg_score"                            
myres "maxOverTimeSX_total + sex + scanner_update + age_scan + ext_avg_score"
```

We get a few nice hits here. That's good as it's our ideal population too, even
though it only has 241 scans.

Now we need to check the association between the brain phenotype and the imputed
expression. We will need to do this without covariates, fixed covariates, and
then just the best fit for DX.


# TODO
 * is the direction of the association correct?
