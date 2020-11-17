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
                write.csv(enrichResult, file=out_fname, quote=F,
                            row.names=F)
            }
            for (db in c('geneontology_Biological_Process_noRedundant',
                            'geneontology_Cellular_Component_noRedundant',
                            'geneontology_Molecular_Function_noRedundant',
                            'pathway_KEGG', 'disease_Disgenet',
                            'phenotype_Human_Phenotype_Ontology',
                            'network_PPI_BIOGRID')) {
                cat(md, score, phen, db, '\n')
                project_name = sprintf('%s_%s_%s_%s_v2', md, score, phen, db)
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
                write.csv(enrichResult, file=out_fname, quote=F,
                            row.names=F)
            }
        }
    }
}
```

I'm actually thinking this will go much faster if I swarm it... let me get some
results with camera and once that's done I'll re-evaluate. I can also read those
two papers Philip sent in while these things run.





# TODO

