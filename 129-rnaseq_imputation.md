# 2020-09-04 07:01:22

Kwangmi handed me the newly cleaned genotyping for all our cohort. It's in
/NCR_SBRB/Ahn/genotying/allfamily_clean2_geno0.05_mind0.02_geno0.02_hwe_finalQC
and it also includes a README file. Let's use that for 1KG imputation.

```bash
cd /data/NCR_SBRB/NCR_genomics/genotyping/
mkdir v3
cd v3
cp /data//NCR_SBRB/Ahn/genotying/allfamily_clean2_geno0.05_mind0.02_geno0.02_hwe_finalQC* .
cp /data//NCR_SBRB/Ahn/genotying/readme.chip.txt .
```

Before I do that, let's change Kwangmi's file because it has MRNs. I'll just use
Labmatrix IDs instead so we can send it to Michigan.

```bash
cut -d" " -f 1,2 allfamily_clean2_geno0.05_mind0.02_geno0.02_hwe_finalQC.fam > mrns.txt
# clean up locally and attached Labmatrix SID
module load plink
plink --bfile allfamily_clean2_geno0.05_mind0.02_geno0.02_hwe_finalQC \
    --update-ids mrns2sids.txt --make-bed --out NCR_v3
```

Going back to setting up the files for imputation

```bash
module load vcftools
module load plink
for i in {1..22}; do
    plink --bfile NCR_v3 --chr $i --recode vcf --out NCR_v3_chr"$i";
    vcf-sort NCR_v3_chr"$i".vcf | bgzip -c > NCR_v3_chr"$i".vcf.gz
done
```

I'll try to just upload these and see if they work, otherwise I might need to flip some alleles.

No go... more than 100 obvious strand flips, so let's use shape it. The first
time I ran it I had issues with duplicates, so let's clean that up first.

```bash
plink --bfile NCR_v3 --list-duplicate-vars ids-only suppress-first
plink --bfile NCR_v3 --exclude plink.dupvar --make-bed --out NCR_v3_noDups
for i in {1..22}; do
    plink --bfile NCR_v3_noDups --chr $i --recode vcf --out NCR_v3_noDups_chr"$i";
    vcf-sort NCR_v3_noDups_chr"$i".vcf | bgzip -c > NCR_v3_noDups_chr"$i".vcf.gz
done
```

```bash
module load shapeit/2.r904
refdir=/fdb/impute2/1000Genomes_Phase3_integrated_haplotypes_Oct2014/1000GP_Phase3/
for c in {1..22}; do
    shapeit -check -T 16 -V NCR_v3_noDups_chr${c}.vcf.gz \
        --input-ref $refdir/1000GP_Phase3_chr${c}.hap.gz \
        $refdir/1000GP_Phase3_chr${c}.legend.gz $refdir/1000GP_Phase3.sample \
        --output-log chr${c}.alignments;
done
# format the files:
for c in {1..22}; do
    grep Strand chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> flip_snps.txt;
    grep Missing chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> missing_snps.txt;
done
plink --bfile NCR_v3_noDups --flip flip_snps.txt \
    --exclude missing_snps.txt --make-bed --out NCR_v3_noDups_flipped
#reconstruct the VCFs as above to send it to the imputation server.
for i in {1..22}; do
    plink --bfile NCR_v3_noDups_flipped --chr ${i} --recode-vcf \
        --out NCR_v3_noDups_flipped_chr${i};
    vcf-sort NCR_v3_noDups_flipped_chr"$i".vcf | bgzip \
        -c > NCR_v3_noDups_flipped_chr"$i".vcf.gz
done
```

And then we try the imputation server again.

The plan now is to use PrediXcanAssociation.py in the ABCD and NCR imputed data.
The question is on what the phenotypes and covariates should be. Let's go for
everADHD, ADHD, age, sex, population code, and the 10 population PCs.

Then we can somewhow analyze whether there is an overlap between variants that
are significant between post-mortem and either ABCD or NCR. It'd have to be a
big enough overlap set, verified by some sort of permutation. And then we should
look for replicates in the other imputed set. The cherry on the cake would be
some sort of neuroimaging finding. We can explore different age groups and
population subgroups to make this work too. Or even use one set of significant
genes from one analysis to filter the other and run a better FDR?

Worst case scenario we do gene sets in all cohorts, using the stats from
PrediXcanAssociation.py for the imputation data sets.

Let's also capture the commands I'm using to run the imputation:

```bash
# bw
module load python
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate imlabtools

mydir=~/data/expression_impute;
python3 $mydir/MetaXcan-master/software/Predict.py \
    --model_db_path $mydir/eqtl/mashr/mashr_Brain_Caudate_basal_ganglia.db \
    --vcf_genotypes /data/NCR_SBRB/ABCD/v201/1KG/chr*.dose.vcf.gz \
    --vcf_mode imputed \
    --prediction_output $mydir/results/ABCD_v201_Caudate_predict_1KG_mashr.txt \
    --prediction_summary_output $mydir/results/ABCD_v201_Caudate_summary_1KG_mashr.txt \
    --verbosity 9 --throw --model_db_snp_key varID \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --liftover $mydir/hg19ToHg38.over.chain.gz

mydir=~/data/expression_impute;
python3 $mydir/MetaXcan-master/software/Predict.py \
    --model_db_path $mydir/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.db \
    --vcf_genotypes /data/NCR_SBRB/ABCD/v201/1KG/chr*.dose.vcf.gz \
    --vcf_mode imputed \
    --prediction_output $mydir/results/ABCD_v201_ACC_predict_1KG_mashr.txt \
    --prediction_summary_output $mydir/results/ABCD_v201_ACC_summary_1KG_mashr.txt \
    --verbosity 9 --throw --model_db_snp_key varID \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --liftover $mydir/hg19ToHg38.over.chain.gz
```

And for NCR data is similar, for example:

```bash
mydir=~/data/expression_impute;
python3 $mydir/MetaXcan-master/software/Predict.py \
    --model_db_path $mydir/eqtl/mashr/mashr_Brain_Caudate_basal_ganglia.db \
    --vcf_genotypes /data/NCR_SBRB/NCR_genomics/genotyping/v3/chr*.dose.vcf.gz \
    --vcf_mode imputed \
    --prediction_output $mydir/results/NCR_v3_Caudate_predict_1KG_mashr.txt \
    --prediction_summary_output $mydir/results/NCR_v3_Caudate_summary_1KG_mashr.txt \
    --verbosity 9 --throw --model_db_snp_key varID \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --liftover $mydir/hg19ToHg38.over.chain.gz

mydir=~/data/expression_impute;
python3 $mydir/MetaXcan-master/software/Predict.py \
    --model_db_path $mydir/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.db \
    --vcf_genotypes /data/NCR_SBRB/NCR_genomics/genotyping/v3/chr*.dose.vcf.gz \
    --vcf_mode imputed \
    --prediction_output $mydir/results/NCR_v3_ACC_predict_1KG_mashr.txt \
    --prediction_summary_output $mydir/results/NCR_v3_ACC_summary_1KG_mashr.txt \
    --verbosity 9 --throw --model_db_snp_key varID \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --liftover $mydir/hg19ToHg38.over.chain.gz
```

Now that imputation results are ready:

```bash
cd /data/NCR_SBRB/NCR_genetics/v3
for f in `/bin/ls *zip`; do unzip -P W5hTW7yDJIQxry $f; done
```

From Philip:

```
[5:08 PM] Shaw, Philip (NIH/NHGRI) [E]
    yup--- the neatest finding would be that genes that we find to have grex for case vs control in the brain overlap with those in the living subject---- and the 'grex' may be gene-set overlap ---enriching the same developmental periods, but not necessarilty wiht exactly the diff exp genes.  
​[5:09 PM] Shaw, Philip (NIH/NHGRI) [E]
    If that holds---then you could see how the 'grex' genes associate with brain dimension in the living subject----then concluding that these grex genes have a trophic impact on the brain.  
```

# 2020-09-07 12:25:06

Just to make things run faster, I did:

```r
a = read.table('~/data/expression_impute/results/NCR_v3_ACC_predict_1KG_mashr.txt', header=1)
saveRDS(a, '~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
```

I'll also use the same data we used for imputation in order to compute the
population PCs, but here I'll need to remove that last sample so it's not
contributing to the PCs:

```bash
#bw
cd /data/NCR_SBRB/NCR_genomics/genotyping/v3
echo "3333333 3239" >> rm_ids.txt
module load plink
plink --bfile NCR_v3_noDups_flipped --remove rm_ids.txt --make-bed \
    --out NCR_v3_noDups_flipped_noExtra
/data/NCR_SBRB/software/KING/king --mds --cpus 30 -b NCR_v3_noDups_flipped_noExtra.bed
```

```r
pcs = read.table('/data/NCR_SBRB/NCR_genomics/genotyping/v3/kingpc.ped')
pcsf = pcs[, c(2, 7:26)]
new_names = c('IID', sapply(1:20, function(x) sprintf('PC%.2d', x)))
colnames(pcsf) = new_names
write.csv(pcsf, file='/data/NCR_SBRB/NCR_genomics/genotyping/v3/pop_pcs.csv', row.names=F)
```

Now, let's try to do some basic analysis, starting by looking at WNH only,
because if I remember correctly that's our strongest result in post-mortem. It's
always possible that the expression imputation is heavily skewed to WNH
populations, so let's start there.

![](images/2020-09-07-13-39-00.png)

Here it'll depend on how many people we want to include. Let's do PC01 < 0 and
PC02 > -.025 for now.

```r
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('/data/NCR_SBRB/NCR_genomics/genotyping/v3/pop_pcs.csv')
keep_me = pcs$PC01 < 0 & pcs$PC02 > -.025
pcs = pcs[keep_me, ]
data = merge(a, pcs, by='IID', all.x=F, all.y=F)
gf = read.csv('~/data/expression_impute/gf_1119_09072020.csv')
data = merge(data, gf, by.x='IID', by.y='Subject.Code...Subjects', all.x=F, all.y=F)
data$Case = factor(data$Case)
```

And the analysis is similar to what we normally do:

```r
library(lme4)
library(car)
fit = lmer(ENSG00000261456.5 ~ Case + Sex...Subjects + PC01 + (1|FAMID), data=data)
Anova(fit)
```

or, if we run for all genes:

```r
library(caret)

grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]

# remove zero or near zero variance genes
set.seed(42)
pp_order = c('zv', 'nzv')
pp = preProcess(data[, grex_vars], method = pp_order)
X = predict(pp, data[, grex_vars])
good_grex = colnames(X)

# set up some heuristic that at least half the subjects should have non-zero imputation
nzeros = colSums(X==0)
good_grex = good_grex[nzeros < (nrow(X)/2)]
pvals = sapply(good_grex, function(x) {
    fm_str = sprintf('%s ~ Case + Sex...Subjects + (1|FAMID)', x)
    fit = lmer(as.formula(fm_str), data=data)
    return(Anova(fit)['Case', 3])})
```

And I could easily parallelize that to go much faster. We just need to figure
out if we are using the right covariates.

I'm seeing scenarios where we have many subjects with 0 for some genes. I'll
need to use whatever heuristic I'm applying for the post-mortem here as well
then.

# 2020-09-08 06:56:34

Let's now focus on how to remove bad genes (not represented in enough subjects)
and how to massage the results to run the gene set analysis.

I'm using filterByExpr for rnaseq, so I'll need to generate some plots to
qualitatively crop the imputed data.

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier that's NOT caudate labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)
data = data[data$Sex=='M',]

imWNH = which(data$C1 > 0 & data$C2 < -.065)
data = data[imWNH, ]

# have 19 super clean WNH subjects for ACC and 20 for Caudate
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
# data matrix goes on a diet...
data = data[, !grepl(colnames(data), pattern='^ENS')]
# remove that weird .num after ENSG
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
rownames(count_matrix) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
count_matrix = count_matrix[!dups, ]

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
# remove any genes without a HUGOID
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
# remove genes that appear more than once
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
# keep only gene counts for genes that we have information
imnamed = rownames(count_matrix) %in% G_list$ensembl_gene_id
count_matrix = count_matrix[imnamed, ]

# some data variables modifications
data$POP_CODE = as.character(data$POP_CODE)
data[data$POP_CODE=='WNH', 'POP_CODE'] = 'W'
data[data$POP_CODE=='WH', 'POP_CODE'] = 'W'
data$POP_CODE = factor(data$POP_CODE)
data$Individual = factor(data$hbcc_brain_id)
data[data$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
data[data$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
data$MoD = factor(data$Manner.of.Death)
data$batch = factor(as.numeric(data$run_date))

library(caret)
pp_order = c('zv', 'nzv')
pp = preProcess(t(count_matrix), method = pp_order)
X = predict(pp, t(count_matrix))
geneCounts = t(X)

# match gene counts to gene info
G_list2 = merge(rownames(geneCounts), G_list, by=1)
colnames(G_list2)[1] = 'ensembl_gene_id'

# keep only autosomal genes
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
geneCounts = geneCounts[imautosome, ]
G_list2 = G_list2[imautosome, ]

library(edgeR)
isexpr <- filterByExpr(geneCounts, group=data$Diagnosis)
genes = DGEList( geneCounts[isexpr,], genes=G_list2[isexpr,] ) 
```

Now we make plots with one line per subject:

```r
plot(density(cpm(geneCounts[,1], log=T)), las=2, main="", xlab="lCPM")
for (s in 2:nrow(data)) {
  lines(density(cpm(geneCounts[, s], log=T)), las=2)
}
# and now we plot in red after cleaning
for (s in 1:nrow(data)) {
  lines(density(cpm(genes[, s], log=T)), las=2, col='red')
}
```

![](images/2020-09-08-09-27-50.png)

Now, let's do the same thing with the imputed data:

```r
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
keep_me = pcs$PC01 < 0 & pcs$PC02 > -.025
pcs = pcs[keep_me, ]
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
gf = read.csv('~/data/expression_impute/gf_1119_09072020.csv')
imp_data = merge(imp_data, gf, by.x='IID', by.y='Subject.Code...Subjects', all.x=F, all.y=F)
imp_data$Case = factor(imp_data$Case)

library(caret)
grex_vars = colnames(imp_data)[grepl(colnames(imp_data), pattern='^ENS')]
# remove zero or near zero variance genes
set.seed(42)
pp_order = c('zv', 'nzv')
pp = preProcess(imp_data[, grex_vars], method = pp_order)
X = predict(pp, imp_data[, grex_vars])
grex_vars = colnames(X)

# data matrix goes on a diet...
imp_data = imp_data[, !grepl(colnames(imp_data), pattern='^ENS')]
# remove that weird .num after ENSG
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
colnames(X) = id_num
grex_vars = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
X = t(X[, !dups])

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
# remove any genes without a HUGOID
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
# remove genes that appear more than once
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
# keep only imputed gene expression for genes that we have information
imnamed = rownames(X) %in% G_list$ensembl_gene_id
X = X[imnamed, ]
grex_vars = grex_vars[imnamed]

# match gene counts to gene info
G_list2 = merge(G_list, X, by.x=1, by.y=0)

# keep only autosomal genes
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
G_list2 = G_list2[imautosome, ]
X = G_list2[, 4:ncol(G_list2)]
rownames(X) = G_list2$ensembl_gene_id
grex_vars = G_list2$ensembl_gene_id
```

And let's make a similar plot:

```r
plot(density(X[, 1]), las=2, main="", xlab="impVal")
for (s in 2:col(X)) {
  lines(density(X[, s]), las=2)
}
```

![](images/2020-09-08-09-45-51.png)

Hum... these values are quite odd. What are we predicting here? According to
Sam,:

```
Negative values represent underexpression of the gene, positive values represent overexpression. I am trying to remember what the parameters were to determine the reference expression, whether that is entirely taken from the reference cohort with genotype and expression data, or if it is in relation to the imputed cohort

Genes were retained if they had an R2 <0.01 for cross-validation between imputed expression and actual expression within the GTEx reference cohort. 12,419 genes were expressed in at least one cohort and 7,392 genes were expressed in all three. Following further quality control removing genes with at least 30% zero expression values in at least one cohort, 6,575 genes were included in the final analyses.
```

But the MASHR models don't give R2, so we're stuck with looking at % of zeros,
which looks funky too:

```r
nzeros = rowSums(X==0)
plot(sort(nzeros, decreasing=T))
```

![](images/2020-09-08-11-33-09.png)

So, there isn't really a good cut-off in terms of subjects (Y axis). For now,
let's stick with 50%, and we can go up or lower later.

```r
nzeros = rowSums(X==0)
good_grex = grex_vars[nzeros < (ncol(X)/2)]
library(lme4)
library(car)
mydata = cbind(imp_data, t(X))
pvals = sapply(good_grex, function(x) {
    fm_str = sprintf('%s ~ Case + Sex...Subjects + (1|FAMID)', x)
    fit = lmer(as.formula(fm_str), data=mydata)
    return(Anova(fit)['Case', 3])})
```

Going back to RNAseq:

```r
genes = calcNormFactors( genes)
form = ~ Diagnosis + scale(RINe) + scale(PMI) + scale(Age) + MoD
design = model.matrix( form, data)
vobj_tmp = voom( genes, design, plot=FALSE)
dupcor <- duplicateCorrelation(vobj_tmp, design, block=data$batch)
vobj = voom( genes, design, plot=FALSE, block=data$batch,
             correlation=dupcor$consensus)
fit <- lmFit(vobj, design, block=data$batch, correlation=dupcor$consensus)
fitDC <- eBayes( fit )
resDC = topTable(fitDC, coef='DiagnosisControl', number=Inf)

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
  # print(head(stat))
  index = ids2indices(gene_sets, names(stat))
  cameraPR( stat, index )
}
load('~/data/rnaseq_derek/adhd_genesets_philip.RDATA')
load('~/data/rnaseq_derek/c5_gene_sets.RData')
load('~/data/rnaseq_derek/brain_disorders_gene_sets.RData')
load('~/data/rnaseq_derek/data_for_alex.RData')
co = .9 
idx = anno$age_category==1 & anno$cutoff==co
genes_overlap = unique(anno[idx, 'anno_gene'])
for (s in 2:5) {
  idx = anno$age_category==s & anno$cutoff==co
  g2 = unique(anno[idx, 'anno_gene'])
  genes_overlap = intersect(genes_overlap, g2)
}
genes_unique = list()
for (s in 1:5) {
  others = setdiff(1:5, s)
  idx = anno$age_category==s & anno$cutoff==co
  g = unique(anno[idx, 'anno_gene'])
  for (s2 in others) {
    idx = anno$age_category==s2 & anno$cutoff==co
    g2 = unique(anno[idx, 'anno_gene'])
    rm_me = g %in% g2
    g = g[!rm_me]
  }
  genes_unique[[sprintf('dev%s_c%.1f', s, co)]] = unique(g)
}
genes_unique[['overlap']] = unique(genes_overlap)

adhd_dream_cameraDC = get_enrich_order2( resDC, t2 ) 
c5_dream_cameraDC = get_enrich_order2( resDC, c5_all)
dis_dream_cameraDC = get_enrich_order2( resDC, disorders)
dev_dream_cameraDC = get_enrich_order2( resDC, genes_unique )
```

And let's check that the enrichment functions run here:

```r
nzeros = rowSums(X==0)
keep_me = nzeros < (ncol(X)/2)
good_grex = grex_vars[keep_me]
library(lme4)
library(car)
mydata = cbind(imp_data, t(X))
chisq = sapply(good_grex, function(x) {
    fm_str = sprintf('%s ~ Case + Sex...Subjects + (1|FAMID)', x)
    fit = lmer(as.formula(fm_str), data=mydata)
    return(Anova(fit)['Case', 1])})
imp_res = cbind(G_list2[keep_me, 1:2], chisq)
colnames(imp_res)[3] = 'F'  # quick hack to work with function

adhd_imp_camera = get_enrich_order2( imp_res, t2 )
c5_imp_camera = get_enrich_order2( imp_res, c5_all)
dis_imp_camera = get_enrich_order2( imp_res, disorders)
dev_imp_camera = get_enrich_order2( imp_res, genes_unique )
```

Or, to make it faster:

```r
library(doParallel)

ncores = 2 #as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
cl = makeCluster(ncores)
nzeros = rowSums(X==0)
keep_me = nzeros < (ncol(X)/2)
good_grex = grex_vars[keep_me]
mydata = cbind(imp_data, t(X))
clusterExport(cl, c("mydata"))
chisq = parSapply(cl, good_grex,
                  function(x) {
                    fm_str = sprintf('%s ~ Case + Sex...Subjects + (1|FAMID)', x)
                    fit = lme4::lmer(as.formula(fm_str), data=mydata)
                    return(car::Anova(fit)['Case', 1])})
imp_res = cbind(G_list2[keep_me, 1:2], chisq)
colnames(imp_res)[3] = 'F'  # quick hack to work with function

adhd_imp_camera = get_enrich_order2( imp_res, t2 )
c5_imp_camera = get_enrich_order2( imp_res, c5_all)
dis_imp_camera = get_enrich_order2( imp_res, disorders)
dev_imp_camera = get_enrich_order2( imp_res, genes_unique )
```

Now it's a matter of computing the intersection of the c5 results, because the
other ones are easy to visualize. For example:

```r
intersect(rownames(c5_dream_cameraDC)[c5_dream_cameraDC$PValue<.01],
          rownames(c5_imp_camera)[c5_imp_camera$PValue<.01])
intersect(resDC[resDC$P.Value<.05, 'hgnc_symbol'],
          imp_res[imp_res[,'Pr(>Chisq)']<.05, 'hgnc_symbol'])
```

## Refining all this

OK, now that the general framework is in place, let's start refining the
results. First, let's clean up the WNH thresholds. Then, we can play with
running adults and kids only, or even just males.

```r
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
gf = read.csv('~/data/expression_impute/gf_1119_09092020.csv')
imp_data = merge(imp_data, gf, by.x='IID', by.y='Subject.Code...Subjects', all.x=F, all.y=F)
```

Let's look for new thresholds:

```r
library(ggplot2)
ggplot(imp_data, aes(x=PC01, y=PC02, col=POP_CODE)) + geom_point()
```

![](images/2020-09-08-13-52-20.png)

I'm going to go even more strict here and go for 0-0:

```r
grex_vars = colnames(imp_data)[grepl(colnames(imp_data), pattern='^ENS')]

keep_me = imp_data$PC01 < 0 & imp_data$PC02 > 0
imp_data = imp_data[keep_me, ]
library(caret)

set.seed(42)
pp_order = c('zv', 'nzv')
pp = preProcess(imp_data[, grex_vars], method = pp_order)
X = predict(pp, imp_data[, grex_vars])
grex_vars = colnames(X)
imp_data = imp_data[, !grepl(colnames(imp_data), pattern='^ENS')]
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
colnames(X) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
grex_vars = id_num
X = t(X[, !dups])

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(X) %in% G_list$ensembl_gene_id
X = X[imnamed, ]
grex_vars = grex_vars[imnamed]
G_list2 = merge(G_list, X, by.x=1, by.y=0)
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
G_list2 = G_list2[imautosome, ]
X = G_list2[, 4:ncol(G_list2)]
rownames(X) = G_list2$ensembl_gene_id
grex_vars = G_list2$ensembl_gene_id

library(doParallel)
ncores = 2 #as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
cl = makeCluster(ncores)
nzeros = rowSums(X==0)
keep_me = nzeros < (ncol(X)/2)
good_grex = grex_vars[keep_me]
mydata = cbind(imp_data, t(X))
clusterExport(cl, c("mydata"))
chisq = parLapply(cl, good_grex,
                  function(x) {
                    fm_str = sprintf('%s ~ Diagnosis + age_base + Sex...Subjects + sample_type + PC01 + PC02 + PC03 + PC04 + PC05 + (1|FAMID)', x)
                    fit = lme4::lmer(as.formula(fm_str), data=mydata)
                    return(car::Anova(fit)['Diagnosis', c(1,3)])})
imp_res = cbind(G_list2[keep_me, 1:2], do.call(rbind, chisq))
colnames(imp_res)[3] = 'F'  # quick hack to work with function

adhd_imp_camera = get_enrich_order2( imp_res, t2 )
c5_imp_camera = get_enrich_order2( imp_res, c5_all)
dis_imp_camera = get_enrich_order2( imp_res, disorders)
dev_imp_camera = get_enrich_order2( imp_res, genes_unique )
sum(c5_imp_camera$FDR<.05)
sum(c5_imp_camera$PValue<.01)
```

![](images/2020-09-09-10-27-58.png)
![](images/2020-09-09-10-28-30.png)

### Males only

```r
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
gf = read.csv('~/data/expression_impute/gf_1119_09092020.csv')
imp_data = merge(imp_data, gf, by.x='IID', by.y='Subject.Code...Subjects', all.x=F, all.y=F)

imp_data = imp_data[imp_data$Sex...Subjects=='Male',]
grex_vars = colnames(imp_data)[grepl(colnames(imp_data), pattern='^ENS')]
keep_me = imp_data$PC01 < 0 & imp_data$PC02 > 0
imp_data = imp_data[keep_me, ]
library(caret)
set.seed(42)
pp_order = c('zv', 'nzv')
pp = preProcess(imp_data[, grex_vars], method = pp_order)
X = predict(pp, imp_data[, grex_vars])
grex_vars = colnames(X)
imp_data = imp_data[, !grepl(colnames(imp_data), pattern='^ENS')]
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
colnames(X) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
grex_vars = id_num
X = t(X[, !dups])

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(X) %in% G_list$ensembl_gene_id
X = X[imnamed, ]
grex_vars = grex_vars[imnamed]
G_list2 = merge(G_list, X, by.x=1, by.y=0)
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
G_list2 = G_list2[imautosome, ]
X = G_list2[, 4:ncol(G_list2)]
rownames(X) = G_list2$ensembl_gene_id
grex_vars = G_list2$ensembl_gene_id

library(doParallel)
ncores = 2 #as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
cl = makeCluster(ncores)
nzeros = rowSums(X==0)
keep_me = nzeros < (ncol(X)/2)
good_grex = grex_vars[keep_me]
mydata = cbind(imp_data, t(X))
clusterExport(cl, c("mydata"))
chisq = parLapply(cl, good_grex,
                  function(x) {
                    fm_str = sprintf('%s ~ Diagnosis + age_base + sample_type + PC01 + PC02 + PC03 + PC04 + PC05 + (1|FAMID)', x)
                    fit = lme4::lmer(as.formula(fm_str), data=mydata)
                    return(car::Anova(fit)['Diagnosis', c(1,3)])})
imp_res = cbind(G_list2[keep_me, 1:2], do.call(rbind, chisq))
colnames(imp_res)[3] = 'F'  # quick hack to work with function

adhd_imp_camera = get_enrich_order2( imp_res, t2 )
c5_imp_camera = get_enrich_order2( imp_res, c5_all)
dis_imp_camera = get_enrich_order2( imp_res, disorders)
dev_imp_camera = get_enrich_order2( imp_res, genes_unique )
sum(c5_imp_camera$FDR<.05)
sum(c5_imp_camera$PValue<.01)
```

![](images/2020-09-09-10-29-47.png)
![](images/2020-09-09-10-30-20.png)

### Adults only

```r
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
gf = read.csv('~/data/expression_impute/gf_1119_09092020.csv')
imp_data = merge(imp_data, gf, by.x='IID', by.y='Subject.Code...Subjects', all.x=F, all.y=F)

imp_data = imp_data[imp_data$age_group=='adult',]
grex_vars = colnames(imp_data)[grepl(colnames(imp_data), pattern='^ENS')]
keep_me = imp_data$PC01 < 0 & imp_data$PC02 > 0
imp_data = imp_data[keep_me, ]
library(caret)
set.seed(42)
pp_order = c('zv', 'nzv')
pp = preProcess(imp_data[, grex_vars], method = pp_order)
X = predict(pp, imp_data[, grex_vars])
grex_vars = colnames(X)
imp_data = imp_data[, !grepl(colnames(imp_data), pattern='^ENS')]
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
colnames(X) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
grex_vars = id_num
X = t(X[, !dups])

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(X) %in% G_list$ensembl_gene_id
X = X[imnamed, ]
grex_vars = grex_vars[imnamed]
G_list2 = merge(G_list, X, by.x=1, by.y=0)
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
G_list2 = G_list2[imautosome, ]
X = G_list2[, 4:ncol(G_list2)]
rownames(X) = G_list2$ensembl_gene_id
grex_vars = G_list2$ensembl_gene_id

library(doParallel)
ncores = 2 #as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
cl = makeCluster(ncores)
nzeros = rowSums(X==0)
keep_me = nzeros < (ncol(X)/2)
good_grex = grex_vars[keep_me]
mydata = cbind(imp_data, t(X))
clusterExport(cl, c("mydata"))
chisq = parLapply(cl, good_grex,
                  function(x) {
                    fm_str = sprintf('%s ~ Diagnosis + age_base + Sex...Subjects + sample_type + PC01 + PC02 + PC03 + PC04 + PC05 + (1|FAMID)', x)
                    fit = lme4::lmer(as.formula(fm_str), data=mydata)
                    return(car::Anova(fit)['Diagnosis', c(1,3)])})
imp_res = cbind(G_list2[keep_me, 1:2], do.call(rbind, chisq))
colnames(imp_res)[3] = 'F'  # quick hack to work with function

adhd_imp_camera = get_enrich_order2( imp_res, t2 )
c5_imp_camera = get_enrich_order2( imp_res, c5_all)
dis_imp_camera = get_enrich_order2( imp_res, disorders)
dev_imp_camera = get_enrich_order2( imp_res, genes_unique )
sum(c5_imp_camera$FDR<.05)
sum(c5_imp_camera$PValue<.01)
```

![](images/2020-09-09-10-49-37.png)
![](images/2020-09-09-10-50-06.png)

### Adult-man only

```r
a = readRDS('~/data/expression_impute/results/NCR_v3_ACC_1KG_mashr.rds')
iid2 = sapply(a$IID, function(x) strsplit(x, '_')[[1]][2])
a$IID = iid2
pcs = read.csv('~/data/expression_impute/pop_pcs.csv')
imp_data = merge(a, pcs, by='IID', all.x=F, all.y=F)
gf = read.csv('~/data/expression_impute/gf_1119_09092020.csv')
imp_data = merge(imp_data, gf, by.x='IID', by.y='Subject.Code...Subjects', all.x=F, all.y=F)

imp_data = imp_data[imp_data$Sex...Subjects=='Male',]
imp_data = imp_data[imp_data$age_group=='adult',]
grex_vars = colnames(imp_data)[grepl(colnames(imp_data), pattern='^ENS')]
keep_me = imp_data$PC01 < 0 & imp_data$PC02 > 0
imp_data = imp_data[keep_me, ]
library(caret)
set.seed(42)
pp_order = c('zv', 'nzv')
pp = preProcess(imp_data[, grex_vars], method = pp_order)
X = predict(pp, imp_data[, grex_vars])
grex_vars = colnames(X)
imp_data = imp_data[, !grepl(colnames(imp_data), pattern='^ENS')]
id_num = sapply(grex_vars, function(x) strsplit(x=x, split='\\.')[[1]][1])
colnames(X) = id_num
dups = duplicated(id_num)
id_num = id_num[!dups]
grex_vars = id_num
X = t(X[, !dups])

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
G_list0 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id",
                 "hgnc_symbol", "chromosome_name"),values=id_num,mart= mart)
G_list <- G_list0[!is.na(G_list0$hgnc_symbol),]
G_list = G_list[G_list$hgnc_symbol!='',]
G_list <- G_list[!duplicated(G_list$ensembl_gene_id),]
imnamed = rownames(X) %in% G_list$ensembl_gene_id
X = X[imnamed, ]
grex_vars = grex_vars[imnamed]
G_list2 = merge(G_list, X, by.x=1, by.y=0)
imautosome = which(G_list2$chromosome_name != 'X' &
                   G_list2$chromosome_name != 'Y' &
                   G_list2$chromosome_name != 'MT')
G_list2 = G_list2[imautosome, ]
X = G_list2[, 4:ncol(G_list2)]
rownames(X) = G_list2$ensembl_gene_id
grex_vars = G_list2$ensembl_gene_id

library(doParallel)
ncores = 2 #as.numeric(Sys.getenv('SLURM_CPUS_PER_TASK'))
cl = makeCluster(ncores)
nzeros = rowSums(X==0)
keep_me = nzeros < (ncol(X)/2)
good_grex = grex_vars[keep_me]
mydata = cbind(imp_data, t(X))
clusterExport(cl, c("mydata"))
chisq = parLapply(cl, good_grex,
                  function(x) {
                    fm_str = sprintf('%s ~ Diagnosis + age_base + sample_type + PC01 + PC02 + PC03 + PC04 + PC05 + (1|FAMID)', x)
                    fit = lme4::lmer(as.formula(fm_str), data=mydata)
                    return(car::Anova(fit)['Diagnosis', c(1,3)])})
imp_res = cbind(G_list2[keep_me, 1:2], do.call(rbind, chisq))
colnames(imp_res)[3] = 'F'  # quick hack to work with function

adhd_imp_camera = get_enrich_order2( imp_res, t2 )
c5_imp_camera = get_enrich_order2( imp_res, c5_all)
dis_imp_camera = get_enrich_order2( imp_res, disorders)
dev_imp_camera = get_enrich_order2( imp_res, genes_unique )

sum(c5_imp_camera$FDR<.05)
sum(c5_imp_camera$PValue<.01)
```

![](images/2020-09-09-10-50-31.png)
![](images/2020-09-09-10-51-01.png)

For comparison, these are the ACC post-mortem results:

![](images/2020-09-09-11-05-53.png)
![](images/2020-09-09-11-06-51.png)

Is there any overlap in gene sets or even genes?

All:
![](images/2020-09-09-11-28-30.png)

Males:
![](images/2020-09-09-11-29-06.png)
![](images/2020-09-09-11-29-25.png)

Adults:
![](images/2020-09-09-11-12-34.png)
![](images/2020-09-09-11-19-56.png)

Adult-male:
![](images/2020-09-09-11-13-33.png)
![](images/2020-09-09-11-13-50.png)
![](images/2020-09-09-11-14-05.png)
![](images/2020-09-09-11-21-04.png)

Here are some number for the Adult-male WNH sample:

![](images/2020-09-09-11-52-07.png)

# 2020-12-16 17:01:15

It looks like the ABCD imputation never finished. I can see that after 10h it
only processed 4 chromossomes. And it uses a single thread, less than 1Gb of
memory. So, it's best to just sbatch it so we can give it more time than what we
can use in a sinteractive session.

So, put this in ~/data/expression_impute/ABCD_ACC_MASHR.sh:

```bash
#!/bin/bash

#SBATCH -o ~/data/expression_impute/ABCD_ACC_MASHR.out
#SBATCH -e ~/data/expression_impute/ABCD_ACC_MASHR.err
module load python
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate imlabtools

mydir=~/data/expression_impute;
python3 $mydir/MetaXcan-master/software/Predict.py \
    --model_db_path $mydir/eqtl/mashr/mashr_Brain_Anterior_cingulate_cortex_BA24.db \
    --vcf_genotypes /data/NCR_SBRB/ABCD/v201/1KG/chr*.dose.vcf.gz \
    --vcf_mode imputed \
    --prediction_output $mydir/results/ABCD_v201_ACC_predict_1KG_mashr.txt \
    --prediction_summary_output $mydir/results/ABCD_v201_ACC_summary_1KG_mashr.txt \
    --verbosity 9 --throw --model_db_snp_key varID \
    --on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" \
    --liftover $mydir/hg19ToHg38.over.chain.gz
```

Then we can do:

```bash
sbatch --mem=8g --time 6-00:00:00 ~/data/expression_impute/ABCD_ACC_MASHR.sh
```

which allocates 2 CPU and 8Gb... plenty.

And repeat the steps above for Caudate and then the 2 EN models.

# TODO
 * try using the elastic net models
 * play with zero threshold
 * look into epiXcan
 * need some sort of permutation test to evaluate intersections?
 * blood samples only?
 * figure out those unknown samples
 * use entire population?
 * what's going on with caudate?
 * what happens to main results when we add all races and women?

# From meeting with Philip
 * remitted should be included as adhd. Typical are really typical. Maybe just
   do ever_ADHD?
 * narrow to age range that's comparable to postmortem
 * only one per family (even for extended family, for healthy, no family history)
 * look at cognition as well
