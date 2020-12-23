# 2020-11-30 21:01:06

The other idea is to see how PRS correlates with our post-mortem results,
similar to what they did in:

https://science.sciencemag.org/content/sci/362/6420/eaat8127.full.pdf?casa_token=fS2Ibw-5XBAAAAAA:VqdYfFH0YZnXUoU9w0EgxTl7y1xzqIFOvHmeJFcQOcsUFkx7kHTEj-3VNSCO3ltMhMbybsPqJCn5mfs


Kwangmi says the PM QCed data for analysis are
Shaw/Family_genotyping/Shaw01_2020/QCed for analysis.

But we always imputed it first, before computing PRS. So let's do that as well.
Following the ENIGMA pipeline:

**NOTE that this genotyping wave is hg38! So, I had to lift over form hg38 to 19
before running shapeit because I couldn't find a shapeit reference for hg38. It
shouldn't matter in the end as the imputed result is always in the reference
inputation, which is hg19 in 1KG Phase3** 

```bash
# bw
cd ~/data/post_mortem/genotyping
module load plink
# Filter SNPs out from your dataset which do not meet Quality Control criteria
# (Minor Allele Frequency < 0.01; Genotype Call Rate < 95%; Hardy­Weinberg
# Equilibrium < 1x10­6)
export datafileraw=Shaw01_2020_QCed_autosome
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb \
      --make-bed --out ${datafileraw}_filtered
export datafile=${datafileraw}_filtered
awk '{ if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig" ; else print $2 ;}' $datafile.bim | grep ambig | awk '{print $1}' > ambig.list;
plink --bfile $datafile --exclude ambig.list --make-founders --out lastQC \
    --maf 0.01 --hwe 0.000001 --make-bed --noweb
# awk '{print $2, $1":"$4}' lastQC.bim > updateSNPs.txt
plink --bfile lastQC --noweb --list-duplicate-vars --out lastQC
# plink --bfile lastQCb37 --exclude lastQCb37.dupvar --out lastQCb37_noduplicates \
    # --make-bed --noweb
plink --bfile lastQC --exclude lastQC.dupvar --out lastQC_noduplicates \
    --make-bed --noweb
module load vcftools
module load crossmap
for i in {1..22}; do
    plink --bfile lastQC_noduplicates --chr $i --recode vcf --out PM_chr"$i";
    # liftover
    crossmap vcf hg38ToHg19.over.chain PM_chr${i}.vcf \
        /fdb/GATK_resource_bundle/hg19/ucsc.hg19.fasta \
        PM_chr${i}.hg19.vcf;
    vcf-sort PM_chr"$i".hg19.vcf | bgzip -c > PM_chr"$i".hg19.vcf.gz
done
# check the strand using shape it before uploading.
module load shapeit/2.r904
refdir=/fdb/impute2/1000Genomes_Phase3_integrated_haplotypes_Oct2014/1000GP_Phase3/
for c in {1..22}; do
    shapeit -check -T 16 -V PM_chr${c}.hg19.vcf.gz \
        --input-ref $refdir/1000GP_Phase3_chr${c}.hap.gz \
        $refdir/1000GP_Phase3_chr${c}.legend.gz $refdir/1000GP_Phase3.sample \
        --output-log chr${c}.alignments;
done
# format the files:
for c in {1..22}; do
    grep Strand chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> flip_snps.txt;
    grep Missing chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> missing_snps.txt;
done
# I'm having some issues with duplicate ID that even the list-duplicate command
# in ENIGMA's protocol is not finding, because they have different alleles. So,
# let's remove them completely from the analysis, before we flip it using ShapeIt results:
module load plink
plink --bfile lastQCb37 --write-snplist --out all_snps
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink --bfile lastQCb37 --exclude duplicated_snps.snplist --make-bed --out lastQCb37_noduplicates
# flip and remove all bad ids
plink --bfile lastQCb37_noduplicates --flip flip_snps.txt \
    --exclude missing_snps.txt --make-bed --out lastQCb37_noduplicates_flipped
#reconstruct the VCFs as above to send it to the imputation server.
module load vcftools
for i in {1..22}; do
    plink --bfile lastQCb37_noduplicates_flipped --chr ${i} --recode-vcf \
        --out PM_chr${i}_flipped;
    vcf-sort PM_chr"$i"_flipped.vcf | bgzip -c > PM_chr"$i"_flipped.vcf.gz
done
```

Then we send it to the imputation server. This didn't work... maybe I'll try
sending Kwangmi's file directly, without any cleaning?

```bash
# bw
cd ~/data/post_mortem/genotyping
module load vcftools
module load plink
for i in {1..22}; do
    plink --bfile Shaw01_2020_QCed_autosome --chr ${i} --recode-vcf \
        --out PMraw_chr${i};
    vcf-sort PMraw_chr"$i".vcf | bgzip -c > PMraw_chr"$i".vcf.gz
done
```

That didn't work. More than 100 obvious strand flips detected. Let's run shape
it then.

```bash
cd ~/data/post_mortem/genotyping
module load plink
export datafileraw=Shaw01_2020_QCed_autosome
plink --bfile $datafileraw --hwe 1e-6 --geno 0.05 --maf 0.01 --noweb \
      --make-bed --out ${datafileraw}_filtered
export datafile=${datafileraw}_filtered
awk '{print $2, $1":"$4}' $datafile.bim > updateSNPs.txt
plink --bfile $datafile --update-name updateSNPs.txt --make-bed \
    --out ${datafile}_ren --noweb --list-duplicate-vars
plink --bfile ${datafile}_ren --exclude ${datafile}_ren.dupvar \
    --out ${datafile}_ren_nodups --make-bed --noweb
for i in {1..22}; do
    plink --bfile ${datafile}_ren_nodups --chr ${i} --recode-vcf \
        --out PMraw_chr${i};
    vcf-sort PMraw_chr"$i".vcf | bgzip -c > PMraw_chr"$i".vcf.gz
done

module load shapeit/2.r904
refdir=/fdb/impute2/1000Genomes_Phase3_integrated_haplotypes_Oct2014/1000GP_Phase3/
for c in {1..22}; do
    shapeit -check -T 16 -V PMraw_chr${c}.vcf.gz \
        --input-ref $refdir/1000GP_Phase3_chr${c}.hap.gz \
        $refdir/1000GP_Phase3_chr${c}.legend.gz $refdir/1000GP_Phase3.sample \
        --output-log chr${c}.alignments;
done
# format the files:
rm -rf flip_snps.txt missing_snps.txt;
for c in {1..22}; do
    grep Strand chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> flip_snps.txt;
    grep Missing chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> missing_snps.txt;
done

plink --bfile Shaw01_2020_QCed_autosome --write-snplist --out all_snps
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink --bfile Shaw01_2020_QCed_autosome --exclude duplicated_snps.snplist --make-bed --out Shaw01_2020_QCed_autosome_noduplicates
# flip and remove all bad ids
plink --bfile Shaw01_2020_QCed_autosome_noduplicates --flip flip_snps.txt \
    --exclude missing_snps.txt --make-bed --out Shaw01_2020_QCed_autosome_noduplicates_flipped
#reconstruct the VCFs as above to send it to the imputation server.
module load vcftools
for i in {1..22}; do
    plink --bfile Shaw01_2020_QCed_autosome_noduplicates_flipped --chr ${i} --recode-vcf \
        --out PMraw_chr${i}_flipped;
    vcf-sort PMraw_chr"$i"_flipped.vcf | bgzip -c > PMraw_chr"$i"_flipped.vcf.gz
done
```

I'm having all kinds of issues here. I asked Kwangmi for help... 

```
hey, can you help me with some imputation? I normally use the Michigan
Imputation server to impute our genotypes to 1KG, but first I need to flip the
strands in the Illumina data before uploading the files for imputation. That's
never a issue, but I'm having trouble with the post-mortem data because it's in
hg38. I normally use shapeit to identify the strands to be flipped, but I can't
find hg38 reference files for shapeit, only hg19. I then tried to liftover the
post-mortem data to hg19, which seemed to work, but shapeit gives an error
because there is massive misalignment between the lifted data and the reference
panel (that never happened when I played with data native to hg19). Any ideas?
```

Maybe try beagle for phasing? Kwangmi suggested MAF .05 as the population is
quite varied in the PM study Maybe calculate PRS on the non-imputed data? what's
the overlap?

# 2020-12-02 09:46:15

Maybe we don't need to impute at all for PRS? Let's see what's the overlap
between the original data and the GWAS base... no way. The base file has 8M
variants, and our original data only has 646K. It needs to be done in the
imputed data...

While I'm reading about strands, here are some useful material:

http://onetipperday.sterding.com/2015/12/snp-allele-coding-schema.html

https://www.illumina.com/documents/products/technotes/technote_topbot.pdf

https://support.illumina.com/bulletins/2017/06/how-to-interpret-dna-strand-and-allele-information-for-infinium-.html

https://support.illumina.com/bulletins/2016/06/simple-guidelines-for-identifying-topbottom-topbot-strand-and-ab-allele.html

This finally worked! I used this tool:
https://github.com/seasky002002/Strandscript

```bash
# bw
cd ~/data/post_mortem/genotyping/
plink --bfile Shaw01_2020_QCed_autosome --out PM_test --recode
cd Strandscript-master
perl ./bin/step1-mismatch.pl -o test/ -n PM -g GRCh38 \
    -in ../InfiniumOmniExpressExome-8v1-4_A2.csv
perl ./bin/step2-flip.pl -o test/ -in test/new_PM.csv -map ../PM_test.map \
    -ped ../PM_test.ped -c 0.2
for i in {1..22}; do
    plink --file test/flipped_PM_test --chr ${i} --recode vcf \
        --out PMraw_chr${i} --real-ref-alleles;
    # replace # by chr#
    awk '{if($0 !~ /^#/) print "chr"$0; else print $0}' \
        PMraw_chr${i}.vcf > PMraw_chr${i}_ren.vcf
    vcf-sort PMraw_chr"$i"_ren.vcf | bgzip -c > PMraw_chr"$i"_ren.vcf.gz
done
```

Then I just ran it through the Imputation Server and it worked! Now the data is
in 1KG Phase 3, and I selected the multi-racial option. In any case, it's now in
hg19.

```bash
# bw
cd ~/data/post_mortem/genotyping/1KG
for f in `/bin/ls *zip`; do unzip -P jXyw3G7FePopCi $f; done
```

Time to create the PLINK version.

```bash
# sinteractive
for c in {1..22}; do
    plink --vcf chr${c}.dose.vcf.gz --double-id --make-bed --out chr${c};
done
rm -rf merge_list.txt;
for c in {2..22}; do
    echo "chr${c}" >> merge_list.txt;
done
plink --bfile chr1 --merge-list merge_list.txt --make-bed --out PM_1KG
# cleaning based on imputation stats
for c in {1..22}; do
    echo $c;
    zcat chr${c}.info.gz | awk '{ print $1,$5,$7 }' - >> r2s.txt;
done
awk '$2 > .01 && $3 > .9 { print }' r2s.txt > rsids_MAFbtp01_rsbtp9.txt
plink --bfile PM_1KG --extract rsids_MAFbtp01_rsbtp9.txt --geno .05 \
    --make-bed --out PM_1KG_genop05MAFbtp01rsbtp9
```

Let's use the HRC variables just for renaming, as I can't find a similar file
for 1KG.

```bash
wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz
awk '{print $1 ":" $2 " " $3}' HRC.r1-1.GRCh37.wgs.mac5.sites.tab | \
    tail -n +2 > rename_ids.txt 
cat rename_ids.txt | sort -u -k 1,1 | uniq > unique_rename_ids.txt
# that took a while...
# we need to do two renaming. First, remove the variant from the name column
cut -f 2 PM_1KG_genop05MAFbtp01rsbtp9.bim > tmp_name.txt;
cut -d":" -f 1,2 tmp_name.txt > tmp_name2.txt;
paste tmp_name.txt tmp_name2.txt > update_snps1.txt
plink --bfile PM_1KG_genop05MAFbtp01rsbtp9 --update-name update_snps1.txt \
    --make-bed --out tmp
plink --bfile tmp --write-snplist --out all_snps
cat all_snps.snplist | sort | uniq -d > duplicated_snps.snplist
plink --bfile tmp --exclude duplicated_snps.snplist --out tmp_nodups --make-bed --noweb
# then replace positional by rs ids
plink --bfile tmp_nodups --update-name unique_rename_ids.txt \
    --make-bed --out PM_1KG_genop05MAFbtp01rsbtp9_renamed
```

```bash
module load R
Rscript /data/NCR_SBRB/software/PRSice_2.2.5/PRSice.R  \
    --prsice /data/NCR_SBRB/software/PRSice_2.2.5/PRSice_linux \
    --base ~/pgc2017/adhd_jun2017  \
    --target ~/data/post_mortem/genotyping/1KG/PM_1KG_genop05MAFbtp01rsbtp9_renamed \
    --all-score \
    --lower 5e-08 --upper .5 --interval 5e-05 \
    --no-regress \
    --out PM_1KG_PRS_adhd_jun2017

Rscript /data/NCR_SBRB/software/PRSice_2.2.5/PRSice.R  \
    --prsice /data/NCR_SBRB/software/PRSice_2.2.5/PRSice_linux \
    --base ~/pgc2017/adhd_eur_jun2017  \
    --target ~/data/post_mortem/genotyping/1KG/PM_1KG_genop05MAFbtp01rsbtp9_renamed \
    --all-score --stat OR \
    --lower 5e-08 --upper .5 --interval 5e-05 \
    --no-regress \
    --out PM_1KG_PRS_adhd_eur_jun2017
```

Now it's just a matter of putting it all together for analysis!

# 2020-12-03 05:54:06

```r
# this takes a while because we're reading in TXT files!
my_dir = '~/data/post_mortem/genotyping/1KG/'
a = read.table(sprintf('%s/PM_1KG_PRS_adhd_jun2017.all.score', my_dir),
               header=1)
b = read.table(sprintf('%s/PM_1KG_PRS_adhd_eur_jun2017.all.score', my_dir),
                header=1)

# keep only some of the PRs columns that were created
mycols = c('IID', 'X0.00010005', 'X0.00100005', 'X0.01', 'X0.1',
            'X5.005e.05', 'X0.00050005', 'X0.00500005', 'X0.0500001',
            'X0.5', 'X0.4', 'X0.3', 'X0.2')
new_names = c('IID', sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                              .5, .4, .3, .2),
                            function(x) sprintf('ADHD_PRS%f', x)))
af = a[, mycols]
colnames(af) = new_names
bf = b[, mycols]
new_names = gsub('ADHD', x=new_names, 'ADHDeur')
colnames(bf) = new_names

m = merge(af, bf, by='IID')

write.csv(m, file='~/data/post_mortem/genotyping/1KG/merged_PM_1KG_PRS_12032020.csv', row.names=F)
```

Now, when we run this analysis, we can do it in different ways. Primarily, we
can do the same thing the Science paper did and just try to evaluate it against
the first PC. But we can also just replace Diagnosis by PRS in the main
regression and see if anything comes out of that.

Within those two options, we can use the WNH PRS on the WNH samples only, then
the general on everyone, then the most appropriate PRS for each group.

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
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
# we're down from 60K to 38K samples by only looking at the ones with hgnc symbol. We might be losing too much here, so it's a step to reconsider in the future

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

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])

data2 = cbind(data, mydata)
# we don't residualize Diagnosis!

form = ~ PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes, design, plot=FALSE)
fit <- lmFit(vobj, design)
fit2 <- eBayes( fit )
resids = residuals(fit2, genes)
```

Let's first break down the 3 different sets of PRS we will try:

```r
fname = '~/data/post_mortem/genotyping/1KG/merged_PM_1KG_PRS_12032020.csv'
prs = read.csv(fname)
prs$hbcc_brain_id = sapply(prs$IID,
                          function(x) {
                              br = strsplit(x, '_')[[1]][2];
                              as.numeric(gsub(br, pattern='BR',
                                              replacement=''))})
imWNH = data$C1 > 0 & data$C2 < -.075
wnh_brains = data[which(imWNH),]$hbcc_brain_id

# using whole population PRS
m = merge(data, prs, by='hbcc_brain_id')
m = m[, 1:63]
colnames(m)[52:63] = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                              .5, .4, .3, .2),
                            function(x) sprintf('PRS%f', x))
data_whole = m

# WNH samples only
m = merge(data, prs, by='hbcc_brain_id')
m = m[, c(1:51, 64:75)]
colnames(m)[52:63] = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                              .5, .4, .3, .2),
                            function(x) sprintf('PRS%f', x))
keep_me = m$hbcc_brain_id %in% wnh_brains
data_wnh = m[keep_me, ]

# using the most appropriate PRS
m = merge(data, prs, by='hbcc_brain_id')
prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
m[, prs_names] = NA
keep_me = m$hbcc_brain_id %in% wnh_brains
m[keep_me, prs_names] = m[keep_me, 64:75]
m[!keep_me, prs_names] = m[!keep_me, 52:63]
data_app = m[, c(1:51, 76:87)]
```

And of course we'll need to try this for all possible PRS values. Let's try the
PCA approach first, as it's easier to quantify the results.

Now that I'm re-reading the paper, they only took the first PC in Fig 3 to
compare the different modalities. For PRS, they used the formula in the
supplement. So, let's just go with that. I can use the formula to obtain the
best PRS threshold, and start examining that, but I'll likely eventually run all
of them just for comparison.

Actually, I think it's both. To come up with a single R2, you need one summary
value, so in that case it has to be the first PC.

```r
library(fmsb)
library(lmtest)
mod.baseline = glm(Diagnosis ~ Sex + Age, family=binomial, data=data_app)
mod.full = glm(Diagnosis ~ PRS0.000100 + Sex + Age, family=binomial, data=data_app)
adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
prs.significance = lrtest(mod.baseline, mod.full)
```

So, if we are goig to try all of them:

```r
mydata = data_app

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
for (prs in prs_names) {
    # fm_root = 'Diagnosis ~ %s Sex + Age'
    fm_root = 'Diagnosis ~ %s Sex + Age + C1 + C2 + C3 + C4 + C5'
    fm_str = sprintf(fm_root, '')
    mod.baseline = glm(as.formula(fm_str), family=binomial, data=mydata)
    fm_str = sprintf(fm_root, sprintf('%s +', prs))
    mod.full = glm(as.formula(fm_str), family=binomial, data=mydata)
    adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
    prs.significance = lrtest(mod.baseline, mod.full)
    myp = prs.significance$"Pr(>Chisq)"[2] 
    cat(sprintf('%s: pval = %.3f\n', fm_str, myp))
}
```

```
(whole)
Diagnosis ~ PRS0.000100 + Sex + Age: pval = 0.980
Diagnosis ~ PRS0.001000 + Sex + Age: pval = 0.424
Diagnosis ~ PRS0.010000 + Sex + Age: pval = 0.226
Diagnosis ~ PRS0.100000 + Sex + Age: pval = 0.706
Diagnosis ~ PRS0.000050 + Sex + Age: pval = 0.933
Diagnosis ~ PRS0.000500 + Sex + Age: pval = 0.289
Diagnosis ~ PRS0.005000 + Sex + Age: pval = 0.285
Diagnosis ~ PRS0.050000 + Sex + Age: pval = 0.209
Diagnosis ~ PRS0.500000 + Sex + Age: pval = 0.715
Diagnosis ~ PRS0.400000 + Sex + Age: pval = 0.834
Diagnosis ~ PRS0.300000 + Sex + Age: pval = 0.753
Diagnosis ~ PRS0.200000 + Sex + Age: pval = 0.501

(appropriate)
Diagnosis ~ PRS0.000100 + Sex + Age: pval = 0.519
Diagnosis ~ PRS0.001000 + Sex + Age: pval = 0.069
Diagnosis ~ PRS0.010000 + Sex + Age: pval = 0.150
Diagnosis ~ PRS0.100000 + Sex + Age: pval = 0.273
Diagnosis ~ PRS0.000050 + Sex + Age: pval = 0.062
Diagnosis ~ PRS0.000500 + Sex + Age: pval = 0.178
Diagnosis ~ PRS0.005000 + Sex + Age: pval = 0.041
Diagnosis ~ PRS0.050000 + Sex + Age: pval = 0.107
Diagnosis ~ PRS0.500000 + Sex + Age: pval = 0.681
Diagnosis ~ PRS0.400000 + Sex + Age: pval = 0.601
Diagnosis ~ PRS0.300000 + Sex + Age: pval = 0.837
Diagnosis ~ PRS0.200000 + Sex + Age: pval = 0.842

(WNH)
Diagnosis ~ PRS0.000100 + Sex + Age: pval = 0.670
Diagnosis ~ PRS0.001000 + Sex + Age: pval = 0.910
Diagnosis ~ PRS0.010000 + Sex + Age: pval = 0.802
Diagnosis ~ PRS0.100000 + Sex + Age: pval = 0.744
Diagnosis ~ PRS0.000050 + Sex + Age: pval = 0.624
Diagnosis ~ PRS0.000500 + Sex + Age: pval = 0.854
Diagnosis ~ PRS0.005000 + Sex + Age: pval = 0.586
Diagnosis ~ PRS0.050000 + Sex + Age: pval = 0.618
Diagnosis ~ PRS0.500000 + Sex + Age: pval = 0.968
Diagnosis ~ PRS0.400000 + Sex + Age: pval = 0.934
Diagnosis ~ PRS0.300000 + Sex + Age: pval = 0.999
Diagnosis ~ PRS0.200000 + Sex + Age: pval = 0.941
```

It looks like using the "Appropriate" dataset has better results in predicting
Diagnosis. But there are other cofounders there too, such as comorbidity and
substance abuse, that are not being taken into this analysis. Let's at least add
the population PCs just to mirror the Science paper a bit more and see if they
have any difference. Just the first 5:

```
(whole)
Diagnosis ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.435
Diagnosis ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.380
Diagnosis ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.117
Diagnosis ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.092
Diagnosis ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.542
Diagnosis ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.585
Diagnosis ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.054
Diagnosis ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.194
Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.015
Diagnosis ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.018
Diagnosis ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.045
Diagnosis ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.138

(appropriate)
Diagnosis ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.676
Diagnosis ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.471
Diagnosis ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.140
Diagnosis ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.057
Diagnosis ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.405
Diagnosis ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.528
Diagnosis ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.013
Diagnosis ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.125
Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.004
Diagnosis ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.009
Diagnosis ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.019
Diagnosis ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.065

(WNH)
Diagnosis ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.356
Diagnosis ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.208
Diagnosis ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.645
Diagnosis ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.244
Diagnosis ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.662
Diagnosis ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.276
Diagnosis ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.019
Diagnosis ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.201
Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.436
Diagnosis ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.537
Diagnosis ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.559
Diagnosis ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.493
```

That made a huge difference. Good thing I tried it. Didn't impact WNH results as
much, but that's fine.

```r
resids.pca = prcomp(t(resids), scale=TRUE)
print(summary(resids.pca)$importance[, 'PC1'][2])
pcs = resids.pca$x
rownames(pcs) = data$hbcc_brain_id

mydata = merge(pcs, data_app, by.x=0, by.y='hbcc_brain_id')
for (prs in prs_names) {
    fm_root = 'PC1 ~ %s Sex + Age'
    fm_root = 'PC1 ~ %s 1'
    # fm_root = 'PC1 ~ %s Sex + Age + C1 + C2 + C3 + C4 + C5'
    fm_str = sprintf(fm_root, '')
    mod.baseline = glm(as.formula(fm_str), family=gaussian, data=mydata)
    fm_str = sprintf(fm_root, sprintf('%s + ', prs))
    mod.full = glm(as.formula(fm_str), family=gaussian, data=mydata)
    adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
    prs.significance = lrtest(mod.baseline, mod.full)
    myp = prs.significance$"Pr(>Chisq)"[2] 
    cat(sprintf('%s: pval = %.3f\n', fm_str, myp))
}
```

```
(whole)
PC1 ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.794
PC1 ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.434
PC1 ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.543
PC1 ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.667
PC1 ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.081
PC1 ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.997
PC1 ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.912
PC1 ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.960
PC1 ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.295
PC1 ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.250
PC1 ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.437
PC1 ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.866

(appropriate)
PC1 ~ PRS0.000100 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.465
PC1 ~ PRS0.001000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.365
PC1 ~ PRS0.010000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.690
PC1 ~ PRS0.100000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.952
PC1 ~ PRS0.000050 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.370
PC1 ~ PRS0.000500 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.620
PC1 ~ PRS0.005000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.525
PC1 ~ PRS0.050000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.863
PC1 ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.471
PC1 ~ PRS0.400000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.387
PC1 ~ PRS0.300000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.528
PC1 ~ PRS0.200000 + Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.924
```

Not much going on there. Maybe there is something if I remove the PCs? They were
already somewhat removed to begin with...

```
(whole)
PC1 ~ PRS0.000100 + Sex + Age: pval = 0.659
PC1 ~ PRS0.001000 + Sex + Age: pval = 0.302
PC1 ~ PRS0.010000 + Sex + Age: pval = 0.169
PC1 ~ PRS0.100000 + Sex + Age: pval = 0.298
PC1 ~ PRS0.000050 + Sex + Age: pval = 0.094
PC1 ~ PRS0.000500 + Sex + Age: pval = 0.527
PC1 ~ PRS0.005000 + Sex + Age: pval = 0.375
PC1 ~ PRS0.050000 + Sex + Age: pval = 0.392
PC1 ~ PRS0.500000 + Sex + Age: pval = 0.787
PC1 ~ PRS0.400000 + Sex + Age: pval = 0.771
PC1 ~ PRS0.300000 + Sex + Age: pval = 0.890
PC1 ~ PRS0.200000 + Sex + Age: pval = 0.525

(appropriate)
PC1 ~ PRS0.000100 + Sex + Age: pval = 0.791
PC1 ~ PRS0.001000 + Sex + Age: pval = 0.271
PC1 ~ PRS0.010000 + Sex + Age: pval = 0.372
PC1 ~ PRS0.100000 + Sex + Age: pval = 0.671
PC1 ~ PRS0.000050 + Sex + Age: pval = 0.318
PC1 ~ PRS0.000500 + Sex + Age: pval = 0.349
PC1 ~ PRS0.005000 + Sex + Age: pval = 0.412
PC1 ~ PRS0.050000 + Sex + Age: pval = 0.210
PC1 ~ PRS0.500000 + Sex + Age: pval = 0.558
PC1 ~ PRS0.400000 + Sex + Age: pval = 0.654
PC1 ~ PRS0.300000 + Sex + Age: pval = 0.916
PC1 ~ PRS0.200000 + Sex + Age: pval = 0.612
```

Not really. But if that's the logic, then I should remove age and sex as well:

```
(whole)
PC1 ~ PRS0.000100 +  1: pval = 0.794
PC1 ~ PRS0.001000 +  1: pval = 0.353
PC1 ~ PRS0.010000 +  1: pval = 0.196
PC1 ~ PRS0.100000 +  1: pval = 0.386
PC1 ~ PRS0.000050 +  1: pval = 0.144
PC1 ~ PRS0.000500 +  1: pval = 0.637
PC1 ~ PRS0.005000 +  1: pval = 0.422
PC1 ~ PRS0.050000 +  1: pval = 0.506
PC1 ~ PRS0.500000 +  1: pval = 0.639
PC1 ~ PRS0.400000 +  1: pval = 0.618
PC1 ~ PRS0.300000 +  1: pval = 0.931
PC1 ~ PRS0.200000 +  1: pval = 0.678

(appropriate)
PC1 ~ PRS0.000100 +  1: pval = 0.644
PC1 ~ PRS0.001000 +  1: pval = 0.254
PC1 ~ PRS0.010000 +  1: pval = 0.349
PC1 ~ PRS0.100000 +  1: pval = 0.605
PC1 ~ PRS0.000050 +  1: pval = 0.288
PC1 ~ PRS0.000500 +  1: pval = 0.277
PC1 ~ PRS0.005000 +  1: pval = 0.393
PC1 ~ PRS0.050000 +  1: pval = 0.188
PC1 ~ PRS0.500000 +  1: pval = 0.701
PC1 ~ PRS0.400000 +  1: pval = 0.805
PC1 ~ PRS0.300000 +  1: pval = 0.753
PC1 ~ PRS0.200000 +  1: pval = 0.769
```

How about WNH?

```r
keep_me = data$hbcc_brain_id %in% wnh_brains
resids.pca = prcomp(t(resids[, keep_me]), scale=TRUE)
print(summary(resids.pca)$importance[, 'PC1'][2])
pcs = resids.pca$x
rownames(pcs) = data[keep_me,]$hbcc_brain_id

mydata = merge(pcs, data_wnh, by.x=0, by.y='hbcc_brain_id')
for (prs in prs_names) {
    # fm_root = 'PC1 ~ %s Sex + Age'
    # fm_root = 'PC1 ~ %s 1'
    fm_root = 'PC1 ~ %s Sex + Age + C1 + C2 + C3 + C4 + C5'
    fm_str = sprintf(fm_root, '')
    mod.baseline = glm(as.formula(fm_str), family=gaussian, data=mydata)
    fm_str = sprintf(fm_root, sprintf('%s + ', prs))
    mod.full = glm(as.formula(fm_str), family=gaussian, data=mydata)
    adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
    prs.significance = lrtest(mod.baseline, mod.full)
    myp = prs.significance$"Pr(>Chisq)"[2] 
    cat(sprintf('%s: pval = %.3f\n', fm_str, myp))
}
```

```
(WNH)
PC1 ~ PRS0.000100 +  1: pval = 0.544
PC1 ~ PRS0.001000 +  1: pval = 0.135
PC1 ~ PRS0.010000 +  1: pval = 0.273
PC1 ~ PRS0.100000 +  1: pval = 0.639
PC1 ~ PRS0.000050 +  1: pval = 0.539
PC1 ~ PRS0.000500 +  1: pval = 0.330
PC1 ~ PRS0.005000 +  1: pval = 0.191
PC1 ~ PRS0.050000 +  1: pval = 0.646
PC1 ~ PRS0.500000 +  1: pval = 0.761
PC1 ~ PRS0.400000 +  1: pval = 0.784
PC1 ~ PRS0.300000 +  1: pval = 0.800
PC1 ~ PRS0.200000 +  1: pval = 0.784

PC1 ~ PRS0.000100 +  Sex + Age: pval = 0.306
PC1 ~ PRS0.001000 +  Sex + Age: pval = 0.043
PC1 ~ PRS0.010000 +  Sex + Age: pval = 0.210
PC1 ~ PRS0.100000 +  Sex + Age: pval = 0.621
PC1 ~ PRS0.000050 +  Sex + Age: pval = 0.330
PC1 ~ PRS0.000500 +  Sex + Age: pval = 0.096
PC1 ~ PRS0.005000 +  Sex + Age: pval = 0.158
PC1 ~ PRS0.050000 +  Sex + Age: pval = 0.546
PC1 ~ PRS0.500000 +  Sex + Age: pval = 0.820
PC1 ~ PRS0.400000 +  Sex + Age: pval = 0.819
PC1 ~ PRS0.300000 +  Sex + Age: pval = 0.842
PC1 ~ PRS0.200000 +  Sex + Age: pval = 0.733

PC1 ~ PRS0.000100 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.814
PC1 ~ PRS0.001000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.111
PC1 ~ PRS0.010000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.875
PC1 ~ PRS0.100000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.221
PC1 ~ PRS0.000050 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.702
PC1 ~ PRS0.000500 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.121
PC1 ~ PRS0.005000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.811
PC1 ~ PRS0.050000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.271
PC1 ~ PRS0.500000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.034
PC1 ~ PRS0.400000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.040
PC1 ~ PRS0.300000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.025
PC1 ~ PRS0.200000 +  Sex + Age + C1 + C2 + C3 + C4 + C5: pval = 0.045
```

Nothing to write home about... maybe there is something there in WNH-only, but
that's only 30 subjects. If we're doing a WNH analysis, I think we'd need to do
it from scratch, including the initial removal of PCs?

Just because I'll wonder about this later, doing resids(fit) = resids(fit2), so
we don't need to worry about doing it in the results of lmFit or Bayes.

But let's check on the individual genes. Anything interesting there? Let's go
with our best PRS and then I can try other things:

```r
lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
rownames(mydata) = data$hbcc_brain_id

data2 = merge(data_app, mydata, by.x='hbcc_brain_id', by.y=0)
genes2 = genes[, data$hbcc_brain_id %in% data2$hbcc_brain_id]
form = ~ PRS0.500000 + PC1 + PC2 + PC7 + PC8 + PC9
design = model.matrix( form, data2)
vobj = voom( genes2, design, plot=FALSE)
prs.fit <- lmFit(vobj, design)
prs.fit2 <- eBayes( prs.fit )
res = topTable(prs.fit2, coef='PRS0.500000', number=Inf)
```

OK, now I have lots of results. What to do? Maybe first we should check if there
is a significant overlap between these genes and the PM_ACC genes?

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
for (t in c(.05, .01, .005, .001)) {
    prs_genes = res[res$P.Value < t, 'hgnc_symbol']
    dx_genes = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
    go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
    go.obj <- testGeneOverlap(go.obj)
    inter = intersect(prs_genes, dx_genes)
    pval = getPval(go.obj)
    cat(sprintf('t=%.3f, prs=%d, pm=%d, in=%d, p=%f\n', t,
                length(prs_genes), length(dx_genes), length(inter), pval))
}
```

This is a potentially cool result. 

```
t=0.050, prs=1504, pm=1335, in=202, p=0.000000
t=0.010, prs=305, pm=325, in=19, p=0.000004
t=0.005, prs=168, pm=165, in=10, p=0.000004
t=0.001, prs=42, pm=38, in=1, p=0.086516
```

So, these are quite cool. Anything in the caudate?

```r
myregion = 'Caudate'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
# remove obvious outlier (that's NOT caudate) labeled as ACC
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)

# at this point we have 55 samples for ACC
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
# we're down from 60K to 38K samples by only looking at the ones with hgnc symbol. We might be losing too much here, so it's a step to reconsider in the future

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

lcpm = cpm(genes, log=T)
set.seed(42)
lcpm.pca <- prcomp(t(lcpm), scale=TRUE)

library(nFactors)
eigs <- lcpm.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
mydata = data.frame(lcpm.pca$x[, keep_me])
rownames(mydata) = data$hbcc_brain_id

data2 = merge(data_app, mydata, by.x='hbcc_brain_id', by.y=0)
genes2 = genes[, data$hbcc_brain_id %in% data2$hbcc_brain_id]
form = ~ PRS0.500000 + PC1 + PC3 + PC5 + PC6 + PC8
design = model.matrix( form, data2)
vobj = voom( genes2, design, plot=FALSE)
prs.fit <- lmFit(vobj, design)
prs.fit2 <- eBayes( prs.fit )
res2 = topTable(prs.fit2, coef='PRS0.500000', number=Inf)
```

Then we run the same over-representation analysis:

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
for (t in c(.05, .01, .005, .001)) {
    prs_genes = res2[res2$P.Value < t, 'hgnc_symbol']
    dx_genes = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
    go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
    go.obj <- testGeneOverlap(go.obj)
    inter = intersect(prs_genes, dx_genes)
    pval = getPval(go.obj)
    cat(sprintf('t=%.3f, prs=%d, pm=%d, in=%d, p=%f\n', t,
                length(prs_genes), length(dx_genes), length(inter), pval))
}
```

And this is the Caudate... interesting:

```
t=0.050, prs=1339, pm=1063, in=72, p=0.860254
t=0.010, prs=264, pm=220, in=1, p=0.964251
t=0.005, prs=125, pm=117, in=0, p=1.000000
t=0.001, prs=27, pm=26, in=0, p=1.000000
```

# 2020-12-04 09:49:11

Might as well save the results for all PRS thresholds:

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC2 + PC7 + PC8 + PC9', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t, 'hgnc_symbol']
        dx_genes = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval = getPval(go.obj)
        this_res = c(p, t, length(prs_genes), length(dx_genes), length(inter),
                     pval)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('PRS', 'nomPvalThresh', 'PRsgenes', 'PMgenes',
                      'overlap', 'pval')
write.csv(all_res, file='~/data/post_mortem/all_acc_prs_overlap_results.csv',
          row.names=F)
```

And run the same thing for Caudate:

```r
all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC3 + PC5 + PC6 + PC8', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t, 'hgnc_symbol']
        dx_genes = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval = getPval(go.obj)
        this_res = c(p, t, length(prs_genes), length(dx_genes), length(inter),
                     pval)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('PRS', 'nomPvalThresh', 'PRsgenes', 'PMgenes',
                      'overlap', 'pval')
write.csv(all_res, file='~/data/post_mortem/all_caudate_prs_overlap_results.csv',
          row.names=F)
```

# 2020-12-22 20:46:32

Do these results change if we classify the genes between over and under
expressed? Or can we at leat quantify the overlap between over and under
expressed?

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC2 + PC7 + PC8 + PC9', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t & res$t > 0, 'hgnc_symbol']
        dx_genes = rnaseq_acc[rnaseq_acc$P.Value < t & rnaseq_acc$t > 0,
                              'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allUp = union(res[res$t > 0, 'hgnc_symbol'],
                      rnaseq_acc[rnaseq_acc$t > 0, 'hgnc_symbol'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(p, t, 'up', length(prs_genes), length(dx_genes), length(inter),
                     pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t & res$t < 0, 'hgnc_symbol']
        dx_genes = rnaseq_acc[rnaseq_acc$P.Value < t & rnaseq_acc$t < 0,
                              'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allDown = union(res[res$t < 0, 'hgnc_symbol'],
                      rnaseq_acc[rnaseq_acc$t < 0, 'hgnc_symbol'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(p, t, 'down', length(prs_genes), length(dx_genes), length(inter),
                     pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('PRS', 'nomPvalThresh', 'direction', 'PRsgenes', 'PMgenes',
                      'overlap', 'pvalWhole', 'pvalDirOnly')
write.csv(all_res, file='~/data/post_mortem/all_accUpDown_prs_overlap_results.csv',
          row.names=F)
```

The results weren't as conclusive as the no-direction results. Similarly, we run
this for the Caudate:

```r
library(GeneOverlap)
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
all_res = c()
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC3 + PC5 + PC6 + PC8', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t & res$t > 0, 'hgnc_symbol']
        dx_genes = rnaseq_caudate[rnaseq_caudate$P.Value < t & rnaseq_caudate$t > 0,
                              'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allUp = union(res[res$t > 0, 'hgnc_symbol'],
                      rnaseq_caudate[rnaseq_caudate$t > 0, 'hgnc_symbol'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allUp))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(p, t, 'up', length(prs_genes), length(dx_genes), length(inter),
                     pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
    for (t in c(.05, .01, .005, .001)) {
        prs_genes = res[res$P.Value < t & res$t < 0, 'hgnc_symbol']
        dx_genes = rnaseq_caudate[rnaseq_caudate$P.Value < t & rnaseq_caudate$t < 0,
                              'hgnc_symbol']
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=nrow(res))
        go.obj <- testGeneOverlap(go.obj)
        inter = intersect(prs_genes, dx_genes)
        pval1 = getPval(go.obj)
        allDown = union(res[res$t < 0, 'hgnc_symbol'],
                      rnaseq_caudate[rnaseq_caudate$t < 0, 'hgnc_symbol'])
        go.obj <- newGeneOverlap(prs_genes, dx_genes, genome.size=length(allDown))
        go.obj <- testGeneOverlap(go.obj)
        pval2 = getPval(go.obj)
        this_res = c(p, t, 'down', length(prs_genes), length(dx_genes), length(inter),
                     pval1, pval2)
        all_res = rbind(all_res, this_res)
    }
}
colnames(all_res) = c('PRS', 'nomPvalThresh', 'direction', 'PRsgenes', 'PMgenes',
                      'overlap', 'pvalWhole', 'pvalDirOnly')
write.csv(all_res,
          file='~/data/post_mortem/all_caudateUpDown_prs_overlap_results.csv',
          row.names=F)
```


# 2020-12-23 06:07:42

Let's then just make Venn diagrams of the significant results:

```r
library(GeneOverlap)
library(VennDiagram)

load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')

prs_names = sapply(c(.0001, .001, .01, .1, .00005, .0005, .005, .05,
                      .5, .4, .3, .2),
                   function(x) sprintf('PRS%f', x))
for (p in prs_names) {
    cat(p, '\n')
    form = as.formula(sprintf('~ %s + PC1 + PC2 + PC7 + PC8 + PC9', p))
    design = model.matrix( form, data2)
    vobj = voom( genes2, design, plot=FALSE)
    prs.fit <- lmFit(vobj, design)
    prs.fit2 <- eBayes( prs.fit )
    res = topTable(prs.fit2, coef=p, number=Inf)

    for (t in c(.05, .01, .005, .001)) {
        up_prs = res[res$P.Value < t & res$t > 0, 'hgnc_symbol']
        up_pm = rnaseq_acc[rnaseq_acc$P.Value < t & rnaseq_acc$t > 0,
                           'hgnc_symbol']
        main_str = sprintf('%s - %.3f - up', p, t)
        out_fname = sprintf('~/tmp/up_%s_%.3f.png', p, t)
        venn.plot = venn.diagram(list(PM = up_pm, PRS = up_prs),
                          euler.d=TRUE, fill=c('red','green'), main=main_str,
                          filename=out_fname, imagetype='png')
        
        down_prs = res[res$P.Value < t & res$t < 0, 'hgnc_symbol']
        down_pm = rnaseq_acc[rnaseq_acc$P.Value < t & rnaseq_acc$t < 0,
                           'hgnc_symbol']
        main_str = sprintf('%s - %.3f - down', p, t)
        out_fname = sprintf('~/tmp/down_%s_%.3f.png', p, t)
        venn.plot = venn.diagram(list(PM = down_pm, PRS = down_prs),
                          euler.d=TRUE, fill=c('red','green'), main=main_str,
                          filename=out_fname, imagetype='png')
    }
}
```

Not much of an overlap in these results, even if I invert the sign of the PRS in
case I coded it in the wrong direction. Actually, let's check this...

```r
library(fmsb)
library(lmtest)
mod.baseline = glm(Diagnosis ~ Sex + Age + C1 + C2 + C3 + C4 + C5, family=binomial, data=data_app)
mod.full = glm(Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + C3 + C4 + C5, family=binomial, data=data_app)
adjustedR2 = NagelkerkeR2(mod.full)$R2 - NagelkerkeR2(mod.baseline)$R2
prs.significance = lrtest(mod.baseline, mod.full)
```

```
r$> summary(mod.full)                                                                                                                          

Call:
glm(formula = Diagnosis ~ PRS0.500000 + Sex + Age + C1 + C2 + 
    C3 + C4 + C5, family = binomial, data = data_app)

Deviance Residuals: 
     Min        1Q    Median        3Q       Max  
-2.57234  -0.73259  -0.07879   0.78486   1.65693  

Coefficients:
              Estimate Std. Error z value Pr(>|z|)  
(Intercept)  3.542e+01  1.518e+01   2.334   0.0196 *
PRS0.500000  3.348e+04  1.388e+04   2.413   0.0158 *
SexM         2.585e+00  1.377e+00   1.877   0.0605 .
Age         -1.363e-01  7.423e-02  -1.836   0.0664 .
C1           1.642e+02  1.527e+02   1.075   0.2823  
C2           2.919e+02  2.543e+02   1.148   0.2510  
C3           1.559e+01  3.312e+02   0.047   0.9624  
C4           4.746e+02  2.440e+02   1.945   0.0518 .
C5           2.093e+02  1.988e+02   1.053   0.2923  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

    Null deviance: 74.192  on 53  degrees of freedom
Residual deviance: 46.660  on 45  degrees of freedom
AIC: 64.66

Number of Fisher Scoring iterations: 6


r$> summary(data_app$Diagnosis)                                                                                                                
Control    Case 
     30      24 
```

Yes, so higher PRS is associated with the higher level, or Case in this
variable.



# TODO
  * 