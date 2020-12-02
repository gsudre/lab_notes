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



```bash
module load R
Rscript /data/NCR_SBRB/software/PRSice_2.2.5/PRSice.R  \
    --prsice /data/NCR_SBRB/software/PRSice_2.2.5/PRSice_linux \
    --base ~/pgc2017/adhd_jun2017  \
    --target /data/NCR_SBRB/NCR_genetics/v2/1KG/NCR_1KG_genop05MAFbtp01rsbtp9_renamed \
    --all-score \
    --lower 5e-08 --upper .5 --interval 5e-05 \
    --no-regress \
    --out NCR_1KG_PRS_adhd_jun2017

Rscript /data/NCR_SBRB/software/PRSice_2.2.5/PRSice.R  \
    --prsice /data/NCR_SBRB/software/PRSice_2.2.5/PRSice_linux \
    --base ~/pgc2017/adhd_eur_jun2017  \
    --target /data/NCR_SBRB/NCR_genetics/v2/1KG/NCR_1KG_genop05MAFbtp01rsbtp9_renamed \
    --all-score --stat OR \
    --lower 5e-08 --upper .5 --interval 5e-05 \
    --no-regress \
    --out NCR_1KG_PRS_adhd_eur_jun2017

```

I'm having all kinds of issues here. I asked KWangmi for help... 

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

Maybe try beagle for phasing?
Maybe calculate PRS on the non-imputed data? what's the overlap?