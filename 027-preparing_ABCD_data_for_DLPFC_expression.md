# 2019-07-02 11:54:06

We will assume the ABCD data released in 2.0 is properly cleaned, at least
according to what they described in their README, which is quite similar to what
we do in our pipeline.

So, the idea is just to prepare it to be imputed in the Michigan server.

```bash
# bw sin
module load plink
cd /lscratch/$SLURM_JOBID
for i in {1..22}; do
    plink --bfile ABCD_release_2.0_r2 --chr ${i} --recode-vcf --out vcf_chr${i};
done
module load vcftools
for f in `ls vcf*.vcf`; do
    echo $f;
    vcf-sort $f | bgzip -c > ${f}.gz &
done
```

And we should do a sneak peak using shapeit to see what kind of flips and
removals we'll need to do:

```bash
# bw sin
module load shapeit
refdir=/fdb/impute2/1000Genomes_Phase3_integrated_haplotypes_Oct2014/1000GP_Phase3/
for c in {1..22}; do
    echo $c >> chrs.txt;
done;
cat chrs.txt | parallel --max-args=1 -j $SLURM_CPUS_PER_TASK \
    shapeit -check -T 16 -V vcf_chr{}.vcf.gz \
        --input-ref $refdir/1000GP_Phase3_chr{}.hap.gz \
        $refdir/1000GP_Phase3_chr{}.legend.gz $refdir/1000GP_Phase3.sample \
        --output-log chr{}.alignments;
```

Now we need to properly flip and remove all bad ids. First, format the files:

```bash
# bw sin
for c in {1..22}; do
    grep Strand chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> flip_snps.txt
done
for c in {1..22}; do
    grep Missing chr${c}.alignments.snp.strand | cut -f 4 | sort | uniq >> missing_snps.txt;
done
plink --bfile ABCD_release_2.0_r2 --flip flip_snps.txt \
    --exclude missing_snps.txt --make-bed --out ABCD_release_2.0_r2_flipped
```

Let's just make sure the set is still clean:

plink --bfile merged_noDups_clean_flipped --missing --out ibd
awk '$6 > .01 {print $0}' ibd.imiss | wc -l
awk '$5 > .01 {print $0}' ibd.lmiss | wc -l

And construct the VCFs as above to send it to the imputation server.

for i in {1..22}; do   plink --bfile merged_noDups_clean_flipped --chr ${i} --recode-vcf --out vcf_chr${i}_flipped; done
for f in `ls vcf*flipped.vcf`; do echo $f; vcf-sort $f | bgzip -c > ${f}.gz; done
