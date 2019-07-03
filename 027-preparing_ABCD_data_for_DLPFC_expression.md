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

```bash
# easier to copy to the imputation server this way
cp vcf_chr*.vcf.gz ~/data/datashare/;
a='';
for m in {1..22}; do
    a=$a' 'https://hpc.nih.gov/~sudregp/vcf_chr${m}.vcf.gz;
done
echo $a
```

It turns out that I didn't need to flip the starnds on this one, and apparently
the imputation server processed them correctly. I'll store the results in
shaw/ABCD after decompressing, so we don't need to keep on storing the password.
