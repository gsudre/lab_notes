# 2018-12-14 11:42:06

I'm working on getting the genotype data ready for submission. Following the
Evernote note and Philip's recommendations, I can start from
merged_noDups_clean, which has all boxes and it's been cleaned. But it has all
996 samples, and Philip said we should only send the 500 or so we used in the
analysis. So, now I need to start trimming them. I first figured out the longNSB
and family ID for the 512 good samples to share in Wendy's sheet. But if I'm
doing that, then I might as well just use the flipped version, which is ready
for imputation. Then:

```bash
cd /home/sudregp/data/prs/geno3
mkdir dbgap
cd dbgap
plink --bfile ../merged_noDups_clean_flipped --keep keep_long_nsbs.txt --make-bed --out dbGaP
```

I need to chat with Jen to see if I need to rename all of them, or if I can just
leave them as is. In any case, here are the next steps to prepare the data.

But I think we might eb able to just submit the PLINK files.

