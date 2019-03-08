# 2019-03-07 09:41:37

Some quick analysis ot check if the change in tract is heritable within
families. We'll start with tract, but if there is a sniff of result there we can
quickly move to anatomy and resting state.

Because it's a quick analysis, I'll start with our usual TORTOISE pipeline. We
can always re-run this later with the new FDT pipeline. So, I started with a
Labmatrix query that exclude everyone that had errors importing and bad data. 
I'm also dealing only with scans that have already gone through
DTI-TK, because we cannot wait for more preprocessing.

So, the result is not_failed_tortoise_DTI.xlsx, which shows 283 individuals from
196 nuclear families (25 extended within that), in a total of 694 scans.

But we need to first cut out any families made up of single individuals, and
individual with only one scan

```r
m = read.csv('~/data/heritability_change/dti_tracts_qc.csv')
m = m[!is.na(m$norm.trans), ]
keep_me = c()
for (s in unique(m$Medical.Record...MRN)) {
    if (sum(m$Medical.Record...MRN == s) > 1) {
        keep_me = c(keep_me, which(m$Medical.Record...MRN == s))
    }
}
m = m[keep_me, ]
```

I noticed that one of the mask IDs had NAs for movement QC. Likely the RPd files
weren't transferred... I won't worry about it as this will unlikely be our final
pipeline, and just remove the mask ID.

So, we go from 693 scans to 647. Next, cut based on families:

```r
keep_me = c()
for (s in unique(m$Nuclear.ID)) {
    fam_mrns = m[m$Nuclear.ID == s, "Medical.Record...MRN"]
    if (length(unique(fam_mrns)) > 1) {
        keep_me = c(keep_me, which(m$Nuclear.ID == s))
    }
}
for (s in unique(m$Extended.ID)) {
    fam_mrns = m[m$Extended.ID == s, "Medical.Record...MRN"]
    if (length(unique(fam_mrns)) > 1) {
        keep_me = c(keep_me, which(m$Extended.ID == s))
    }
}
m = m[unique(keep_me), ]
```

OK, now we're down to 316 scans, for 125 individuals, in 68 nuclear families (21
extended).

Within that, we can finally pick the two first good time points to avoid any
sort of outcome analysis.

Let's take a look at the distributions first:

![](images/2019-03-07-17-40-41.png)

They look quite normal, so our mu+-3Sd thresholds should be fine:

```r
for (p in c('mean_fa', 'mean_ad', 'mean_rd')) {
    mu=mean(m[, p])
    s = sd(m[, p])
    print(p)
    print(mu - 3 * s)
    print(mu + 3 * s)
}
```

```
[1] "mean_fa"
[1] 0.3326389
[1] 0.4107519
[1] "mean_ad"
[1] 1.096487
[1] 1.292673
[1] "mean_rd"
[1] 0.5669469
[1] 0.7447344
```

I don't want to normalize the other metrics because it doesn't make much sense
to have both tails discarded. We can also do it by eye. So, our final thresholds
were:

```r
m = m[!is.na(m$norm.trans), ]
idx = (m$mean_fa < 0.4107519 & m$mean_fa > 0.3326389 &
       m$mean_ad < 1.292673 & m$mean_ad > 1.096487 &
       m$mean_rd < 0.7447344 & m$mean_rd > 0.5669469) # &
    #    m$norm.trans < 2.5 & m$norm.rot < 0.04 & m$missingVolumes < 4)
m = m[idx,]
```

That leaves us with 309 out of the 316 scans. Even worse if I include the other
QC thresholds, so maybe use them as covariates later. Now we just need to select the
first two for each person:

```r
keep_me = c()
for (s in unique(m$Medical.Record...MRN)) {
    subj_scans = m[m$Medical.Record...MRN==s, ]
    dates = as.Date(as.character(subj_scans$"record.date.collected...Scan"),
                                 format="%m/%d/%Y")
    if (length(dates) >= 2) {
        sdates = sort(dates)  # not sure why index.return is not working...
        # make sure there is at least 6 months between scans
        next_scan = 2
        while (((sdates[next_scan] - sdates[1]) < 180) && (next_scan < length(sdates))) {
            next_scan = next_scan + 1
        }
        dob = as.Date(as.character(m$Date.of.Birth), format="%m/%d/%Y")
        first_scan_age = (sdates[1] - dob)/365.25
        if (((sdates[next_scan] - sdates[1]) >= 180) && (first_scan_age < 26)) {
            idx1 = which(dates == sdates[1])
            keep_me = c(keep_me, which(m$Mask.ID == subj_scans[idx1, 'Mask.ID']))
            idx2 = which(dates == sdates[next_scan])
            keep_me = c(keep_me, which(m$Mask.ID == subj_scans[idx2, 'Mask.ID']))
        }
    }
}
m2 = m[keep_me, ]
```

So, with clean scans we have 122 people (244 scans). Then, we also need to
remove people who are by themselves, after cleaning.

```r
> table(m2$Extended.ID)
9000 9002 9003 9005 9007 9008 9010 9011 9014 9016 9019 9022 9025 9028
   2    2    4    4    2    4    4    2    2    2    2    2    2    6
> sum(table(m2$Nuclear.ID)>2)
[1] 17
```

So, in the end we need either extended ID >= 4 (2 people * 2 scans), or the same
for nuclear ID.

```r
good_nuclear = names(table(m2$Nuclear.ID))[table(m2$Nuclear.ID) >= 4]
good_extended = names(table(m2$Extended.ID))[table(m2$Extended.ID) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, which(m2$Nuclear.ID == f))
}
for (f in good_extended) {
    keep_me = c(keep_me, which(m2$Extended.ID == f))
}
m3 = m2[unique(keep_me), ]
```

Oh wow, now we're down to only 224 scans, on 122 people (67 nuclear families)!
Not sure we can go through with it... but let's try it anyways.

# 2019-03-08 09:48:34

First step is to calculate the tracts for each of the scans.

```bash
sed "s/$/_tensor_diffeo\.nii/g" ids224.txt > tensors224.txt
cd /Volumes/Shaw/dti_robust_tsa/analysis_may2017/
export DTITK_ROOT=/Applications/dtitk-2.3.3-Darwin-x86_64/
/Applications/dtitk-2.3.3-Darwin-x86_64/scripts/tsa_sampling ~/data/heritability_change/ids224.txt ../ixi_aging_template_v3.0/tsa/ mean
# edit them first
python ~/research_code/lab_mgmt/convert_dti_sampling.py
Rscript ~/research_code/dti/compile_tract_table.R
```

Then, we residualize the brain before computing the slope on that:

```r
a = read.csv('~/data/heritability_change/dti_tracts_qc_clean.csv')
b = read.csv('~/data/heritability_change/dti_mean_phenotype_224.csv')
m = merge(a, b, by.x='Mask.ID', by.y='file')
tract_names = colnames(b)[2:ncol(b)]
library(MASS)

res.lm <- lm(inatt_vol_lh.x ~ Sex...Subjects + int_avg_freesurfer5.3 + age_at_scan + I(age_at_scan^2) + ext_avg_freesurfer5.3 + int_avg_freesurfer5.3 + mprage_QC, data = df_clean)
step <- stepAIC(res.lm, direction = "both", trace = F)
print(summary(step))
mycluster = residuals(step)
```

Note that I don't want to use the familial relationships at this point yet (i.e.
not use lme for residuals)... not sure if it'd bias the residuals in sending to
SOLAR later. Also, even though we could potentially say that residualizing based
on age is not necessary, as we'll take the slope later, we do need to do it for
QC. Maybe we could even try both?