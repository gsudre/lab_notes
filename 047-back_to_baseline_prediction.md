# 2019-10-11 13:54:04

After all the heritability work, I have a few new ideas for baseline prediction.
I think it makes sense to start with a clinical application, basically with only
people who have ADHD at baseline. We can do this regardless of whether they have
a second clinical assessment or not, as we'll use it just to get a normative
distribution of ADHD scans. 

Then, we'll do the two-step outlier detection QC, and finally get the earliest
good scan for each subject. We'll then keep only people that have one or more
clinical assessments after the earliest good scan, and compute the slope to the
next assessment, and also the last assessment. We can also binarize them into
improvers and non-improvers, and give a probability of getting better score. So,
in the end we'd have 4 different outcomes (2 continuous, 2 binary).

For sanity check, we should keep the QC variables in the dataset and check if
they provide a better accuracy than the actual data.

When scripts are all working, we could play with the OD threshold.

The trick here will be combining the datasets. How do we merge across domains,
if the slope is calculated with respect to the modality? I think in that
situation I'll need to re-run the QC using all modalities simultaneously,
fixating on dates that have all 3 modalities. Then we can simply merge the
closest neuropsych to it.

I'll start with DTI just to get the scripts ready...

```r
source('~/research_code/baseline_prediction/prep_dti_JHUtracts_data.R')
```