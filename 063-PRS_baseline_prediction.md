# 2020-01-06 09:55:15

Let's see if we can use the PRS for baseline prediction. We'll define the
possible group of people as only who has PRS, and go from there. We'll pick the
best PRS target (as is significant, but not necessarily clinically meaningful),
and then will tyr to improve it by throwing in new phenotypes.

If PRS fails from the get-go, we'll try using the phenotypes shown in prrevious
papers to be good baseline predictors. Those phenotypes are the easy ones to get
(e.g. initial symptom count, co-morbidities), and then, using the same approach,
we'll check how much our deep phenotyping adds to it.

A few things to keep in mind:

* Use DICA_on entries judiciously. Check the symptom history, whether the
  medication is making a difference, etc. 
* Philip does not want to use longitudinal data here.
* We'll need to do some data descriptives as well. 
* Need to work in the PCs and possibly do WNH-only analysis when playing with PRS.
* Might need to also add in the easy initial phenotypes to compare to PRS. The
  downside there is that those phenotypes were studied in the context off
  outcome, and our population might not be old enough to have those phenotypes
  as the best predictors to begin with.
* Might need to play with different thresholds of initial symptoms to see who to
  keep (e.g. include NVs?)

# Prepping the data

The first step is to get the symptom trajectories and subset them to the people
that have PRS. Then, we can clean up the DICA_on entries.

```r
clin = read.csv('~/data/baseline_prediction/prs_start/clinical_01062020.csv')
prs = read.csv('/Volumes/NCR/reference/merged_NCR_1KG_PRS_12192019.csv')
keep_me = clin$MRN %in% prs$MRN
clin2 = clin[keep_me, ]
long = names(table(clin2$MRN))[table(clin2$MRN)>1]
keep_me = clin2$MRN %in% long
clin3 = clin2[keep_me, ]
clin3$entry = sapply(1:nrow(clin3),
                     function(x) sprintf('%s_%s', clin3$MRN[x], clin3$DOA[x]))
clin4 = clin3[!duplicated(clin3$entry), ]
long = names(table(clin4$MRN))[table(clin4$MRN)>1]
keep_me = clin4$MRN %in% long
clin5 = clin4[keep_me, ]

write.csv(clin5,
          file='~/data/baseline_prediction/prs_start/long_clin_01062020.csv',
          row.names=F)
```

I then removed any DICA_on entries that seem to have affected either SX count,
and kept only people with baseline <= 16. Now, we need to decide if we include
evryone or just people with some SX at baseline. 
