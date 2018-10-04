# TODO

**2018-09-06 13:00:00**

These are from last chat with Philip:

- [x] Do the within ADHD analysis, using remission and persistent labels. But also in a separate analysis use ADHD_NOS (anyone with 3, 4, or 5 SX), and group them into improvers and nonimprovers (persistent and remission doesn't apply to NOS). Of course, using slopes would work for both analysis.
- [x] Try redoing the analysis by ignoring the people that get worse. They could be a group on their own. Should we even include them?
- [ ] Try re-balancing the classes in the CV to see if the results are any different.
- [ ] Plot results in the brain just to see if they make sense. 
- [x] Push the univariate analysis even further using FDR.
- [x] Use SNPs from PGC when selecting within our sample, instead of univariate, otherwise it would be a mini-GWAS that people won't like.

**2018-09-13 13:15:21**

From WIP idle mind:

- [x] GO ahead with cognitive data that Jen prepaered
- [ ] Try to get results with a test group in case Wendy's calling project
  doesn't yield too many people 
- [ ] Do basic univariate analysis for comparison
- [x] Start using rsfMRI data Jen and Esteban finished processing

**2018-09-18 15:41:14**

From quick chat with Philip today:

- [ ] Run everything on ever had ADHD vs NV, because genes should show something
  these, and also in SX at assessment
- [ ] Run genome-wide significance for SNPs, just to show it, even though nothing
  will come out of it
- [x] Run dummy classifiers to assess significance (i.e. not show only the best ones)
- [x] Use new DTI motion parameters to assess if we're using correct IDs.
  (will need to remove scans to reduce correlation! check R notebook note)
- [x] Try best uni models getting features from training set! (not necessarily
  bad if using replication set from Wendy, but definitely needed otherwise)

**2018-09-25 13:44:15**

- [ ] Maybe give Wendy more names of people to interview back, if modalities
  other than DTI perform well?
- [x] Maybe add sex, age, and movement variables to the models to see if they
  share some of the loadings
- [ ] Combine datasets. Here, we can try to use only people that have all
  variables, or go for allowing NAs. Most algorithms in autoML can handle them
  (other than GNB, I think). The main (empirical) question is to whether we will
  join all datasets, or only combinations of them. I imagine we could test every
  possible combination (using and not using NAs), but it makes more sense to
  just pick the best pipeline in each modality (e.g. only one of the rsFMRI
  aparcs, one one DTI modality, etc). In the worst case scenario, we could just
  clump together the good modalities, compared to everything, and show that it's
  not just a matter of having more features.
- [x] Run subgroups in latent groups as well (see main notes)
- [x] Change it to only remove NA columns for PCA (and maybe univariate) pipeline

**2018-09-27 10:25:38**

Another idle mind during WIP:

- [x] Try using the actual clinical data as future predictor
- [x] Check if good variables are at all correlated to sex, motion, or age, or
  if we could use just those variables to predict the outcome.
 

**2018-09-28 16:26:21**

Suggestion from Philip:

- [ ] Try decoding just using age, handedness, IQ, and sex, to simulate the best
  results in the ADHD-200 competition 
- [ ] Are these kids already remitted/persistent by the time we're using as baseline?
- [ ] Should we calculate slopes using as baseline the day we use as baseline in the phenotype? (That will be an issue when combining the data...).

