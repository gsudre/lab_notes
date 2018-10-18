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

- [x] Try decoding just using age, handedness, IQ, and sex, to simulate the best
  results in the ADHD-200 competition 
- [ ] Are these kids already remitted/persistent by the time we're using as baseline?
- [ ] Should we calculate slopes using as baseline the day we use as baseline in the phenotype? (That will be an issue when combining the data...).

**2018-10-05 15:29:52**

Just brainstorming...

- [ ] What if we boosted our numbers by using our findings that the remitted
  brain looks more like the NV brain? So, we'd train using NV=remission to
  collapse it to only 2 classes, but only test on kids who actually have ADHD.
  The idea there is that, if the kid's brain (and SX?) will eventually look more
  like NV, then we can force them at this earlier time point to look the same as
  well. It would also boost our numbers, where remitted is usually the smallest.
  
  
  **2018-10-10 14:25:35**

  Just had a chat with Philip on next directions. He wants me to stick with a
  single model, but it's OK to find model parameters within the validation set.
  If I want, I can later test something else, or even a combination. For
  example, stick to DeepLearning, then try just XGB (or something else) for
  comparison. Or even both together. But for main results, he thinks it'll be
  easier to for people to understand a single model approach.

  He also thinks I shoulkd check right away how many of the people we're
  classifying have already remitted.

  Making ROC plots is a good idea because people are familiar with it. But
  the barplots I'm planning to show are also OK. He said there are guidelines on
  how to show these results (PRISM? PRISMA?). 

  Finally, we should focus on pairwise comparisons, mostly to include the NVs.
  We can do it for ADHDonly and ADHD_NOS, where the later the groups would be
  improvers and not improvers, based on slopes. We could also do that approach
  for ADHD_only. But it looks like slope prediction is a burst for now.
  
**2018-10-18 15:02:47**

Some more ideas after playing a bit with the function to plot the results:

- [ ] make ROC curves
- [ ] generate boxes for all other datasets. the main plot should only be the
  main phenotype for each data domain. But we should probably generate again the
  results within domain just to get a better idea.
- [ ] check if we should be doing worse in nvVSrem type of comparisons (maybe
  not if we can do nvVSadhd at baseline, baselically showing that the brain is
  not remitted yet at the time of scan.
- [ ] check that the brain regions are similiar for similar decoding tasks, but
  different for differen ones (e.g. nvVSper should be similar to nvVSnonimp, but
  different than perVSrem. Also, nvVSper and nvVSrem should be similar to
  nvVSadhd, when we run the latter.