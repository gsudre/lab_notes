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
- [ ] re-run nuisance imaging variables just to make the points again using the
  current testing framework
  
  **2018-10-29 13:12:25**

  Had a quick chat with Philip today. We should use some network parameters for
  fMRI, see if they help. Also, it might be interesting to check variability
  within the scan session? Basically, split the session into two and check the
  absolute difference in correlation (could be over the whole brain or specific
  regions). We could even derive a difference distribution, and test the mean
  and standard deviation of that distribution, comparing them between NV and
  ADHDs.

  Something else Philip suggested was to look at the heritability of network
  parameters in fMRI (for a future project, using the family data), and also to
  look at https://www.stats.ox.ac.uk/~snijders/siena/ to look at the evolution
  of our networks. We could do this within a scan session, or across
  longitudinal progressions of a correlation matrix, or even using individual
  levels of a family tree as different time points in the analysis. 

  **2018-10-30 16:17:25**

  Alright, so our current situation is that there are some (but unimpressive) results from DTI, and same thing but even worse from rsFMRI. Nothing in structural ROIs, SNPs, or an other domain. So, here are some of the plates I have spinning:

 * Downsampling voxelwise structural
 * Same-space analysis (Calhoun) for all brain data 
 * ICASSO analysis on data matrices (only AD so far)
 * Combine raw results that are working to see if we get any boost
 * Network analysis on fMRI data to use instead of correlation matrix

I don't want to mess too much with ICASSO as I'm seeing some recent literature not as much in favor of it anymore. I do want to explore the downsampling issue, as well as the same space analysis a bit more. Downsampling will be a good compliment to the other results, and will go well with the combininng goal for later.

Of course, network analysis will be worth a try not only for this project, but also for the family heritability project later.

Increasing smoothing didn't help with downsampling the structural data... what if I use something from here? https://brainder.org/2016/05/31/downsampling-decimating-a-brain-surface/

**2018-10-31 09:25:38**

Just had a chat with Philip, and he suggested I should try FDR again and also a
spatial filter. I think they're good ideas, and I particularly like the spatial
filter option that I had previously discarded because it would be somewhat hard
to implement. Let's see how hard it would actually be.

**2018-11-01 09:42:03**

Philip sent me this paper (https://www.nature.com/articles/nn.4179), so I need
to check if they have an atlas of the network that I can use, or at least I
should try to reproduce their methods.

**2018-11-07 10:32:54**

For the spatial filter, what if we averaged within blobs? That would reduce the
number of variables, and maybe still perform well?

**2018-11-09 15:14:21**

Not sure what to do with rsFMRI yet. I don't have a good analogy for spatial
averages there. I also don't know if using network metrics won't be better.
There doesn't seem to be a consensus in the field. Maybe I just concatenate
everything as raw first? Well, struct won't work. Or, the other option is to do
ICA there as well. 

**2018-11-20 14:27:26**

Just had a catch up chat with Philip. I should keep plowing away with fMRI, see
if anything else comes out of it. But I should maybe focus on diagnostic, and
then use those variables for prognostic. Or use our own old results, as a way of
reducing variables. The story could be: these are good for diagnosis, but not
for prognosis. These other regions we have found before are also not good. And
even doing a major ML effort is also not good. But do visualize the diagnostic
results. Also, try the more simple stuff (not brain), like regressions and
groupwise tests, maybe even logistic regression, presenting odds-ratio tabls.

I should still explore individual symptoms and voting to reduce variance. Maybe
use just first and last symptom and define categories within the symptoms. Just
pick the one or two in each category (inatt, H, I) that has the most variance.

Philip is OK with the voting idea, and also fMRI variability, but be sure for
controlling for movement there.

I also thought: could I use GLM instead of DL in the spatial results, just to
reduce overfitting, variance and dependence in the test set size?

