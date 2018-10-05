# 2018-10-05 13:06:52

I've been struggling with what would be the best way to show results under a
random scenario. For classification, we could run the exact same algorithm after
data randomization, then show results using simply the majority class, and then
finally using the class probabilities for voting.

In the prediction scenario, we could use random data (as above), but also the
mean and the median.

Then, the idea is that we'd have a series of plots, and in each plot would have
a series of boxplots (one per phenotype), and horizontal bars showing the 95th
percentile of each of the random approaches above. There would be a plot for
accuracy, one for AUC, F1-ratio, specificity, and sensitivity. Or maybe get rid
of F1 just to make it symmetric, or add precision to the plots (recall = sensitivity).

For prediction, we can play with RMS and R2.

---

I then realized that we can simply use the probabilities based on the training
set, and not have to do majority and then probabilities. That takes care of
classification (just need to run the random stuff).

So, the code in automl/random_autoValidation.R is working. We just need to run
it a whole bunch of times to get some distributions. Note, a few things:

* Random AUC should be aroiund .5 anyways, so there shouldn't be much variation there.
* Most of the variation will be in the classification accuracy.