# 2020-02-19 15:51:03

After going over summary of results with Philip, a few more things to try:

* covary age and sex in everything
* try WNH only for PRS univariate. Is it better?
* check the difference between clinicals, make sure the minimal difference is reasonable
* add sex as a possible univariate, but then try excluding and including it for FDR
* do univariate selection for 2-class, 3-class and 4 class individually
* remove externalizing, internalizing and med rates for 4 groups univariate,
  possibly even 3 groups FDR
* try adding base_hi and base_inatt together
* show proportion of variance for the 2 group models with base_sx included
* code medication to have 2 groups only, combining none and nostim

Alright, let's get going then.