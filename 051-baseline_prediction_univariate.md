# 2019-11-04 14:45:06

Let's then explore the univariate filter idea for baseline prediction. If
anything, it'll better resemble usual neuroimaging papers.

Teh idea is that we'll split the data into training, validation, and testing. We
will pick univariately inside training, and do multiple combinations of
linearSVC, SVC, or Ensembles with the top X features.

