# 2018-12-06 10:20:19

Given the latest results with the descriptives, I got curious about whether
using a robust lme method wouldn't be better, as it would protect against
outliers. Not that all our results have outliers, but some do and it would be
better to avoid those, if still keeping the same results that are good. 

First, let's check whether that solves our issue.

I played a bit with robustlmm and rlme. The latter didn't quite work. And
robustlme did work, but getting p-values for the fixed estimates is non-trivial.
We could just approximate the t-stats and get a p-value from that. But I think
that at this time it might be just easier to ignore results that seem driven by
outliers. 

I tried running MDS, but the maskids that came out as outliers weren't the same
as the ones in the structural result. The other option is to filter outliers as
NAs, and use the regular LME, as it omits NAs...

Let's try that.

It didn't quite work, as many of the p-values were being stored as NAs,
regardless of how I chose to deal with NAs. Even if I just send the good data,
maybe it doesn't have enough data to compute the p-values?

Another option is to also look at the subjScale pipeline, because there's a
chance that they wouldn't be outliers after that transform, and those results
mostly tracked the no transformation results. Just a bit smaller.

So, two questions:

1) Are subjScale results around the same region as None?
2) Are subjScale less succeptible to outliers, when compared to None?



