# 2018-12-14 13:02:58

Let's figure out who we can use for testing (test_candidates.xslx). Within modality, I'll start by
looking at everyone who has a scan younger than 18, then remove the ones who are
already in the training set. That should leave out a group of people who we have
either already seen since I froze the data, or that Wendy could potentially call
back.

The way I'm shaping this analysis, we're only interested in ADHDNOS or ADHD, but
let's grab everyone and then we can decide.

For the people already in the study, instead of pulling from the clinical file
I'll pull from the scan files (in_train), because it's possible that a kid is in the
clinical file but not in the scans because he didn't scan well before freezing,
but did scan well afterwards. So, we want to include him this time.

The next step is to remove any scans considered bad under the current QC
thresholds (bad_now).

## structural
As it currently stands, we have 163 structural scans to go in, plus another 38
to be QCed (anything <=2.5). Have to check how many of those are ADHDs, but it's not a bad
number. Let's trim it a bit though. Many of them actually are duplicates. So,
I'll remove scans so that we only keep the first good baseline (keep_scan).
Then, as of now, we have 156 subjects for testing and 24 that could be added
after QC. The question now is how many are ADHD, how many still need to be
called, etc.

## DTI

The values I used for mean mask in the past are not current anymore, so I'll
need to re-run the DTI QC to see what's the actual min and maximum used in each
metric. But I don't want to change the dataset I'm playing with, so I just
regenerated all QC values, and used the min and max I had used before to trim
the new test data. Then, I ran the same procedure as in structural, to keep the
earliest good scan. Now, we have 178 for testing and 43 with NAs that could potentially go in.

These number seem high, but after I start splitting them into NV vs ADHD, and
checking which ones would still need to be called, they'll go down very fast.

## resting

I started looking at the test resting scans, but then I realized it would be
best if I gathered all the NAs first (For all modalities) before deciding which
scan to use. The reason there is that I was choosing the best scan if it was
later than an NA, but if the NA is actually good, then I'd rather choose the
earlier NA scan.







