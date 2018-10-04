# 2018-10-04 14:15:21

Per PHilip's suggestion, I'm running the set of variables that got the best
ADHD200 results: IQ, age, handedness and sex. For that, I'm keeping NAs, because
we didn't have handedness or IQ for several folks.

Then, it's just a matter of properly setting up the RData file and the raw script
to take in the factors.

```r
data = read.csv('~/data/baseline_prediction/adhd200_10042018.csv')
colnames(data) = c('MRN', 'vCateg_sex', 'v_Age', 'vCateg_handedness', 'v_IQ')
save(data, file='~/data/baseline_prediction/adhd200_10042018.RData.gz', compress=T)
```

