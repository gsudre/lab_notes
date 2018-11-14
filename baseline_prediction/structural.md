# 2018-09-12 14:39:51

We're having some issues running even the limited version of the autoML scripts.
Not sure what's going on, but it doesn't even look like we're getting to run
automl, as we're not stopping after some time. Let's run a few tests...

# 2018-09-13 11:52:26

It looks like it takes a long time to just prepare the data. But using all the
variables in the structural (voxels) we already go up to 115Gb, and then the
models start running. I'm sure it will break as soon as they run longer... ye,
it did.

That's fine, as our best models were all using univariate, so we could go with
that, or even PCA. We just need to make sure the code is working fine for that.

# 2018-11-14 15:32:22

I managed to downsample the structural data. So, let's recreate the data files and then run the within domain tests again:

```r
library(gdata)
clin = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'mprage')
clin = clin[, c(1,4)]
colnames(clin) = c('MRN', 'mask.id')
for (m in c('area', 'thickness', 'volume')) {
    print(m)
    load(sprintf('lh.%s.ico4.gzip', m))
    cnames = sapply(1:ncol(data), function(x) sprintf('v_lh_%04d', x))
    colnames(data) = cnames
    lh=data
    load(sprintf('rh.%s.ico4.gzip', m))
    cnames = sapply(1:ncol(data), function(x) sprintf('v_rh_%04d', x))
    colnames(data) = cnames
    data = cbind(lh, data)
    mask.id = read.table('subjects_list.txt')[,1]
    data = cbind(mask.id, data)
    data = merge(clin, data, by='mask.id')
    save(data, file=sprintf('~/data/baseline_prediction/struct_%s_11142018_260timeDiff12mo.RData.gz', m), compress=T)
}
```
