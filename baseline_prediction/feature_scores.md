# 2018-09-27 12:44:16

Let's play a bit with plotting the DTI results back in the brain. I downloaded a
couple successful models form the cluster (AUC close to 1), in the n272 dataset.
So, to retrieve the feature scores, we do:

```r
library(h2o)
h2o.init()
m1 = h2o.loadModel('/Users/sudregp/tmp/DeepLearning_grid_0_AutoML_20180927_103833_model_50')
m2 = h2o.loadModel('/Users/sudregp/tmp/DeepLearning_grid_0_AutoML_20180927_103854_model_243')
a = read.table('~/data/baseline_prediction/dti_voxels_272/0968_ad.txt')
cnames = sapply(1:nrow(a), function(d) sprintf('v%05d', d))
rownames(a) = cnames
a[, 4] = 0
a[m1@model$variable_importances$variable, 4] = m1@model$variable_importances$relative_importance
write.table(a, file='~/tmp/m1.txt', row.names=F, col.names=F)
a[, 4] = 0
a[m2@model$variable_importances$variable, 4] = m2@model$variable_importances$relative_importance
write.table(a, file='~/tmp/m2.txt', row.names=F, col.names=F)
```

Then we need to put it into a .nii:

```bash
cd ~/tmp
cat m1.txt | 3dUndump -master ~/data/baseline_prediction/mean_272_fa_skeleton_mask.nii.gz -ijk -datum float -prefix m1 -overwrite -
cat m2.txt | 3dUndump -master ~/data/baseline_prediction/mean_272_fa_skeleton_mask.nii.gz -ijk -datum float -prefix m2 -overwrite -
```

So, the results are very sparse as one would imagine. After uni01, only about 135 out of
12K voxels are at play. Still, we get some 3-5 voxel clusters, in regions that
are similar for the nonew and without that restriction. They also make some
sense on where they are, even though I didn't put it into a template. They're
the usual parietal regions, as also some ventral/inferior candidates (uncinate)?

Now, does the n223 somewhat mirror the n272 result? 




