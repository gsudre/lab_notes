# 2018-10-24 15:23:18

The idea here is to optimize as much as we can the autoencoders to find an optimal lower respresentation space of our data. We can easily think of our fMRI correlation matrix as a 2D image, and the DTI as a 3D image. I'm not sure what to make of the structural image, but if it's not in a 3D format we can easily shape it to be that (or even 2D).

So, I'm following code from:

https://www.datacamp.com/community/tutorials/autoencoder-keras-tutorial

but also 

https://blog.keras.io/building-autoencoders-in-keras.html

So, first thing is to save the fMRI correlation matrices as 2D, or the DTI as 3D. Let's go with DTI because we're benchmarking much of our other stuff like that.

```python
import nibabel as nb
import pandas as pd
import numpy as np

fname = '/Users/sudregp/data/baseline_prediction/dti_gf_09212018_223timeDiff12mo.csv'
pheno = pd.read_csv(fname)
all_data = []
for m in pheno['mask.id']:
    print(m)
    data_fname = '/Volumes/Shaw/dti_robust_tsa/analysis_may2017/%04d_tensor_diffeo_ad.nii.gz' % m
    img = nb.load(data_fname)
    subj_data = img.get_data()
    all_data.append(subj_data)
d = np.array(all_data)
np.savez('/Users/sudregp/tmp/ad_223_3d.npz', data=d)
```