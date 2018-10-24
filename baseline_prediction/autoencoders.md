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

I've been playing with the optimizing autoencoder code (non-convulaiton), but I'm reaching an asymptote of .0317, not overfitting. But that's using a very simple sigmoid model with 50 hidden neurons. I've trying much more complex stuff, deeper and bigger, but nothing is doing better. Maybe try playing with other parameters?

Also, it's true that using GPU makes it muuuuuuch faster.

Just note that we're not making use of the spatial information in the images, so that's why it might be better to use the convolution networks. We just need to understand how to best capture the lower dimensional model.

Use this example instead!!! https://blog.keras.io/building-autoencoders-in-keras.html

they show how to get the low-dimensional representation, and we can always adapt the training later.