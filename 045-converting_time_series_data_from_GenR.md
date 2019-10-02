# 2019-10-02 10:23:02

Qhick note on creating correlation matrices from fMRI data GenR sent us:

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob


fid = open('/Users/sudregp/data/philip/GenR_FCTimeSeries/rois.txt', 'r')
rois = [line.rstrip() for line in fid]
fid.close()

cfiles = np.sort(glob.glob('/Users/sudregp/data/philip/GenR_FCTimeSeries/Cortical/genR_to_HOCortical_*.mat'))
scfiles = np.sort(glob.glob('/Users/sudregp/data/philip/GenR_FCTimeSeries/SubCortical/genR_to_HOSub_*.mat'))

for cfname, scfname in zip(cfiles, scfiles):
    subj = cfname.split('_')[-1].replace('.mat', '')
    xc = loadmat(cfname)['x']
    xsc = loadmat(scfname)['x']
    y = np.hstack([xc, xsc])
    cc = np.corrcoef(y.T)
    df = pd.DataFrame(cc, index=rois, columns=rois)
    out_fname = '/Users/sudregp/data/philip/GenR_FCTimeSeries/crosscorr/%s.csv' % subj
    df.to_csv(out_fname)
```