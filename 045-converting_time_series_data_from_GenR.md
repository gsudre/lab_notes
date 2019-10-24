# 2019-10-02 10:23:02

Quick note on creating correlation matrices from fMRI data GenR sent us:

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob
from sklearn import preprocessing


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
    ys = preprocessing.scale(y)
    cc = np.corrcoef(ys.T)
    df = pd.DataFrame(cc, index=rois, columns=rois)
    out_fname = '/Users/sudregp/data/philip/GenR_FCTimeSeries/crosscorr/%s.csv' % subj
    df.to_csv(out_fname)
```

# 2019-10-04 10:47:04

Then Pat said they got some more data, so here's how I parsed them:

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob
from sklearn import preprocessing


rois = ['roi%03d' % (i+1) for i in range(160)]

cfiles = np.sort(glob.glob('/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/genR_to_dos160_*.mat'))

for cfname in cfiles:
    subj = cfname.split('_')[-1].replace('.mat', '')
    xc = loadmat(cfname)['tc_filt']
    ys = preprocessing.scale(xc)
    cc = np.corrcoef(ys.T)
    df = pd.DataFrame(cc, index=rois, columns=rois)
    out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/crosscorr/%s.csv' % subj
    df.to_csv(out_fname)
```

# 2019-10-17 11:17:59

Philip received Global Signal files from Andrew to generate GSR-based
correlation matrices. So, following the math in
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2694109/, I do this for each
dataset. He also asked me to add the amygdala from DK to the DS data atlas:

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob
import os


rois = ['roi%03d' % (i+1) for i in range(162)]

cfiles = np.sort(glob.glob('/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/genR_to_dos160_*.mat'))

for cfname in cfiles:
    subj = cfname.split('_')[-1].replace('.mat', '')
    gs_fname = '/Users/sudregp/Downloads/GlobalSignal/genR_GS_%s.mat' % subj
    if os.path.exists(gs_fname):
        xc = loadmat(cfname)['tc_filt']
        extra_fname = '/Users/sudregp/data/philip/GenR_FCTimeSeries/SubCortical/genR_to_HOSub_%s.mat' % subj
        ex = loadmat(extra_fname)['x']
        # left and right amygdala
        B = np.hstack([xc, ex[:, [12, 13]]])
        cc = np.corrcoef(B.T)
        df = pd.DataFrame(cc, index=rois, columns=rois)
        out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/crosscorr_plusAmygdala/%s.csv' % subj
        df.to_csv(out_fname)

        g = loadmat(gs_fname)['gs']
        gp = np.linalg.pinv(g)
        beta_g = np.dot(gp, B)
        Bp = B - np.dot(g, beta_g)

        cc = np.corrcoef(Bp.T)
        df = pd.DataFrame(cc, index=rois, columns=rois)
        out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/crosscorr_plusAmygdala_GSR/%s.csv' % subj
        df.to_csv(out_fname)
    else:
        print('Could not find GS file for %s' % subj)
```

Then, do the same for the other set of ROIs.

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob
import os


fid = open('/Users/sudregp/data/philip/GenR_FCTimeSeries/rois.txt', 'r')
rois = [line.rstrip() for line in fid]
fid.close()

cfiles = np.sort(glob.glob('/Users/sudregp/data/philip/GenR_FCTimeSeries/Cortical/genR_to_HOCortical_*.mat'))
scfiles = np.sort(glob.glob('/Users/sudregp/data/philip/GenR_FCTimeSeries/SubCortical/genR_to_HOSub_*.mat'))

for cfname, scfname in zip(cfiles, scfiles):
    subj = cfname.split('_')[-1].replace('.mat', '')
    gs_fname = '/Users/sudregp/Downloads/GlobalSignal/genR_GS_%s.mat' % subj
    if os.path.exists(gs_fname):
        xc = loadmat(cfname)['x']
        xsc = loadmat(scfname)['x']
        B = np.hstack([xc, xsc])
        g = loadmat(gs_fname)['gs']

        gp = np.linalg.pinv(g)
        beta_g = np.dot(gp, B)
        Bp = B - np.dot(g, beta_g)

        cc = np.corrcoef(Bp.T)
        df = pd.DataFrame(cc, index=rois, columns=rois)
        out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/rsfMRI_desikan_killany/crosscorr_gsr/%s.csv' % subj
        df.to_csv(out_fname)
    else:
        print('Could not find GS file for %s' % subj)
```

Philip also asked me to implement a band-pass filter. Just like the GSR
analysis, this is highly unorthodox. Both are normally done in voxel space, not
in ROI space. Still, it's theoretically possible. So, let's do it. According to
the GenR papers, the TR is 2s
(https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24064). So, let's create
the filter:

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob
import os


fid = open('/Users/sudregp/data/philip/GenR_FCTimeSeries/rois.txt', 'r')
rois = [line.rstrip() for line in fid]
fid.close()

cfiles = np.sort(glob.glob('/Users/sudregp/data/philip/GenR_FCTimeSeries/Cortical/genR_to_HOCortical_*.mat'))
scfiles = np.sort(glob.glob('/Users/sudregp/data/philip/GenR_FCTimeSeries/SubCortical/genR_to_HOSub_*.mat'))

for cfname, scfname in zip(cfiles, scfiles):
    subj = cfname.split('_')[-1].replace('.mat', '')
    gs_fname = '/Users/sudregp/Downloads/GlobalSignal/genR_GS_%s.mat' % subj
    if os.path.exists(gs_fname):
        xc = loadmat(cfname)['x']
        xsc = loadmat(scfname)['x']
        B = np.hstack([xc, xsc])
        g = loadmat(gs_fname)['gs']

        gf = mne.filter.filter_data(g.T, .5, .01, .1).T
        Bf = mne.filter.filter_data(B.T, .5, .01, .1).T

        gp = np.linalg.pinv(gf)
        beta_g = np.dot(gp, Bf)
        Bp = Bf - np.dot(gf, beta_g)

        cc = np.corrcoef(Bp.T)
        df = pd.DataFrame(cc, index=rois, columns=rois)
        out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/rsfMRI_desikan_killany/crosscorr_BP_gsr/%s.csv' % subj
        df.to_csv(out_fname)
    else:
        print('Could not find GS file for %s' % subj)
```

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob
import os
import mne


rois = ['roi%03d' % (i+1) for i in range(162)]

cfiles = np.sort(glob.glob('/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/genR_to_dos160_*.mat'))

for cfname in cfiles:
    subj = cfname.split('_')[-1].replace('.mat', '')
    gs_fname = '/Users/sudregp/Downloads/GlobalSignal/genR_GS_%s.mat' % subj
    if os.path.exists(gs_fname):
        xc = loadmat(cfname)['tc_filt']
        extra_fname = '/Users/sudregp/data/philip/GenR_FCTimeSeries/SubCortical/genR_to_HOSub_%s.mat' % subj
        ex = loadmat(extra_fname)['x']
        # left and right amygdala
        B = np.hstack([xc, ex[:, [12, 13]]])

        Bf = mne.filter.filter_data(B.T, .5, .01, .1).T

        cc = np.corrcoef(Bf.T)
        df = pd.DataFrame(cc, index=rois, columns=rois)
        out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/crosscorr_plusAmygdala_BP/%s.csv' % subj
        df.to_csv(out_fname)

        g = loadmat(gs_fname)['gs']
        gf = mne.filter.filter_data(g.T, .5, .01, .1).T

        gp = np.linalg.pinv(gf)
        beta_g = np.dot(gp, Bf)
        Bp = Bf - np.dot(gf, beta_g)

        cc = np.corrcoef(Bp.T)
        df = pd.DataFrame(cc, index=rois, columns=rois)
        out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/crosscorr_plusAmygdala_BP_GSR/%s.csv' % subj
        df.to_csv(out_fname)
    else:
        print('Could not find GS file for %s' % subj)
```
