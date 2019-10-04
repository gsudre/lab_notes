# 2019-10-02 10:23:02

Quick note on creating correlation matrices from fMRI data GenR sent us:

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

# 2019-10-04 10:47:04

Then Pat said they got some more data, so here's how I parsed them:

```python
from scipy.io import loadmat
import numpy as np
import pandas as pd
import glob


rois = ['roi%03d' % (i+1) for i in range(160)]

cfiles = np.sort(glob.glob('/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/genR_to_dos160_*.mat'))

for cfname in cfiles:
    subj = cfname.split('_')[-1].replace('.mat', '')
    xc = loadmat(cfname)['tc_filt']
    cc = np.corrcoef(xc.T)
    df = pd.DataFrame(cc, index=rois, columns=rois)
    out_fname = '/Volumes/Shaw/GenR_peer_networks_brain_morphology/DOS_160/crosscorr/%s.csv' % subj
    df.to_csv(out_fname)
```


```bash
for line in `cat ~/tmp/cpt.csv`; do
    name=`echo $line | cut -d"," -f 2`;
    dt=`echo $line | cut -d"," -f 1`;
    grep $dt ~/tmp/RESPONSES_10032019.txt;
done
```