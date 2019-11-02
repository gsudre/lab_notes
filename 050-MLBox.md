# 2019-11-01 16:14:59

Just because I'm not getting any good results with TPOT, let's give a try with
MLBox, before we go back to H2O.

Here's some code I ran that worked (had to use mlbox environment in BW):

```python
import mlbox
import os
import multiprocessing
home = os.path.expanduser('~')
phen_fname = home + '/data/baseline_prediction/dti_JHUtracts_ADRDonly_OD0.95.csv'
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
import sys
data = pd.read_csv(phen_fname)
target = 'SX_HI_groupStudy'
features_fname = home + '/data/baseline_prediction/ad_rd_vars.txt'
fid = open(features_fname, 'r')
feature_names = [line.rstrip() for line in fid]
fid.close()
data.rename(columns={target: 'class'}, inplace=True)
data['class'] = data['class'].map({'improvers': 1, 'nonimprovers': 0})
myseed=42
training_indices, validation_indices = train_test_split(data.index,
                                                        stratify = target_class,
                                                        train_size=0.8,
                                                        test_size=0.2,
                                                        random_state=myseed)
df = {'train': data.loc[training_indices, feature_names],
      'test': data.loc[validation_indices, feature_names],
      'target': data.loc[training_indices, 'class']}
df = mlbox.Drift_thresholder().fit_transform(df)
space = {
        'ne__numerical_strategy' : {"space" : [0, 'mean']},

                'ce__strategy' : {"space" : ["label_encoding", "random_projection", "entity_embedding"]},

                        'fs__strategy' : {"space" : ["variance", "rf_feature_importance"]},
                                'fs__threshold': {"search" : "choice", "space" : [0.1, 0.2, 0.3]},

                                        'est__strategy' : {"space" : ["LightGBM"]},
                                                'est__max_depth' : {"search" : "choice", "space" : [5,6]},
                                                        'est__subsample' : {"search" : "uniform", "space" : [0.6,0.9]}

                                                                }
opt = mlbox.Optimiser(scoring='roc_auc', n_folds=3)
best = opt.optimise(space, df, max_evals = 5)
model = mlbox.Predictor().fit_predict(best, df)
```

