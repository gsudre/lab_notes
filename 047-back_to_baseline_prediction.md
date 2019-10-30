# 2019-10-11 13:54:04

After all the heritability work, I have a few new ideas for baseline prediction.
I think it makes sense to start with a clinical application, basically with only
people who have ADHD at baseline. We can do this regardless of whether they have
a second clinical assessment or not, as we'll use it just to get a normative
distribution of ADHD scans. 

Then, we'll do the two-step outlier detection QC, and finally get the earliest
good scan for each subject. We'll then keep only people that have one or more
clinical assessments after the earliest good scan, and compute the slope to the
next assessment, and also the last assessment. We can also binarize them into
improvers and non-improvers, and give a probability of getting better score. So,
in the end we'd have 4 different outcomes (2 continuous, 2 binary).

For sanity check, we should keep the QC variables in the dataset and check if
they provide a better accuracy than the actual data.

When scripts are all working, we could play with the OD threshold.

The trick here will be combining the datasets. How do we merge across domains,
if the slope is calculated with respect to the modality? I think in that
situation I'll need to re-run the QC using all modalities simultaneously,
fixating on dates that have all 3 modalities. Then we can simply merge the
closest neuropsych to it.

I'll start with DTI just to get the scripts ready...

```r
source('~/research_code/baseline_prediction/prep_dti_JHUtracts_data.R')
```

# 2019-10-29 12:33:51

I should try the classification with and without adding age and sex to the
features, and then just checking how good our age and sex predictions can get by
themselves. For now, I'll play a bit with TPOT just to see how good this can
get.

Of course, by following the examples on the website, we should have multiple
validation sets, so that we can have a distribution. But maybe we can just set
the seeds, and go from there? Like, a list of seeds (starting with 42)?

For now, let's make sure a minimal example using TPOT is working. Then we can
start tweaking it a bit and running different seeds.

```python
from tpot import TPOTClassifier, config
from sklearn.model_selection import train_test_split
import pandas as pd
import numpy as np
from dask.distributed import Client
client = Client(n_workers=32, threads_per_worker=1)
myseed=1234

data = pd.read_csv('~/data/baseline_prediction/dti_JHUtracts_ADRDonly_OD0.95.csv')
data.rename(columns={'SX_HI_groupStudy': 'class'}, inplace=True)
data['class'] = data['class'].map({'improvers': 1, 'nonimprovers': -1})
print(data['class'].value_counts())

feature_names = [fname for fname in data.columns if fname[:2] in ['ad', 'rd']]
target_class = data['class'].values
training_indices, validation_indices = train_test_split(data.index,
														stratify = target_class, train_size=0.75,
														test_size=0.25,
														random_state=myseed,
														scoring='roc_auc')

# removing some warnings by hard coding parameters in the dictionary
my_config = config.classifier_config_dict
my_config['sklearn.linear_model.LogisticRegression']['solver'] = 'lbfgs'
preproc = [v for v in my_config.keys() if v.find('preprocessing') > 0]
for p in preproc:
	my_config[p]['validate'] = [False]

tpot = TPOTClassifier(n_jobs=-1, random_state=myseed, verbosity=2,
						config_dict=my_config, use_dask=True)

X = data[feature_names].values
tpot.fit(X[training_indices], target_class[training_indices])

### after
tpot.export('tpot_adrd_pipeline.py')
tpot.score(X[validation_indices], data.loc[validation_indices, 'class'].values)
```

I'm not going to worry too much about the warnings for now. Let's just make sure
we can save our best models for testing and go from there.

Now let's see if we can run this in Biowulf. I had to create a tpot environment
so make sure all packages were setup properly.

```bash
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate tpot
```

I removed xgboost package from conda to check whether it was repsonsible for
some of the pipelines breaking. There were some inconsistencies in DASK about
output parameters it received from CV results, so I'm checking whether that can
be originating from xgb.




# TODO
* determine random chance classifiers / regressors
* remove warnings?
* play with OD threshold
* construct distributions in the validation set
* run for all 6 classification and 6 regression tasks
* 