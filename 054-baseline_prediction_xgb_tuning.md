# 2019-11-15 13:44:51

Let's do some tunning for SX_inatt_groupStudy, RD:

```python
#Import libraries:
import pandas as pd
import numpy as np
import xgboost as xgb
from xgboost.sklearn import XGBClassifier
from sklearn import metrics
from sklearn.model_selection import GridSearchCV, train_test_split
import os

home = os.path.expanduser('~')
phen_fname = home + '/data/baseline_prediction/dti_rd_OD0.95_11052019.csv'
target = 'SX_inatt_groupStudy'
myseed = 42

# 16 if running it locally
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', '16'))

data = pd.read_csv(phen_fname)

# remove columns that are all NaNs
data.dropna(axis=1, how='all', inplace=True)

data.rename(columns={target: 'class'}, inplace=True)
data['class'] = data['class'].map({'improvers': 1, 'nonimprovers': 0})
scoring='roc_auc'

# it's a feature to be used if it starts with v_
predictors = [f for f in data.columns if f.find('v_') == 0]

y = data['class'].values
training_indices, testing_indices = train_test_split(data.index,
                                                            stratify = y, train_size=0.65,
                                                            test_size=0.35,
                                                            random_state=myseed)

train = data.iloc[training_indices, :]
```

Now we define the reporting function:

```python
def modelfit(alg, dtrain, predictors,useTrainCV=True, cv_folds=3, early_stopping_rounds=50):
    if useTrainCV:
        xgb_param = alg.get_xgb_params()
        xgtrain = xgb.DMatrix(dtrain[predictors].values, label=dtrain['class'].values)
        cvresult = xgb.cv(xgb_param, xgtrain, num_boost_round=alg.get_params()['n_estimators'], nfold=cv_folds,
            metrics='auc', early_stopping_rounds=early_stopping_rounds)
        alg.set_params(n_estimators=cvresult.shape[0])
    # #Fit the algorithm on the data
    # alg.fit(dtrain[predictors], dtrain['class'],eval_metric='auc')
    # #Predict training set:
    # dtrain_predictions = alg.predict(dtrain[predictors])
    # dtrain_predprob = alg.predict_proba(dtrain[predictors])[:,1]
    #Print model report:
    print("\nModel Report")
    # print("Accuracy : %.4g" % metrics.accuracy_score(dtrain['class'].values, dtrain_predictions))
    # print "AUC Score (Train): %f" % metrics.roc_auc_score(dtrain['Disbursed'], dtrain_predprob)
    idx = np.argmax(cvresult.values[:, 2])
    print("AUC Score (Train): %.4f +- %.4f" % (cvresult.values[idx, 2], cvresult.values[idx, 3]))
    print('n_estimators = %d' % cvresult.shape[0])
    alg.fit(dtrain[predictors], dtrain['class'],eval_metric='auc')
    feat_imp = pd.Series(alg.get_booster().get_fscore()).sort_values(ascending=False)
    print(feat_imp.head(10))
    # feat_imp.plot(kind='bar', title='Feature Importances')
    # plt.ylabel('Feature Importance Score')
```

Let's start with a very basic model just to find a good learning rate and number
of estimator:

```python
weight = (np.sum(y[training_indices]==0) / np.sum(y[training_indices]==1))
xgb1 = XGBClassifier(
 learning_rate =0.1,
 n_estimators=1000,
 max_depth=5,
 min_child_weight=1,
 gamma=0,
 subsample=0.8,
 colsample_bytree=0.8,
 objective= 'binary:logistic',
 nthread=ncpus,
 scale_pos_weight=weight,
 seed=myseed)
modelfit(xgb1, train, predictors)
```

```
Model Report
AUC Score (Train): 0.5792 +- 0.0905
n_estimators = 108
v_v07810    8
v_v03553    7
v_v10932    5
v_v05059    5
v_v06784    5
v_v02649    3
v_v11657    3
v_v01281    3
v_v00711    3
v_v10931    3
```

Now, using that learning rate and the number of estimators, let's see how well
we can do by tweaking other parameters. Note that we can always change the
learning rate and get a different number of estimators... it'll just take longer
to run.

```python
param_test1 = {
 'max_depth':range(3,10,2),
 'min_child_weight':range(1,6,2)
}
gsearch1 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, n_estimators=108,
 gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=ncpus, scale_pos_weight=weight, seed=myseed),
 param_grid = param_test1, scoring='roc_auc',n_jobs=1,iid=False, cv=3)
gsearch1.fit(train[predictors],train['class'])
gsearch1.best_params_, gsearch1.best_score_
```

We're up oto .589 using {'max_depth': 3, 'min_child_weight': 5}. Let's go more
fine grained:

```python
param_test1 = {
 'max_depth':range(3,6,1),
 'min_child_weight':range(1,10,1)
}
gsearch1 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, n_estimators=108,
 gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=ncpus, scale_pos_weight=weight, seed=myseed),
 param_grid = param_test1, scoring='roc_auc',n_jobs=1,iid=False, cv=3)
gsearch1.fit(train[predictors],train['class'])
gsearch1.best_params_, gsearch1.best_score_
```

Nope, that's still our best. Let's see if we get anything by tuning gamma:

```python
param_test3 = {
 'gamma':[i/10.0 for i in range(0,5)]
}
gsearch3 = GridSearchCV(estimator = XGBClassifier( learning_rate =0.1, n_estimators=108, max_depth=3,
 min_child_weight=5, gamma=0, subsample=0.8, colsample_bytree=0.8,
 objective= 'binary:logistic', nthread=ncpus, scale_pos_weight=weight,seed=myseed), 
 param_grid = param_test3, scoring='roc_auc',n_jobs=1,iid=False, cv=3)
gsearch3.fit(train[predictors],train['class'])
gsearch3.best_params_, gsearch3.best_score_
```

It turns out that gamma 0 is still our best option here... the tutorial at this
point re-runs the first step to check on the optimal number of estimators...
see, I don't like this, because there's an optimal parameter combination that is
not a linear search. One depends on the other. So, it's better to have a search
space that we can navigate through... just like what we're currently doing! 

Maybe we just need to tweak those values a bit to have a good enough search
space.

There I'm starting with:

```python
params = {
    "clf__colsample_bytree": uniform(0.9, 0.1),
    "clf__gamma": uniform(0, 0.5),
    "clf__learning_rate": uniform(0.1, 0.2), # default 0.1 
    "clf__max_depth": randint(3, 6), # default 3
    "clf__n_estimators": randint(100, 150), # default 100
    "clf__subsample": uniform(0.6, 0.2),
}
```

So, let's tweak that a bit:

```python
params = {
    "clf__colsample_bytree": uniform(0.5, 0.45),
    "clf__gamma": uniform(0, 0.5),
    "clf__learning_rate": uniform(0.01, 0.2), # default 0.1 
    "clf__max_depth": randint(3, 6), # default 3
    "clf__n_estimators": randint(100, 150), # default 100
    "clf__subsample": uniform(0.6, 0.3),
    'clf__min_child_weight': randint(1,10),
    'clf__scale_pos_weight': [1, (np.sum(y[training_indices]==0) / np.sum(y[training_indices]==1))],
    'clf__reg_lamba': [1e-5, 1e-2, 0.1, 1, 100],
    'clf__reg_alpha': [1e-5, 1e-2, 0.1, 1, 100],
}
```

This gave me .656 in training and .660 in testing, so we're learning
something... maybe I'll just switch to this now and get actual results on
everything. Then, it's just a matter of combining datasets and running
descriptives.

FYI, I tried increasing the number of iterations to 1000, but the best result
was still the same...



