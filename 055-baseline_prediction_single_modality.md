# 2019-11-19 11:01:38

I think we're at a point to establish single-modality results. I'll run them
interactively, create a table with the results, and save the predictions and
models as well. Here we go:

```python
cd ~/data/baseline_prediction
code = '/home/sudregp/research_code/baseline_prediction/xgb_classifier_65-35.py'
res = '/home/sudregp/data/baseline_prediction/manual_results_xgb_tuned'
s = 42
phen_file = '~/data/baseline_prediction/dti_fa_OD0.95_11052019.csv'
%run $code $phen_file SX_inatt_groupStudy $res $s

import joblib
joblib.dump(my_search, '/home/sudregp/data/baseline_prediction/xgb_models/fa_inatt.sav')
```

And then I ran the groupStudy variables along with DX and lastPersistent for all
phenotypes.

To calculate the confidence interval, I'm just reloading the model. Just note
that the AUC here is calculated a bit differently than what we're getting
before, so the results might be slightly different:

```python
# just run it once
%run ~/research_code/baseline_prediction/proc_ci.py
import joblib
import pandas as pd
from sklearn.model_selection import train_test_split

phen_file = '/home/sudregp/data/baseline_prediction/dti_fa_OD0.95_11052019.csv'
s = 42
target = 'SX_inatt_groupStudy'
data = pd.read_csv(phen_file)
data.dropna(axis=1, how='all', inplace=True)
data.rename(columns={target: 'class'}, inplace=True)
data['class'] = data['class'].map({'improvers': 1, 'nonimprovers': 0})
feature_names = [f for f in data.columns if f.find('v_') == 0]
X = data[feature_names].values
y = data['class'].values
training_indices, testing_indices = train_test_split(data.index,
                                                    stratify = y, train_size=0.65,
                                                    test_size=0.35,
                                                    random_state=s)

my_search = joblib.load('/home/sudregp/data/baseline_prediction/xgb_models/fa_inatt.sav')
y_probs = my_search.predict_proba(X[testing_indices])
ypos_prob = np.array([i[1] for i in y_probs])
auc, auc_cov = delong_roc_variance(y[testing_indices], ypos_prob)
auc_std = np.sqrt(auc_cov)
lower_upper_q = np.abs(np.array([0, 1]) - (1 - alpha) / 2)
ci = stats.norm.ppf(lower_upper_q, loc=auc, scale=auc_std)
ci[ci > 1] = 1
print('AUC test (95pct CI): %.2f (%.2f, %.2f)' % (auc, ci[0], ci[1]))
```

The code above works well when we need to retrieve models, but it doesn't make
sense to run things twice when that's not necessary, so I added bits of it to
the actual decoding script, so that we have the 95% CI right of the bat.

The table below summarizes the current results. I was also looking for patterns
to let me know if I needed to increase or decrease the search space.

pheno | inatt train | inatt test | HI train | HI test | DX train | DX test | outcome train | outcome test
--- | --- | --- | --- | --- | --- | --- | --- | --- |
FA | 
AD | 
RD | 

Now that I fixed the script to use ROC for scoring as well, looks like I'm
really overfitting? Might be worth trying it with the PCA code again.

Or even go back to logreg or SVM... at least it's a simpler model, with a smaller
search space. I'd prefer not to do that so we can remain robust to NAs, but we
might not have a choice.Maybe we could have an ensemble on top of linear
classifiers within modalities, instead of just majority voting, or combined
probability voting. But we need to be careful not to overfit.

Yeah, even with RD we're at .65 testing in inatt and .6 for HI. Using PCA was
actually worse. And that's RD, our best result.

I'm gonna go back to linear models for now. Maybe logreg and SVC, with FPR and
without it, and then sprinkle some PCA on those as well. Let's see what we get.
And then continue the stuff I'm doing to include cognitive, rsfmri and PRS
below. Should probably do ADHD-200 at one point as well. 



## Cognitive

I need to prepare the cognitive domain data. But I'll just go ahead and use
whatever I had accumulated from before, and make sure the variables match our
current framework:

```r
cog = read.csv('~/data/baseline_prediction/cog_baseline_merge_clean.csv')
vars = c(# cpt
         "N_of_omissions", "PC_omissions", "N_commissions", "PC_commissions",
         "hit_RT", "hit_RT_SE", "variability_of_SE", "D.prime", "beta",
         "N_perservations", "PC_perservations", "hit_RT_block_change",
         "hit_RT_SE_block_change", "hit_RT_ISI_change", "hit_RT_SE_ISI_change",
         # WISC raw
         "RAW_DS_total", "RAW_DSB", "Raw_DSF", "Raw_SS_total", "Raw_SSB",
         "Raw_SSF",
         # WISC STD
         "Standard_DS_total", "Standard_DSB", "Standard_DSF",
         "Standard_SS_total", "Standard_SSB", "Standard_SSB.1",
         # wj
         "PS", "SSDS_WJ", "SSVM_WJ",
         # iq
         "FSIQ")
# pick the earliest day for each subject to compare with clinicals
date_vars = c('DOA_CPT', 'date_WISC', 'date_WJ', 'date_IQ')
cog$base_date = NA
for (r in 1:nrow(cog)) {
    subj_dates = c()
    for (v in date_vars) {
        subj_dates = c(subj_dates, as.Date(as.character(cog[r, v]),
                                           format="%m/%d/%y"))
    }
    idx = which.min(subj_dates)
    if (length(idx) > 0) {
        cog[r, 'base_date'] = as.character(cog[r, date_vars[idx]])
    }
}
cog = cog[!is.na(cog$base_date), ]

clin = read.csv('~/data/baseline_prediction/clinical_09182019_clean.csv')
# make sure the SX columns are numeric!
clin$SX_HI = as.numeric(as.character(clin$SX_hi))
clin$SX_inatt = as.numeric(as.character(clin$SX_inatt))
# let's remove on-medication items to not confuse the slopes
clin = clin[clin$source != 'DICA_on', ]

# we'll keep anyone that starts off with 4 or more symptoms
keep_me = c()
for (s in unique(clin$MRN)) {
    subj_idx = which(clin$MRN==s & !is.na(clin$SX_inatt) & !is.na(clin$SX_HI))
        if (length(subj_idx) > 0) {
        subj_data = clin[subj_idx, ]
        dates = as.Date(as.character(subj_data$DOA), format="%m/%d/%y")
        base_DOA = which.min(dates)
        if (subj_data[base_DOA,]$SX_HI >= 4 || subj_data[base_DOA,]$SX_inat >= 4) {
            keep_me = c(keep_me, subj_idx)
        }
    }
}
adhd_clin = clin[keep_me, ]

# crop imaging data to only keep the people with cognitive data
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = mergeOnClosestDate(cog, adhd_clin,
                        unique(adhd_clin$MRN),
                         x.date='base_date',
                         x.id='MRN')
# changing variable names to make it easier to find them later
new_names = c()
for (v in colnames(df)) {
    if (v %in% vars) {
        new_names = c(new_names, sprintf('v_%s', v))
    } else {
        new_names = c(new_names, v)
    }
}
colnames(df) = new_names
data_base = df

data_base[, c('SX_HI_slopeNext', 'SX_inatt_slopeNext',
              'SX_HI_slopeLast', 'SX_inatt_slopeLast',
              'SX_HI_slopeStudy', 'SX_inatt_slopeStudy')] = NA
for (r in 1:nrow(data_base)) {
    subj = data_base[r,]$MRN
    subj_clin = adhd_clin[which(adhd_clin$MRN==subj), ]
    clin_dates = as.Date(as.character(subj_clin$DOA), format="%m/%d/%y")
    dob = as.Date(as.character(data_base[r, 'DOB']), format="%m/%d/%Y")
    age_clinical = as.numeric((clin_dates - dob)/365.25)

    ordered_ages = sort(age_clinical, index.return=T)
    subj_clin = subj_clin[ordered_ages$ix, ]
    cur_age_clin = as.numeric((as.Date(as.character(data_base[r,]$DOA),
                                       format="%m/%d/%Y") - dob)/365.25)
    # in case there are duplicates, take the last one
    date_idx = max(which(ordered_ages$x==cur_age_clin))
        
    for (t in c('SX_inatt', 'SX_HI')) {
        # the easiest one is overall study slope
        fm_str = sprintf('%s ~ age_clinical', t)
        fit = lm(as.formula(fm_str), data=subj_clin, na.action=na.exclude)
        data_base[r, sprintf('%s_slopeStudy', t)] = coefficients(fit)[2]

        # next time point has to be at least min_time from current scan
        cur = 1
        while ((date_idx + cur) < nrow(subj_clin) &&
               ((ordered_ages$x[date_idx + cur] - ordered_ages$x[date_idx]) < min_time/365.25)) {
            cur = cur + 1
        }
        slope = ((subj_clin[date_idx + cur, t] - subj_clin[date_idx, t]) / 
                 (ordered_ages$x[date_idx + cur] - ordered_ages$x[date_idx]))
        data_base[r, sprintf('%s_slopeNext', t)] = slope

        slope = ((subj_clin[nrow(subj_clin), t] - subj_clin[date_idx, t]) / 
                 (ordered_ages$x[nrow(subj_clin)] - ordered_ages$x[date_idx]))
        data_base[r, sprintf('%s_slopeLast', t)] = slope
    }
}
for (s in c('Next', 'Last', 'Study')) {
    for (t in c('SX_inatt', 'SX_HI')) {
        data_base[, sprintf('%s_group%s', t, s)] = NA
        idx = data_base[, sprintf('%s_slope%s', t, s)] < 0
        data_base[idx, sprintf('%s_group%s', t, s)] = 'improvers'
        data_base[!idx, sprintf('%s_group%s', t, s)] = 'nonimprovers'
    }
}
if (with_qc) {
    # change the names of QC variables, sex, and age to be included in
    # prediction
    for (v in c(qc_vars, 'age_at_scan')) {
        cidx = which(colnames(data_base) == v)
        colnames(data_base)[cidx] = sprintf('v_%s', v)
    }
    # make sure all non-binary variables are in the same scale
    my_vars = grepl(colnames(data_base), pattern='^v_')
    data_base[my_vars] = scale(data_base[my_vars])
    data_base$v_isMale = 0
    data_base[data_base$Sex == 'Male',]$v_isMale = 1 
    suffix = 'withQC_'
} else {
    suffix = ''
}
today = format(Sys.time(), "%m%d%Y")
out_fname = sprintf('~/data/baseline_prediction/struct_%s_OD%.2f_%s%s', prop,
                    qtile, suffix, today)
write.csv(data_base, file=sprintf('%s.csv', out_fname), row.names=F, na='', 
          quote=F)

```



# TODO:
* compute 95% CI using saved models and results
* try age prediction again, based on models established in the literature?