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

# 2019-11-20 09:28:26

OK, so I'll run a few experiment, first with AD and RD, to see what variation of
FPR works best. We'll sprinkle in PCA before or after, along with no PCA. Note
that I had to remove the univariate selector from the pipeline when doing PCA
first because more often than not PCA components were not associated with the
outcome, even at .1.

target (AD)| raw_train | raw_test | PCAfirst_train | PCAfirst_test | PCA_last train| PCA_last test
--- | --- | --- | --- | --- | --- | --- |
SX_inatt_groupStudy | 0.61 (0.11) | 0.49 (0.32, 0.66) |  0.60 (0.11) | 0.46 (0.30, 0.63) | 0.61 (0.11) | 0.49 (0.32, 0.65) |
SX_HI_groupStudy |0.58 (0.14)|0.41 (0.24, 0.58)|0.47 (0.14)|0.68 (0.51, 0.86)|0.57 (0.14)|0.41 (0.24, 0.58)

target (RD)| raw_train | raw_test | PCAfirst_train | PCAfirst_test | PCA_last train| PCA_last test
--- | --- | --- | --- | --- | --- | --- |
SX_inatt_groupStudy | 0.63 (0.12) | 0.71 (0.57, 0.86) |  0.58 (0.12) | 0.77 (0.64, 0.90) | 0.63 (0.12) | 0.71 (0.56, 0.86)
SX_HI_groupStudy |0.52 (0.12)|0.36 (0.15, 0.57)!!WRONG LABELS|0.47 (0.12)|0.38 (0.17, 0.58)|0.51 (0.12)|0.65 (0.45, 0.85)|

Results not looking great. Slight preference for PCA first result.

I'll run similar tests using LR, while I play again with TPOT. again, removing
FPR selector in PCAfirst case:

target (AD)| raw_train | raw_test | PCAfirst_train | PCAfirst_test | PCA_last train| PCA_last test
--- | --- | --- | --- | --- | --- | --- |
SX_inatt_groupStudy|0.59 (0.11)|0.61 (0.45, 0.77)|0.61 (0.10)|0.54 (0.37, 0.70)|0.62 (0.11)|0.57 (0.41, 0.73)

target (RD)| raw_train | raw_test | PCAfirst_train | PCAfirst_test | PCA_last train| PCA_last test
--- | --- | --- | --- | --- | --- | --- |
SX_inatt_groupStudy | ||0.59 (0.11)|0.76 (0.63, 0.89)|0.64 (0.12)|0.69 (0.54, 0.84)

While those are running, start descriptive results!!!


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
min_time = 30*9
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
data_base = data_base[!is.na(data_base[,'SX_HI_slopeNext']), ]
for (s in c('Next', 'Last', 'Study')) {
    for (t in c('SX_inatt', 'SX_HI')) {
        data_base[, sprintf('%s_group%s', t, s)] = NA
        idx = data_base[, sprintf('%s_slope%s', t, s)] < 0
        data_base[idx, sprintf('%s_group%s', t, s)] = 'improvers'
        data_base[!idx, sprintf('%s_group%s', t, s)] = 'nonimprovers'
    }
}
today = format(Sys.time(), "%m%d%Y")
out_fname = sprintf('~/data/baseline_prediction/cog_%s', today)
write.csv(data_base, file=sprintf('%s.csv', out_fname), row.names=F, na='', 
          quote=F)
```

Then, for adhdDX and lastPersistent, we do it very similarly to above, except
that we do this instead:

```r
adhd_clin = clin[clin$source != 'DICA_on', ]
# redo everything else, but the ending is
data_base[, 'adhdDX'] = NA
data_base[, 'threeWay'] = NA
for (r in 1:nrow(data_base)) {
    subj = data_base[r,]$MRN
    subj_clin = adhd_clin[which(adhd_clin$MRN==subj), ]
    clin_dates = as.Date(as.character(subj_clin$DOA), format="%m/%d/%y")
    dob = as.Date(as.character(data_base[r, 'DOB']),
                  format="%m/%d/%Y")
    age_clinical = as.numeric((clin_dates - dob)/365.25)

    ordered_ages = sort(age_clinical, index.return=T)
    subj_clin = subj_clin[ordered_ages$ix, ]
    cur_age_clin = as.numeric((as.Date(as.character(data_base[r,]$DOA),
                                       format="%m/%d/%Y") - dob)/365.25)
    # in case there are duplicates, take the last one
    date_idx = max(which(ordered_ages$x==cur_age_clin))
    # strictly DSM5 categorization at the time of scan
    if ((subj_clin[date_idx, 'SX_inatt'] >= 6) ||
        (subj_clin[date_idx, 'SX_HI'] >= 6)) {
        data_base[r, 'adhdDX'] = 'improvers'  # just for back-compatibility
        if ((subj_clin[nrow(subj_clin), 'SX_inatt'] >= 6) ||
            (subj_clin[nrow(subj_clin), 'SX_HI'] >= 6)) {
            data_base[r, 'threeWay'] = 'per'
        } else {
            data_base[r, 'threeWay'] = 'rem'
        }
    } else {
        data_base[r, 'adhdDX'] = 'nonimprovers'
        data_base[r, 'threeWay'] = 'NV'
    }
}
today = format(Sys.time(), "%m%d%Y")
out_fname = sprintf('~/data/baseline_prediction/cog_baseDX_%s', today)
write.csv(data_base, file=sprintf('%s.csv', out_fname), row.names=F, na='', 
          quote=F)
```

And for DSM5Outcome:

```r
keep_me = c()
for (s in unique(clin$MRN)) {
    subj_idx = which(clin$MRN==s & !is.na(clin$SX_inatt) & !is.na(clin$SX_HI))
        if (length(subj_idx) > 0) {
        subj_data = clin[subj_idx, ]
        dates = as.Date(as.character(subj_data$DOA), format="%m/%d/%y")
        base_DOA = which.min(dates)
        if (subj_data[base_DOA,]$SX_HI >= 6 || subj_data[base_DOA,]$SX_inat >= 6) {
            keep_me = c(keep_me, subj_idx)
        }
    }
}
adhd_clin = clin[keep_me, ]
# everything the same until
data_base[, c('lastPersistent')] = NA
for (r in 1:nrow(data_base)) {
    subj = data_base[r,]$MRN
    subj_clin = adhd_clin[which(adhd_clin$MRN==subj), ]
    clin_dates = as.Date(as.character(subj_clin$DOA), format="%m/%d/%y")
    dob = as.Date(as.character(data_base[r, 'DOB']),
                  format="%m/%d/%Y")
    age_clinical = as.numeric((clin_dates - dob)/365.25)

    ordered_ages = sort(age_clinical, index.return=T)
    subj_clin = subj_clin[ordered_ages$ix, ]
    cur_age_clin = as.numeric((as.Date(as.character(data_base[r,]$DOA),
                                       format="%m/%d/%Y") - dob)/365.25)
    # in case there are duplicates, take the last one
    date_idx = max(which(ordered_ages$x==cur_age_clin))
    
    if ((subj_clin[nrow(subj_clin), 'SX_inatt'] >= 6) ||
        (subj_clin[nrow(subj_clin), 'SX_HI'] >= 6)) {
        data_base[r, 'lastPersistent'] = 'improvers'  # just for back-compatibility
    } else {
        data_base[r, 'lastPersistent'] = 'nonimprovers'
    }
}
today = format(Sys.time(), "%m%d%Y")
out_fname = sprintf('~/data/baseline_prediction/cog_DSM5Outcome_%s', today)
write.csv(data_base, file=sprintf('%s.csv', out_fname), row.names=F, na='', 
          quote=F)
```

# 2019-11-21 09:29:00

The TPOT results are taking a long time to run, but at least the training CV is
looking promising, at .81 so far, with AD HI_groupStudy. Not sure if it's just overfitting everything,
but I wonder if I can make it go faster by playing with OMP thread. FYI, the MDS
pipeline is not working at all, so that's using the regular classify pipeline.

But it's taking a veeery long time. So, I'll swarm it and keep track offline:

```bash
# bw
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate tpot
export OMP_NUM_THREADS=1
cd ~/data/baseline_prediction/tpot_swarms

jname=vox42;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/tpot_classify.py;
res=~/data/baseline_prediction/tpot_results_42;
s=42;
for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
    phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
    for i in Next Last Study; do
        for j in SX_inatt SX_HI; do
            target=${j}_group${i};
            echo "python $code $phen $target $res $s | tee -a ${res}/${p}_${target}.txt" >> $swarm_file;
        done;
    done;
    target=lastPersistent;
    phen2=~/data/baseline_prediction/${p}_OD0.95_DSM5Outcome_11052019.csv;
    echo "python $code $phen2 $target $res $s | tee -a ${res}/${p}_${target}.txt" >> $swarm_file;
    target=adhdDX;
    phen2=~/data/baseline_prediction/${p}_OD0.95_baseDX_11072019.csv;
    echo "python $code $phen2 $target $res $s | tee -a ${res}/${p}_${target}.txt" >> $swarm_file;
done;
phen2=~/data/baseline_prediction/clinics_binary_baseDX_11072019.csv;
echo "python $code $phen2 $target $res $s | tee -a ${res}/clinics_binary_${target}.txt" >> $swarm_file;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 --logdir=trash_${jname} \
    --job-name ${jname} --time=4-00:00:00 --merge-output --partition norm
```

So, let's see how long this takes to run. I'm also running some option locally
and interactively just in case.

# 2019-11-22 11:28:20

While I'm waiting for all that stuff above finish, I'm running the cog domain
locally, and generating the degree centrality files. I decided to go with that
for fMRI because it will be easier to clusterize later for the descriptives
analysis. Let's hope there's something there, though. Otherwise I'll have to go
with the connectivity metric I used for the slope heritability paper.

In the meanwhile, I'm trying a more conservative version of TPOT that only has
LinearSVC and LR in its pipeline. Not even any selection or transformation.
Shouldn't take very long to run. I'm running that interactively because even
though the CV results for TPOT were looking good, test results weren't,
indicating overfitting. I might add some transformations or selection later,
depending on the simpler results. I'm basically just trimming down the default
classification setup I found in
https://github.com/EpistasisLab/tpot/blob/master/tpot/config/classifier.py

I can actually swarm that to make things run faster:

```bash
# bw
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate tpot
export OMP_NUM_THREADS=1
cd ~/data/baseline_prediction/tpot_swarms

jname=TPOTcon;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/tpot_classify_conservative.py;
res=~/data/baseline_prediction/tpot_results_con;
s=42;
for tpt in "Selector-Transformer-Classifier" "Transformer-Selector-Classifier" "Classifier"; do
    for i in Next Last Study; do
        for j in SX_inatt SX_HI; do
            target=${j}_group${i};
            for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
                phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
                echo "python $code $phen $target $res $s $tpt | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
            done;
            phen=~/data/baseline_prediction/cog_11202019.csv;
            echo "python $code $phen $target $res $s $tpt | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
        done;
        target=lastPersistent;
        for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
            phen2=~/data/baseline_prediction/${p}_OD0.95_DSM5Outcome_11052019.csv;
            echo "python $code $phen2 $target $res $s $tpt  | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
        done
        phen2=~/data/baseline_prediction/cog_DSM5Outcome_11202019.csv;
        echo "python $code $phen2 $target $res $s $tpt  | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
        target=adhdDX;
        for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
            phen2=~/data/baseline_prediction/${p}_OD0.95_baseDX_11072019.csv;
            echo "python $code $phen2 $target $res $s $tpt  | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
        done
        phen2=~/data/baseline_prediction/cog_baseDX_11202019.csv;
        echo "python $code $phen2 $target $res $s $tpt  | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
    done;
    phen2=~/data/baseline_prediction/clinics_binary_baseDX_11072019.csv;
    echo "python $code $phen2 $target $res $s $tpt | tee -a ${res}/clinics_binary_${target}_${tpt}.txt" >> $swarm_file;
done
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 --logdir=trash_${jname} \
    --job-name ${jname} --time=2-00:00:00 --merge-output --partition norm
```

How about using generalized estimator (LOOCV optimized), within training set?

For now, I'll run the rsfmri results as well to see what version of TPOT works
best, and what are our best results within that modality:

```bash
# bw
source /data/$USER/conda/etc/profile.d/conda.sh
conda activate tpot
export OMP_NUM_THREADS=1
cd ~/data/baseline_prediction/tpot_swarms

jname=fmriCon;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/tpot_classify_conservative.py;
res=~/data/baseline_prediction/tpot_results_con;
s=42;
for tpt in "Selector-Transformer-Classifier" "Transformer-Selector-Classifier" "Classifier"; do
    for i in Next Last Study; do
        for j in SX_inatt SX_HI; do
            target=${j}_group${i};
            for p in s sz; do
                phen=~/data/baseline_prediction/rsfmri_${p}DC_OD0.95_11222019.csv;
                echo "python $code $phen $target $res $s $tpt | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
            done;
        done;
        target=lastPersistent;
        for p in s sz; do
            phen2=~/data/baseline_prediction/rsfmri_${p}DC_OD0.95_DSM5Outcome_11222019.csv;
            echo "python $code $phen2 $target $res $s $tpt  | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
        done
        target=adhdDX;
        for p in s sz; do
            phen2=~/data/baseline_prediction/rsfmri_${p}DC_OD0.95_baseDX_11222019.csv;
            echo "python $code $phen2 $target $res $s $tpt  | tee -a ${res}/${p}_${target}_${tpt}.txt" >> $swarm_file;
        done
    done;
done
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 --logdir=trash_${jname} \
    --job-name ${jname} --time=2-00:00:00 --merge-output --partition norm
```

In the meanwhile I'm also trying 2vs2 ideas. I'm running a couple using XGB in
interactive sessions, but that takes a while to run, so let me see if I can have
a few linear solutions with a smaller search grid.

Nothing there... I'm doing about 56% with LR in RD, SX_HI_study... Not much
better using XGB. DSM5Outcome didn't work either, nor did adhdDX.

# TODO:
* check if there's anything good in conservative TPOT results. at least we
  should be able to decide what targets to use for descriptive work?
* 2vs2 idea XGB results?
* combine spaces through voting to see if there is any sort of improvement.
* what if my semantic space is the symptom counts? We'd have 18 different
  features we're predicting, pottentially binary (but not necessarily).
* try age prediction again, based on models established in the literature?