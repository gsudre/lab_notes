# 2020-03-09 09:36:08

For ML, Philip suggested I could try 2 different approaches: if the AUC is
impressive, I could try including all domains and just listing their
contributions. That's the approach I was doing before. Another option would be
to give the AUC of machines trained in each domain exclusively, and then try
higher order combinations of them, using age and gender all along.

Let's re-run our ML pipelines, keeping in mind that we don't need to run groups
3 and 4 using clinicals, only cd==1, and no need to run medication groups.

```bash
g=2
cd ~/data/baseline_prediction/prs_start
my_script=~/research_code/baseline_prediction/nonstacked_${g}group_dataImpute.R;
out_file=swarm.${g}group_impInter
rm $out_file
for clf in `cat all_clf.txt`; do
    for sx in inatt hi; do
        if [[ $g == 2 ]]; then
            for cm in T F; do
                echo "Rscript $my_script $sx $clf 1 $cm F ~/tmp/residsNOCO_${g}group_impInter.csv;" >> $out_file;
            done;
        else
            echo "Rscript $my_script $sx $clf 1 F F ~/tmp/residsNOCO_${g}group_impInter.csv;" >> $out_file;
        fi;
    done;
done

swarm -g 10 -t 1 --job-name inter${g} --time 20:00 -f $out_file \
    -m R --partition quick --logdir trash
```

and then:

```bash
g=2
cd ~/data/baseline_prediction/prs_start
my_script=~/research_code/baseline_prediction/stacked_${g}group_dataImpute.R;
out_file=swarm.${g}group_impStack
rm $out_file
for clf in hdda rda stepLDA glmStepAIC dwdLinear bayesglm earth LogitBoost \
    kernelpls cforest; do
    for ens in rpart glm glmStepAIC rpart2 C5.0Tree; do
        for sx in inatt hi; do
            if [[ $g == 2 ]]; then
                for cm in T F; do
                    echo "Rscript $my_script $sx $clf $ens 1 $cm F ~/tmp/residsNOCO_${g}group_impStack.csv;" >> $out_file;
                done;
            else
                echo "Rscript $my_script $sx $clf $ens 1 F F ~/tmp/residsNOCO_${g}group_impStack.csv;" >> $out_file;
            fi;
        done
    done;
done

swarm -g 10 -t 1 --job-name stack${g} --time 20:00 -f $out_file \
    -m R --partition quick --logdir trash
```

And time to check what the best results are:

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsNOCO_2group_impStack.csv', header=F)
colnames(res) = c('sx', 'model', 'ensemble', 'clin_diff', 'use_clinical',
                  'use_meds', 'num_groups', 'train_AUC', 'test_AUC')
for (clf in unique(res$model)) {
    for (ens in unique(res$ensemble)) {
        for (cd in unique(res$clin_diff)) {
            for (uc in unique(res$use_clinical)) {
                for (um in unique(res$use_meds)) {
                    idx = (res$model == clf & res$ensemble == ens &
                            res$clin_diff == cd & res$use_clinical == uc &
                            res$use_meds == um)
                    pos = which(idx)
                    if (length(pos) == 2) {
                        my_str = paste(c(clf, ens, cd, uc, um), collapse='_')
                        params = c(params, my_str)
                        scores = c(scores, mean(res[pos, 'test_AUC']))
                    }
                }
            }
        }
    }
}
a = sort(scores, decreasing=T, index.return=T)
print(params[a$ix[1]])
```

The stacked model gives us  "stepLDA_C5.0Tree_1_TRUE_FALSE":

```
> res[res$model=='stepLDA' & res$ensemble=='C5.0Tree' & res$clin_diff==1 & res$use_clinical==T & res$use_meds==F,]
       sx   model ensemble clin_diff use_clinical use_meds num_groups train_AUC
172    hi stepLDA C5.0Tree         1         TRUE    FALSE          2   0.87259
176 inatt stepLDA C5.0Tree         1         TRUE    FALSE          2   0.94769
    test_AUC
172 0.737037
176 0.766667
```

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsNOCO_2group_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'clin_diff', 'use_clinical',
                  'use_meds', 'num_groups', 'train_AUC', 'test_AUC')
for (clf in unique(res$model)) {
    for (cd in unique(res$clin_diff)) {
        for (uc in unique(res$use_clinical)) {
            for (um in unique(res$use_meds)) {
                idx = (res$model == clf &
                        res$clin_diff == cd & res$use_clinical == uc &
                        res$use_meds == um)
                pos = which(idx)
                if (length(pos) == 2) {
                    my_str = paste(c(clf, cd, uc, um), collapse='_')
                    params = c(params, my_str)
                    scores = c(scores, mean(res[pos, 'test_AUC']))
                }
            }
        }
    }
}
a = sort(scores, decreasing=T, index.return=T)
print(params[a$ix[1]])
```

Here our best is "glmnet_1_TRUE_FALSE":

```
> res[res$model=='glmnet' & res$clin_diff==1 & res$use_clinical==T & res$use_meds==F,]
      sx  model clin_diff use_clinical use_meds num_groups train_AUC test_AUC
25 inatt glmnet         1         TRUE    FALSE          2  0.758832 0.829630
88    hi glmnet         1         TRUE    FALSE          2  0.754821 0.766667
```

The glmnet classifier for inatt is all based on base_inatt, so that's not great.
Same thing for hi.

If we go for second place, we have: "cforest_1_TRUE_FALSE":

```
> res[res$model=='cforest' & res$clin_diff==1 & res$use_clinical==T & res$use_meds==F,]
      sx   model clin_diff use_clinical use_meds num_groups train_AUC test_AUC
95    hi cforest         1         TRUE    FALSE          2  0.979339 0.600000
97 inatt cforest         1         TRUE    FALSE          2  0.980978 0.911111
```

Would we have an argument for two different classifiers depending on sx? That's
hard to argue...

Third place is: "glmboost_1_TRUE_FALSE":

```
> res[res$model=='glmboost' & res$clin_diff==1 & res$use_clinical==T & res$use_meds==F,]
      sx    model clin_diff use_clinical use_meds num_groups train_AUC test_AUC
19    hi glmboost         1         TRUE    FALSE          2  0.896694 0.651852
20 inatt glmboost         1         TRUE    FALSE          2  0.930707 0.807407
```

The variable importance here is a bit more distributed, which is nice. What if
we choose based on the classifiers using clinical and not?

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsNOCO_2group_impStack.csv', header=F)
colnames(res) = c('sx', 'model', 'ensemble', 'clin_diff', 'use_clinical',
                  'use_meds', 'num_groups', 'train_AUC', 'test_AUC')
for (clf in unique(res$model)) {
    for (ens in unique(res$ensemble)) {
        idx = (res$model == clf & res$ensemble == ens)
        pos = which(idx)
        if (length(pos) == 4) {
            my_str = paste(c(clf, ens, cd), collapse='_')
            params = c(params, my_str)
            scores = c(scores, mean(res[pos, 'test_AUC']))
        }
    }
}
a = sort(scores, decreasing=T, index.return=T)
print(params[a$ix[1]])
```

Still "stepLDA_C5.0Tree_1" for stacked:

```
> res[res$model=='stepLDA' & res$ensemble=='C5.0Tree' & res$clin_diff==1,]
       sx   model ensemble clin_diff use_clinical use_meds num_groups train_AUC
162    hi stepLDA C5.0Tree         1        FALSE    FALSE          2  0.936639
167 inatt stepLDA C5.0Tree         1        FALSE    FALSE          2  0.938179
172    hi stepLDA C5.0Tree         1         TRUE    FALSE          2  0.872590
176 inatt stepLDA C5.0Tree         1         TRUE    FALSE          2  0.947690
    test_AUC
162 0.507407
167 0.666667
172 0.737037
176 0.766667
```

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsNOCO_2group_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'clin_diff', 'use_clinical',
                  'use_meds', 'num_groups', 'train_AUC', 'test_AUC')
for (clf in unique(res$model)) {
    for (cd in unique(res$clin_diff)) {
        for (uc in unique(res$use_clinical)) {
            for (um in unique(res$use_meds)) {
                idx = (res$model == clf)
                pos = which(idx)
                if (length(pos) == 4) {
                    my_str = paste(c(clf, cd), collapse='_')
                    params = c(params, my_str)
                    scores = c(scores, mean(res[pos, 'test_AUC']))
                }
            }
        }
    }
}
a = sort(scores, decreasing=T, index.return=T)
print(params[a$ix[1]])
```

And glmnet is still the winner:

```
> res[res$model=='glmnet' & res$clin_diff==1,]
      sx  model clin_diff use_clinical use_meds num_groups train_AUC test_AUC
25 inatt glmnet         1         TRUE    FALSE          2  0.758832 0.829630
88    hi glmnet         1         TRUE    FALSE          2  0.754821 0.766667
93    hi glmnet         1        FALSE    FALSE          2  0.771350 0.496296
98 inatt glmnet         1        FALSE    FALSE          2  0.902174 0.518519
```

Just playing around, some other classifiers are worth trying, such as pam:

[1] "inatt,pam,1,TRUE,FALSE,2,0.797554,0.929630"
[1] "hi,pam,1,TRUE,FALSE,2,0.805785,0.740741"

Class distribution is very off here, only using a couple variables, which of
couse include base_sx.

# 2020-03-10 09:13:08

I'm thinking about those two papers Philip sent out:

Multimodal Neuroimaging-based Prediction of Adult Outcomes in Childhood-onset
ADHD using Ensemble Learning Techniques

Prediction Models of Functional Outcomes for Individuals
in the Clinical High-Risk State for Psychosis
or With Recent-Onset Depression

And how they did their ensemble. Let's try to do something similar, just playing
with our 2-class case for now.

Let's run it in the cluster and see what we get:

```bash
g=2
cd ~/data/baseline_prediction/prs_start
my_script=~/research_code/baseline_prediction/stacked_${g}group_dataImpute_LOOCV.R;
out_file=swarm.${g}group_impStackLOOCV
rm $out_file
for clf in hdda rda stepLDA glmStepAIC dwdLinear bayesglm earth LogitBoost \
    kernelpls cforest; do
    for ens in rpart glm glmStepAIC rpart2 C5.0Tree plr; do
        for sx in inatt hi; do
            if [[ $g == 2 ]]; then
                for cm in T F; do
                    echo "Rscript $my_script $sx $clf $ens 1 $cm F ~/tmp/residsNOCO_${g}group_impStackLOOCV.csv;" >> $out_file;
                done;
            else
                echo "Rscript $my_script $sx $clf $ens 1 F F ~/tmp/residsNOCO_${g}group_impStackLOOCV.csv;" >> $out_file;
            fi;
        done
    done;
done

swarm -g 10 -t 1 --job-name stack${g} --time 2:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

That's using 5-fold 5x. I could also try 10-10 to check if results change much?
Also running some plr ensembles on an interactive node.

```bash
grep plr $out_file | parallel --max-args=1 -j 32 {1};
```

Maybe other classifier should also be tested... but let's wait on these results
first.

```bash
g=2
cd ~/data/baseline_prediction/prs_start
my_script=~/research_code/baseline_prediction/stacked_${g}group_dataImpute_LOOCV1010.R;
out_file=swarm.${g}group_impStackLOOCV
rm $out_file
for clf in hdda rda stepLDA glmStepAIC dwdLinear bayesglm earth LogitBoost \
    kernelpls cforest; do
    for ens in rpart glm glmStepAIC rpart2 C5.0Tree plr; do
        for sx in inatt hi; do
            if [[ $g == 2 ]]; then
                for cm in T F; do
                    echo "Rscript $my_script $sx $clf $ens 1 $cm F ~/tmp/residsNOCO_${g}group_impStackLOOCV1010.csv;" >> $out_file;
                done;
            else
                echo "Rscript $my_script $sx $clf $ens 1 F F ~/tmp/residsNOCO_${g}group_impStackLOOCV1010.csv;" >> $out_file;
            fi;
        done
    done;
done

swarm -g 10 -t 1 --job-name stack${g} --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```
