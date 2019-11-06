# 2019-11-04 14:45:06

Let's then explore the univariate filter idea for baseline prediction. If
anything, it'll better resemble usual neuroimaging papers.

Teh idea is that we'll split the data into training, validation, and testing. We
will pick univariately inside training, and do multiple combinations of
linearSVC, SVC, or Ensembles with the top X features.

# 2019-11-05 09:21:28

Her eI'm doing PCA, then univariate:

```bash
# locally for now
code=~/research_code/baseline_prediction/univariate_classifier.py;
phen=~/data/baseline_prediction/dti_JHUtracts_ADRDonly_OD0.95.csv;
vars=~/data/baseline_prediction/ad_rd_vars.txt;
res=~/data/tmp/;
for i in Next Last Study; do
    for j in SX_inatt SX_HI; do
        for s in `cat random25.txt`; do
            python3 $code $phen ${j}_group${i} $vars $res $s;
        done;
    done;
done;
```

And of course it might be the case that other data domains will do better. So,
let's try them in the cluster. But first, let's try to do the PCA oin the
training set only.

I made a few changes to the pipeline. Now we get the variable names
automatically, and do the PCA after univariate selection, to try to disembiguate
clusters of significance. I also only had 56 paremeters to optimize, so I went
for a grid search instead of randomized. Let's see how that goes:

```bash
# bw
export OMP_NUM_THREADS=1
code=~/research_code/baseline_prediction/univariate_classifier.py;
phen=~/data/baseline_prediction/dti_rd_OD0.95_11052019.csv;
res=~/data/tmp/;
for i in Next Last Study; do
    for j in SX_inatt SX_HI; do
        for s in `cat random25.txt`; do
            python3 $code $phen ${j}_group${i} $res $s;
        done;
    done;
done;
```

I also needed new functions:

```bash
Rscript ~/research_code/baseline_prediction/prep_dti_voxel_data.R
Rscript ~/research_code/baseline_prediction/prep_dti_voxel_data_DSM5Outcome.R
Rscript ~/research_code/baseline_prediction/prep_struct_voxel_data.R
Rscript ~/research_code/baseline_prediction/prep_struct_voxel_data_DSM5Outcome.R
```

Each run of DTI voxel is taking about 2min in the cluster. So, say 5min to be
safe. Which means we can fit 48 in 4h. We can do it all using:

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=dti;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/univariate_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in fa ad rd JHUtracts_ADRDonly; do
    phen=~/data/baseline_prediction/dti_${p}_OD0.95_11052019.csv;
    for s in `cat ../random25.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/dti_${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

And for structural we run something similar:

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=struct;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/univariate_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in area volume thickness; do
    phen=~/data/baseline_prediction/struct_${p}_OD0.95_11052019.csv;
    for s in `cat ../random25.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/struct_${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

I'm also trying an RFE classifier. Let's see what we get there:

```bash
# local
code=~/research_code/baseline_prediction/rfe_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in area volume thickness; do
    phen=~/data/baseline_prediction/struct_${p}_OD0.95_11052019.csv;
    for s in `cat random25.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                python3 $code $phen ${j}_group${i} $res $s;
            done;
        done;
        phen2=~/data/baseline_prediction/struct_${p}_OD0.95_DSM5Outcome_11052019.csv;
        python3 $code $phen2 lastPersistent $res $s;
    done;
done;
```

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=dtiFPR;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/univariate_fpr_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in fa ad rd JHUtracts_ADRDonly; do
    phen=~/data/baseline_prediction/dti_${p}_OD0.95_11052019.csv;
    for s in `cat ../random25.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/dti_${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

And for structural we run something similar:

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=structFPR;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/univariate_fpr_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in area volume thickness; do
    phen=~/data/baseline_prediction/struct_${p}_OD0.95_11052019.csv;
    for s in `cat ../random25.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/struct_${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

```bash
# interactive
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=dtiRFE;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/rfe_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in fa ad rd JHUtracts_ADRDonly; do
    phen=~/data/baseline_prediction/dti_${p}_OD0.95_11052019.csv;
    for s in `cat ../random25.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/dti_${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

# 2019-11-06 12:17:39

Let's plot all our results so far, using this manual approach.

```r
mydir='~/data/baseline_prediction/manual_results'
par(mfrow=c(2, 3))
for (j in c('SX_inatt', 'SX_HI')) {
    for (i in c('Next', 'Last', 'Study')) {
        target = sprintf('%s_group%s', j, i)
        tmp = c()
        for (phen in c('dti_fa', 'dti_ad', 'dti_rd', 'struct_area',
                       'struct_volume', 'struct_thickness')) {
            data = read.csv(sprintf('%s/classification_results_RFE_%s_OD0.95_11052019.csv',
                                    mydir, phen), header=0)
            res_rows = which(grepl(data$V1, pattern=target))
            tmp = rbind(tmp, data.frame(group=phen, val=data[res_rows, 'V3']))
            tmp = rbind(tmp, data.frame(group=sprintf('%s_dumb', phen),
                                        val=data[res_rows, 'V5']))
        }
        mytitle = sprintf('%s', target)
        boxplot(as.formula('val ~ group'), data=tmp, main=mytitle, ylim=c(0,1),
                ylab='ROC', las=2, xlab='')
    }
}
```

![](images/2019-11-06-12-47-23.png)

Those were the results for our first attempt using the manual approach. It uses
these parameters:

```python
    params = {"clf__kernel": ['linear', 'rbf'],
            'clf__C': [.001, .01, .1, 1, 10, 100, 1000],
            'selector__percentile': [5, 10, 15, 20]}
    
    estimators = [('some_variace', VarianceThreshold(threshold=0)),
                  ('unit_variance', StandardScaler()),
                  ('selector', SelectPercentile(f_classif)),
                  ('reduce_dim', PCA()),
                  ('clf', SVC(gamma='scale'))]
```

There seems to be something there, especially in RD. Worth running a few more
seeds, if anything. As usual, there was nothing outstanding, just better than
random.

![](images/2019-11-06-12-52-20.png)

The results above use our FPR framework:

```python
    params = {"clf__kernel": ['linear', 'rbf'],
            'clf__C': [.001, .01, .1, 1, 10, 100, 1000],
            # 'selector__percentile': [5, 10, 15, 20],
            'selector__alpha': np.arange(.01, .11, .02)}
    
    estimators = [('some_variace', VarianceThreshold(threshold=0)),
                  ('unit_variance', StandardScaler()),
                #   ('selector', SelectPercentile(f_classif)),
                    ('selector', SelectFpr(f_classif)),
                #   ('reduce_dim', PCA()),
                  ('clf', SVC(gamma='scale'))]
```

We can go up to .7 with RD, so that's exciting. Certainly worth running more
seeds there. 

Finally, this is what I get using RFECV:

![](images/2019-11-06-13-04-36.png)

This is tricky, because it has many good results for HI_group_study. Worth
running more seeds as well, but I wonder if there could be some sort of
overfitting? Worth running with negative seeds as well just to double check...

Let's check just the phenotypes, no dummy results. Then we can look at outcome
results. We'll also change the limits to better appreciate the differences.

```r
mydir='~/data/baseline_prediction/manual_results'
par(mfrow=c(2, 3))
for (j in c('SX_inatt', 'SX_HI')) {
    for (i in c('Next', 'Last', 'Study')) {
        target = sprintf('%s_group%s', j, i)
        tmp = c()
        for (phen in c('dti_fa', 'dti_ad', 'dti_rd', 'struct_area',
                       'struct_volume', 'struct_thickness')) {
            for (pipe in c('', '_RFE', '_FPR')) {
                data = read.csv(sprintf('%s/classification_results%s_%s_OD0.95_11052019.csv',
                                mydir, pipe, phen), header=0)
                res_rows = which(grepl(data$V1, pattern=target))
                tmp = rbind(tmp, data.frame(group=sprintf('%s_%s', phen, pipe),
                                            val=data[res_rows, 'V3']))
            }
        }
        mytitle = sprintf('%s', target)
        boxplot(as.formula('val ~ group'), data=tmp, main=mytitle, ylim=c(0.2,.9),
                ylab='ROC', las=2, xlab='')
    }
}
```

There wasn't a striking difference, but it seems like RFE is doing best,
followed by FPR. FPR is certainly better for RD, but I think that's the only
one... why dont we run a few more seeds for FPR and RFE, just for groupStudy?
Before that, let's look at outcome:

```r
mydir='~/data/baseline_prediction/manual_results'
target = 'lastPersistent'
tmp = c()
for (phen in c('dti_fa', 'dti_ad', 'dti_rd', 'struct_area',
                'struct_volume', 'struct_thickness')) {
    for (pipe in c('', '_RFE', '_FPR')) {
        data = read.csv(sprintf('%s/classification_results%s_%s_OD0.95_DSM5Outcome_11052019.csv',
                        mydir, pipe, phen), header=0)
        res_rows = which(grepl(data$V1, pattern=target))
        tmp = rbind(tmp, data.frame(group=sprintf('%s_%s (%d)', phen, pipe, length(res_rows)),
                                    val=data[res_rows, 'V3']))
    }
}
mytitle = sprintf('%s', target)
boxplot(as.formula('val ~ group'), data=tmp, main=mytitle, ylim=c(0.2,.9),
        ylab='ROC', las=2, xlab='')
```

![](images/2019-11-06-13-26-58.png)

RFE does better as well, but those distributions look funky. They all have 25
seeds, so I'm not sure what the deal is. Worth running more seeds, though.

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=moreFPR;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/univariate_fpr_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
    phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
    for s in `cat ../random125.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=moreRFE;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/rfe_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
    phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
    for s in `cat ../random125.txt`; do
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $res $s" >> $swarm_file;
            done;
        done;
        phen2=~/data/baseline_prediction/${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res $s" >> $swarm_file;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

And let me run RFE with random data just to see if I get those weird
distributions as well.

```bash
# bw
export OMP_NUM_THREADS=4
cd ~/data/baseline_prediction/manual_swarms

jname=rndRFE;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/rfe_classifier.py;
res=~/data/baseline_prediction/manual_results;
for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
    phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
    for s in `cat ../random25.txt ../random125.txt`; do
        i=Study;
        for j in SX_inatt SX_HI; do
            echo "python3 $code $phen ${j}_group${i} $res -$s" >> $swarm_file;
        done;
        phen2=~/data/baseline_prediction/${p}_OD0.95_DSM5Outcome_11052019.csv;
        echo "python3 $code $phen2 lastPersistent $res -$s" >> $swarm_file;
    done;
done;
# swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
#     --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```


# TODO
* plot results (RFE, FPR, percentile) with more seeds
* see how sex and age and qc variables might help results
* do NV vs ADHD
* do 3 group comparison