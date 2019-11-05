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

# TODO
* plot results (RFE, FPR, percentile)