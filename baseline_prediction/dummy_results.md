# 2018-10-05 13:06:52

I've been struggling with what would be the best way to show results under a
random scenario. For classification, we could run the exact same algorithm after
data randomization, then show results using simply the majority class, and then
finally using the class probabilities for voting.

In the prediction scenario, we could use random data (as above), but also the
mean and the median.

Then, the idea is that we'd have a series of plots, and in each plot would have
a series of boxplots (one per phenotype), and horizontal bars showing the 95th
percentile of each of the random approaches above. There would be a plot for
accuracy, one for AUC, F1-ratio, specificity, and sensitivity. Or maybe get rid
of F1 just to make it symmetric, or add precision to the plots (recall = sensitivity).

For prediction, we can play with RMS and R2.

---

I then realized that we can simply use the probabilities based on the training
set, and not have to do majority and then probabilities. That takes care of
classification (just need to run the random stuff).

So, the code in automl/random_autoValidation.R is working. We just need to run
it a whole bunch of times to get some distributions. Note, a few things:

* Random AUC should be around .5 anyways, so there shouldn't be much variation there.
* Most of the variation will be in the classification accuracy.

---

So, let's start running some of the random data. I don't have results(dti_all)
for the autoValidation script yet, but I have a feeling this will be the
ultimate pipeline for test data, so let's run it using randomized data as well:

```bash
job_name=rndDTIautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ALL_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/rnd_models/${target} -$RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

And because there's a long weekend coming, I'm gonna go ahead and queue the random
stuff for thickness, AD and rsFMRI as well:

```bash
job_name=rndDTIautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ALL_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/rnd_models/${target} -$RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=rndADautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/rnd_models/${target} -$RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=rndThicknessAutoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/struct_thickness_09192018_260timeDiff12mo.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/rnd_models/${target} -$RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=rndrsFMRIautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/rnd_models/${target} -$RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

# 2018-10-18 11:40:42

I changed the uni_test_autoValidation_DL.R code to create a random data matrix
when the seed is negative. Let's run 100 of them just to get a flavor, for each
of the datasets, and see what we get. Also, just to save some time, I won't run
the nonew parameters, but we can always do it later.

```bash
job_name=rnd_rsFMRI_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper nvVSrem perVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
        done; 
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER};
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 60 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
job_name=rnd_thickness_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/struct_thickness_09192018_260timeDiff12mo.RData.gz;
rm -rf $swarm_file;
for target in nvVSper nvVSrem perVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
        done; 
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER};
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 60 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
job_name=rnd_dtiAD_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper nvVSrem perVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
        done; 
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER};
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 60 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
job_name=rnd_dtiALL_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ALL_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper nvVSrem perVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
        done; 
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER};
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 60 -t 16 --time 3:00:00 --partition norm --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

# 2018-10-22 09:20:37

Let's grab some of the random results:

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,spec,sens,prec,ratio" > rndAutoframeDL_summary.csv;
for dir in dtiAD_DL dtiALL_DL rsFMRI_DL thickness_DL; do
    echo $dir
    for f in `ls trash_rnd_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 5`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        spec=`grep -A 8 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        sens=`grep -A 7 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        prec=`grep -A 6 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`
        echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$spec,$sens,$prec,$ratio >> rndAutoframeDL_summary.csv;
    done;
done
```

Well, it turns out that the model is still classifying well, even with the
random data. The only explanation I have here is that by pumping in random data
and doing the univariate selection, we're still selecting features that will
help, even thoiugh they come from random data. In fact, the numbers of features
we're selecting makes sense based on the total number of variables and how many
we'd get by random chance, so the algorithm is working in that sense. But the
trick is that of those features selected, some of them are still going to be
good features by chance, even if there should be not really any relationship
between train and test set. For example, in the DTI dataset, we have abut 12.1K
features, and using random data we're selecting about 600 in univariate filter
(about 5%). Of those, about 30 would still be good features just by chance,
which could be enough to drive the algorithm.

So, I think we have two options. One would be to change the threshold of the
univariate filter to not allow for any features to go through just by chance.
That would need to be adaptive to the number of features going in though, which
gets tricky as this is all based on expectations. The other option is to pump
random data after the univariate selection. However, I don't quite think this
would work as it's the same thing as we're doing now, right?

Let's play with a .01 threshold, as that would give us way less features, and in
the AD case almost no chance of getting good features at random (121 then 1.2
features left), but we should also play with the size of the test set, as the
more random data we sprinkle into the test set, the less likely those random
features will be to remain random.

```bash
job_name=rnd_uni01;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni01_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

For the size of train set, I'll split them up so I can run in different
accounts:

```bash
job_name=rnd_train1;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
target=nvVSper;
for i in {1..100}; do
    myseed=-$RANDOM;
    for train in .5 .6 .7 .8 .9; do
        echo "Rscript --vanilla ~/research_code/automl/uni_varTest_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $train" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

```bash
job_name=rnd_train2;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
target=perVSrem;
for i in {1..100}; do
    myseed=-$RANDOM;
    for train in .5 .6 .7 .8 .9; do
        echo "Rscript --vanilla ~/research_code/automl/uni_varTest_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $train" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

Another option would be to select the top X features foloowing an univariate
filter, instead of doing it based on p-value. Maybe worth a try?

Let's frist collect the uni01 results:

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,spec,sens,prec,ratio" > rndUni01DL_summary.csv;
dir=uni01;
for f in `ls trash_rnd_${dir}/*o`; do
    phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 5`;
    target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
    seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
    var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
    model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
    auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
    nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
    ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
    f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    spec=`grep -A 8 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    sens=`grep -A 7 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    prec=`grep -A 6 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`
    echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$spec,$sens,$prec,$ratio >> rndUni01DL_summary.csv;
done
```

Another option I hadn't considered before is to use Lasso to get the variables,
instead of the univariate filter. It might be worth a try, as we wouldn't be
putting as much weight in the univariate selection. Along the same lines of
thought, I coudl potentially just try a no-CV model in the entire training data
and select the top X variables in the model. Then, those would be the variables
we'd look at.

I don't think either idea would make the algorithms less susceptible to needing
a bigger test set, but they might improve things a bit overall.

Also, I found an error in grabbing sens and other metrics from the output.
Basically, that was the maximum over the entire range, which is not what we
want. We want that value in the actual predictions, or get the actual ROC curve
(not the maximum on each axis). 

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,spec,sens,prec,ratio" > rndTestSize_summary.csv;
for dir in train1 train2; do
    echo $dir;
    for f in `ls trash_rnd_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $11}'`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$spec,$sens,$prec,$ratio >> rndTestSize_summary.csv;
    done;
done
```

Yeah, as I expected, the closer to .5 split I get the closer to random results I
get, whichis what I wanted to see. However, that also means leaving on the floor
a whole bunch of our data. In fact, even with real data we might get a lot os
spurious voxels by doing that approach. So, should we use .01 to be even more
selective?

```bash
job_name=rnd_trainp01p5;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
target=nvVSper;
for i in {1..100}; do
    myseed=-$RANDOM;
    train=.5;
    echo "Rscript --vanilla ~/research_code/automl/uni01_varTest_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $train" >> $swarm_file; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

Also, the more I think o fit, if we're not using the stacked ensembles then thee
is no reason do do CV within the training data. Let's remove it for now.

So, we could go back to raw data and not worry about any of this, or try other
methods (e.g. p<.01, or DRF for feature selection) to see if we are less
succeptible to these split artifacts but still maintain a reasonable level of
results. Also, keep in mind that some algorithms might be better than others in
handling these variations, so it's not a bad idea to consider going back to
full-on AutoML for now.

```bash
job_name=rnd_trainp01;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
target=nvVSper;
for i in {1..100}; do
    myseed=-$RANDOM;
    for train in .5 .6 .7 .8 .9; do
        echo "Rscript --vanilla ~/research_code/automl/uni01_varTest_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $train" >> $swarm_file;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

```bash
job_name=rnd_trainRF;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
target=nvVSper;
for i in {1..100}; do
    myseed=-$RANDOM;
    for train in .5 .6 .7 .8 .9; do
        echo "Rscript --vanilla ~/research_code/automl/rfFilter_varTest_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $train" >> $swarm_file;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

Compiling p01 filter results, using only .5 split and DL:

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,ratio" > rndP01P5DL_summary.csv;
dir=trainp01p5;
echo $dir;
for f in `ls trash_rnd_${dir}/*o`; do
    phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $11}'`;
    target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
    seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
    var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
    model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
    auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
    nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
    ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
    f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$ratio >> rndP01P5DL_summary.csv;
done
```

Even using a p01 filter and .5 ratio DL gives us results as good as .7, on an
almost even split class ratio. That's rough... there's no point in looking at
other algorithms here, because we already know that using univariate filter is
giving us artifacts. Let's look at raw data, as that should give us
the expected results, and wait on the RF results in case having a multivariate
selection might give us better results.

```bash
job_name=rnd_trainRaw;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
target=nvVSper;
for i in {1..100}; do
    myseed=-$RANDOM;
    echo "Rscript --vanilla ~/research_code/automl/raw_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $train" >> $swarm_file;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 32 --time 1-0:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

As a side note, I've tried Lasso for feature selection (lambda=1), but I keep
geeting results with all coefficients 0. If I try the automl model, I can't
extract coefficients...

# 2018-10-23 09:23:17

Let's se eif the results using RF for filtering are more encouraging:

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,ratio" > rndRF_summary.csv;
dir=trainRF;
echo $dir;
for f in `ls trash_rnd_${dir}/*o`; do
    phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $11}'`;
    target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
    seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
    var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
    model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
    auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
    nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
    ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
    f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$ratio >> rndRF_summary.csv;
done
```

Yep, we're still seeing inflated results with random data... do at least the raw
results look as expected?

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,ratio" > rndRaw_summary.csv;
dir=trainRaw;
echo $dir;
for f in `ls trash_rnd_${dir}/*o`; do
    phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $11}'`;
    target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
    seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
    var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
    model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
    auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
    nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
    ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
    f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
    echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$ratio >> rndRaw_summary.csv;
done
```

Let's play a bit with dimensionality reduction.

```bash
job_name=dimRedTests;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem; do
    for i in {1..20}; do
        myseed=$RANDOM;
        for dr in PCA ICA AutoEncoder; do
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $dr" >> $swarm_file;
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$myseed $dr" >> $swarm_file;
        done;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 32 --time 1-0:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

And a few variations of PCA:

```bash
job_name=PCATests;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem; do
    for i in {1..20}; do
        myseed=$RANDOM;
        for dr in PCA PCA-elbow PCA-kaiser; do
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $dr" >> $swarm_file;
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$myseed $dr" >> $swarm_file;
        done;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 32 --time 1-0:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

# 2018-10-24 10:55:27

Now, let's see if there is any difference among those:

```bash
echo "target,pheno,var,seed,nfeat,model,auc,f1,acc,ratio" > dimRed_summary.csv;
for dir in PCATests dimRedTests; do
    echo $dir;
    for f in `ls trash_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $11}'`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        seed=`head -n 2 $f | tail -1 | awk '{FS=" "; print $10}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        auc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        ratio=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        f1=`grep -A 2 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        acc=`grep -A 5 "Maximum Metrics:" $f | tail -1 | awk '{FS=" "; print $5}'`;
        echo $target,$phen,$var,$seed,$nfeat,$model,$auc,$f1,$acc,$ratio >> dimRed_summary.csv;
    done;
done
```

So, it looks like for PCA it's a tie between elbow and kaiser, btu I like kaiser more because elbow seems to fail for random data. If we make the scree plots, many times it doesn't have a discernible elbow, while we can always find a kaiser limit as it only depends on the mean. 

![](2018-10-24-11-18-38.png)

However, it doesn't look like these PCA results are better than just using the raw features, right? (for ICA and AE, I'm still playing with the best way to encode the data)

Let's run a few more iterations of PCA to have a better idea.

```bash
job_name=PCATestsDTI;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem; do
    for i in {1..80}; do
        myseed=$RANDOM;
        for dr in PCA PCA-elbow PCA-kaiser; do
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $dr" >> $swarm_file;
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$myseed $dr" >> $swarm_file;
        done;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 32 --time 1-0:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

```bash
job_name=PCATestsFMRI;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem; do
    for i in {1..100}; do
        myseed=$RANDOM;
        for dr in PCA PCA-elbow PCA-kaiser; do
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $myseed $dr" >> $swarm_file;
            echo "Rscript --vanilla ~/research_code/automl/dimRed_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} -$myseed $dr" >> $swarm_file;
        done;
    done;
done;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 60 -t 32 --time 1-0:00:00 --partition norm --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```