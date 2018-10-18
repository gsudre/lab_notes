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
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
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
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
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
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
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
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
    done; 
done;
for sx in inatt HI total; do
    for target in nvVSimp nvVSnonimp impVSnonimp; do
        for i in {1..100}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${USER} $RANDOM" >> $swarm_file;
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