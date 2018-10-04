# 2018-10-04 10:58:57

After chatting with Philip, let's run the best of each validation one-seed test,
looking back at the results in compiling_results from 10/03.

Let's start with DTI, keeping in mind we need to run at least 150 different
seeds to get a somewhat decent distribution:

```bash
job_name=DTItest;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ALL_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 32 --time 4:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

As this seems to be working, let's also run it for struct and rsFMRI:

```bash
job_name=structTest;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/struct_thickness_09192018_260timeDiff12mo.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 32 --time 4:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=rsFMRItest;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 32 --time 4:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

Except that I'm running the struct tests in interactive nodes just to use up as
much processing as I can. I split it in staa..stad. I also reduced the
requirements ot t 24 and g 18 to be able to access a different pool in the quick
partition. Let's hope it doesn't do memory overflow.

Then, I started running struct under Jen's account, but with the old memory
setting because it's structural afterall.

Then, I realized it was not writing any logs, os I might as well just swarm it.
