# 2018-10-22 17:12:04

As I wait for the rndTFilter results (see dummy note), let's re-run some of the
raw variables, but now using the pairwise comparisons. 

```bash
job_name=adRaw;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem nvVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/raw_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_raw/${USER} $RANDOM" >> $swarm_file;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 32 --time 1-0:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

```bash
job_name=rsfmriRaw;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz;
rm -rf $swarm_file;
for target in nvVSper perVSrem nvVSrem; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/raw_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_raw/${USER} $RANDOM" >> $swarm_file;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 60 -t 32 --time 1-0:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10;
```

So, I left rsfmriRaw running with Jen, who is also running rnd_trainRaw for AD.
Philip is running adRaw, and I'm still running rnd_trainRF.