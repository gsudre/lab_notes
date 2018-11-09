# 2018-11-09 15:51:33

This is the first try combining data cross domains. Let's reduce the variables before we do anything, but we should try it differently later. Also, let's combine everything we have, and we can tailor it down later.

```bash
job_name=combineStructDTI_spatialAvgTestDL;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.automl_${job_name};
rm -rf $swarm_file;
f="${mydir}/struct_thickness_09192018_260timeDiff12mo.RData.gz,${mydir}/struct_area_09192018_260timeDiff12mo.RData.gz,${mydir}/struct_volume_09192018_260timeDiff12mo.RData.gz,${mydir}/dti_fa_voxelwise_n223_09212018.RData.gz,${mydir}/dti_ad_voxelwise_n223_09212018.RData.gz,${mydir}/dti_rd_voxelwise_n223_09212018.RData.gz";
mincluster='35,35,35,8,8,8'
for target in nvVSper nvVSrem perVSrem nvVSadhd; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_spatialAverage_multiDomain_test_autoValidation_DL.R $f ${mydir}/long_clin_0918.csv ${target} ${mydir}/models_spatial_across_DL/${USER} $RANDOM $mincluster" >> $swarm_file;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 6:00:00 --partition norm --logdir trash_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:10 2> swarm_wait_${USER};
```

But of course we should look into combining just within domain for comparison:

```bash
job_name=combineStruct_spatialAvgTestDL;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.automl_${job_name};
rm -rf $swarm_file;
f="${mydir}/struct_thickness_09192018_260timeDiff12mo.RData.gz,${mydir}/struct_area_09192018_260timeDiff12mo.RData.gz,${mydir}/struct_volume_09192018_260timeDiff12mo.RData.gz";
mincluster='35,35,35'
for target in nvVSper nvVSrem perVSrem nvVSadhd; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_spatialAverage_multiDomain_test_autoValidation_DL.R $f ${mydir}/long_clin_0918.csv ${target} ${mydir}/models_spatial_across_DL/${USER} $RANDOM $mincluster" >> $swarm_file;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 6:00:00 --partition norm --logdir trash_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:10 2> swarm_wait_${USER};
```

```bash
job_name=combineDTI_spatialAvgTestDL;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.automl_${job_name};
rm -rf $swarm_file;
f="${mydir}/dti_fa_voxelwise_n223_09212018.RData.gz,${mydir}/dti_ad_voxelwise_n223_09212018.RData.gz,${mydir}/dti_rd_voxelwise_n223_09212018.RData.gz";
mincluster='8,8,8'
for target in nvVSper nvVSrem perVSrem nvVSadhd; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_spatialAverage_multiDomain_test_autoValidation_DL.R $f ${mydir}/long_clin_0918.csv ${target} ${mydir}/models_spatial_across_DL/${USER} $RANDOM $mincluster" >> $swarm_file;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 6:00:00 --partition norm --logdir trash_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:10 2> swarm_wait_${USER};
```