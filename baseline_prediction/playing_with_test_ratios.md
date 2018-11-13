# 2018-11-13 15:33:10

I'm going to check how different test ratios affect the spatialAverage results.
We run the 90/10 we normally do, but go all the way down to 50/50. Let's stick
to a given target and property for now.

```bash
job_name=DTIAD_spatialAvgTestRatiosDL;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.automl_${job_name};
rm -rf $swarm_file;
f="${mydir}/dti_fa_voxelwise_n223_09212018.RData.gz";
mincluster='8';
target=nvVSper;
for train in .5 .6 .7 .8 .9; do
    for i in {1..100}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_spatialAverage_multiDomain_varTest_autoValidation_DL.R $f ${mydir}/long_clin_0918.csv ${target} ${mydir}/tmp/${USER} $RANDOM $mincluster $train" >> $swarm_file;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 -t 16 --time 6:00:00 --partition norm --logdir trash_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:10 2> swarm_wait_${USER};
```