# 2018-11-15 13:20:34

Running some tests to see if data transformations make things look better. I coded a within-subject normalization, and also a within-variable normalization. Let's see if it makes any difference for our best resul in each of the 2 brain modalities. Maybe try different algorithms too. We can try it with fMRI later.

```bash
job_name=dataTransforms_rawCV;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.automl_${job_name};
rm -rf $swarm_file;
for f in struct_volume_11142018_260timeDiff12mo.RData.gz dti_ad_voxelwise_n223_09212018.RData.gz; do
    for target in nvVSper perVSrem; do
        for pp in subjScale dataScale subjScale,dataScale None; do
            for algo in DeepLearning DRF GBM GLM; do
                for i in {1..100}; do
                    myseed=$RANDOM;
                    echo "Rscript --vanilla ~/research_code/automl/raw_multiDomain_autoValidation_oneAlgo.R ${mydir}/$f ${mydir}/long_clin_0918.csv ${target} ${mydir}/models_raw_dataTransforms/${USER} $myseed $algo $pp" >> $swarm_file;
                    echo "Rscript --vanilla ~/research_code/automl/raw_multiDomain_autoValidation_oneAlgo.R ${mydir}/$f ${mydir}/long_clin_0918.csv ${target} ${mydir}/models_raw_dataTransforms/${USER} -$myseed $algo $pp" >> $swarm_file;
                done;
            done;
        done;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```