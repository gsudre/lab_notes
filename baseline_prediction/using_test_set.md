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

Then, I realized it was not writing any logs, so I might as well just swarm it.

Now it's time to collect the results:

```bash
echo "target,pheno,var,nfeat,model,metric,val,dummy" > autoTest_summary.csv;
old_phen=''
for dir in rsFMRItest structTest DTItest; do
    echo $dir
    for f in `ls trash_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 5`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        if [[ $var == 'raw' ]] && [[ $phen == 'dti_ALL' ]]; then
            nfeat=36066;
        elif [[ $var == 'raw' ]] && [[ $phen = *'dti_'* ]]; then
            nfeat=12022;
        elif grep -q "Running model on" $f; then  # newer versions of the script
            nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        else # older versions were less verbose
            nfeat=`grep -e "26. *[1-9]" $f | grep "\[1\]" | tail -1 | awk '{ print $3 }'`;
        fi
        if grep -q "Class distribution" $f; then
            dummy=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        else
            dummy=`grep -A 1 "MSE prediction mean" $f | tail -1 | awk '{FS=" "; print $2}'`;
        fi
        echo $target,$phen,$var,$nfeat,$model,$metric,$acc,$dummy >> autoTest_summary.csv;
    done;
done
```

The results using dti_ALL weren't as good as I had expected. 115 out of the 150
ran, but we have a mean of .83, with a reasonable spread. Maybe we need to run
AD? I could also set the validation frame to auto, which would use it also in
the univariate selection?

```bash
job_name=ADtest;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 18 --partition quick -t 24 --time 2:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=DTIautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ALL_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 18 --partition quick -t 24 --time 2:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

I should also find a better way to estimate the dummy results. Just class
majority isn't great...

# 2018-10-05 10:33:55

Let's put in the results of the two tests above... luckily our average goes up:

```bash
echo "target,pheno,var,nfeat,model,metric,val,dummy" > autoTest_summary.csv;
old_phen=''
for dir in ADtest DTIautoframe structTest; do
    echo $dir
    for f in `ls trash_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 5`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        if [[ $var == 'raw' ]] && [[ $phen == 'dti_ALL' ]]; then
            nfeat=36066;
        elif [[ $var == 'raw' ]] && [[ $phen = *'dti_'* ]]; then
            nfeat=12022;
        elif grep -q "Running model on" $f; then  # newer versions of the script
            nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        else # older versions were less verbose
            nfeat=`grep -e "26. *[1-9]" $f | grep "\[1\]" | tail -1 | awk '{ print $3 }'`;
        fi
        if grep -q "Class distribution" $f; then
            dummy=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        else
            dummy=`grep -A 1 "MSE prediction mean" $f | tail -1 | awk '{FS=" "; print $2}'`;
        fi
        echo $target,$phen,$var,$nfeat,$model,$metric,$acc,$dummy >> autoTest_summary.csv;
    done;
done
```

Unfortunately lots of them died again... let's rerun at least the diag group
ones. We should then reconsider if running AD and volume is better.

My guess (but without much evidence) is that running these nodes with smaller
number of cores makes more than one job enter the same node, and that generates
issues with the Java server. I'm re-running everything with t32, but it might
take a while to get the resources.

Wolfgang got back to me and said that's right, so I changed the initial port for
the server to be random. This way we can have more than one job per node.





