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

But it turned out not to be the problem. It was memory, even though jobhist said
it was fine. I'm running them now on norm partition under -t 16 and g10 (jen). Let's
see if that runs. Also, under --quick as Philip. I'm running ADtest under
myself.

Just as a reminder, in the autoValidation function we take the univariate over
90% of the data, which then gets split by autoML into 80-20 for train/split. In
the regular unit_test, the entire dataset is split into 80-10-10, and I used
just the train data to get the univariates.

Because there's a long weekend coming, I'm going to go ahead and queue the
autoValidation for rsfmri, thickness, and dti_ad:

```bash
job_name=rsFMRIautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=ADautoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

```bash
job_name=thicknessAutoframe;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/struct_thickness_09192018_260timeDiff12mo.RData.gz; 
for target in ADHDonly_diag_group2 ADHDonly_OLS_inatt_slope ADHDonly_OLS_HI_slope; do
    for i in {1..150}; do
        echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test/${target} $RANDOM" >> $swarm_file;
    done; 
done;
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
swarm -f $swarm_file -g 40 --partition quick -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10
```

And just before I leave, here's what everyone is running:

```
User     JobId            JobName       Part      St  Reason    Runtime   Walltime    Nodes  CPUs  Memory      Dependency Nodelist
====================================================================================================================================
sudregp  10617195_[109-1  ADtest        quick     PD  Priority      0:00     2:00:00      1  1312   30GB/node                      
sudregp  10623633         sinteractive  interact  R   ---        1:11:14  1-12:00:00      1    32  120GB/node              cn3239  
sudregp  10628584_[0-449  thicknessAut  quick     PD  ---           0:00     3:00:00      1  7200   40GB/node                      
sudregp  10628691_[0-449  rndThickness  quick     PD  ---           0:00     3:00:00      1  7200   40GB/node                      
====================================================================================================================================
cpus running = 32
cpus queued = 15712
jobs running = 1
jobs queued = 941
```

```
[shawp@biowulf research_code]$ sjobs
User   JobId            JobName       Part   St  Reason    Runtime  Walltime  Nodes  CPUs  Memory     Dependency Nodelist
===========================================================================================================================
shawp  10627819_[0-449  rndDTIautofr  quick  PD  Priority     0:00   3:00:00      1  7200  40GB/node                      
shawp  10628593_[0-449  ADautoframe   quick  PD  ---          0:00   3:00:00      1  7200  40GB/node                      
shawp  10628657_[0-449  rndADautofra  quick  PD  ---          0:00   3:00:00      1  7200  40GB/node                      
===========================================================================================================================
cpus running = 0
cpus queued = 21600
jobs running = 0
jobs queued = 1350
```

```
...
frederickja  10626220_99      DTIautoframe  norm   R   ---       16:57   2:00:00      1    16  40GB/node              cn3524  
frederickja  10628602_[0-449  rsFMRIautofr  quick  PD  ---        0:00   3:00:00      1  7200  40GB/node                      
frederickja  10628724_[0-449  rndrsFMRIaut  quick  PD  ---        0:00   3:00:00      1  7200  40GB/node                      
===============================================================================================================================
cpus running = 2400
cpus queued = 14400
jobs running = 150
jobs queued = 900
```

# 2018-10-09 14:50:45

Alright, let's collect the autoframe results we generated over the long weekend.

```bash
echo "target,pheno,var,nfeat,model,metric,val,dummy" > autoframe_summary.csv;
for dir in rsFMRIautoframe DTIautoframe ADautoframe thicknessAutoframe; do
    echo $dir
    for f in `ls trash_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 5`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        if grep -q "Class distribution" $f; then
            dummy=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        else
            dummy=`grep -A 1 "MSE prediction mean" $f | tail -1 | awk '{FS=" "; print $2}'`;
        fi
        echo $target,$phen,$var,$nfeat,$model,$metric,$acc,$dummy >> autoframe_summary.csv;
    done;
done
```

And we can also check the results with randomized data:

```bash
echo "target,pheno,var,nfeat,model,metric,val,dummy" > rndAutoframe_summary.csv;
for dir in rndrsFMRIautoframe rndDTIautoframe rndADautoframe rndThicknessAutoframe; do
    echo $dir
    for f in `ls trash_${dir}/*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 5`;
        target=`head -n 2 $f | tail -1 | awk '{FS=" "; print $8}'`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        if grep -q "Class distribution" $f; then
            dummy=`grep -A 1 "Class distribution" $f | tail -1 | awk '{FS=" "; {for (i=2; i<=NF; i++) printf $i ";"}}'`;
        else
            dummy=`grep -A 1 "MSE prediction mean" $f | tail -1 | awk '{FS=" "; print $2}'`;
        fi
        echo $target,$phen,$var,$nfeat,$model,$metric,$acc,$dummy >> rndAutoframe_summary.csv;
    done;
done
```

It looks like the randomized results still do well. Not sure why though. Is it finding structure in the data for meaningless labels? Maybe... let's plot the other results using the normal dummy parameters, and see what we get.

It looks like the best approach here will be to create an R function that loads
a series of models. For each model, grab the seed ml@parameters$seed, spli the
data, and predict(). We can then use those predictions to compute the exact same
metrics we're using in the random_autoValidation function using caret. It's
probably less costly to load the data once, and split it for each model based on
the seed. The list of models can be gathered by doing grep gpfs on the output
files. 

# 2018-10-12 14:39:52

After my meeting with Philip (see TODO), I made changes to the autovalidate
code. For now, we'll run it just for deep learning, but then we can run other
algorithms. 

```bash
job_name=rsFMRI_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/aparc.a2009s_trimmed_n215_09182018.RData.gz;
for nn in '' nonew_; do
    for target in nvVSper nvVSrem perVSrem; do
        for i in {1..150}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
        done; 
    done;
    for sx in inatt HI total; do
        for target in nvVSimp nvVSnonimp impVSnonimp; do
            for i in {1..150}; do
                echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
            done; 
        done;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait
    while grep -q ERROR swarm_wait; do
        echo "Trying $f"
        swarm -f $f -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait;
        if grep -q ERROR swarm_wait; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Running that on my account. And then we have to run it for AD, ALL_DTI, and struct. We can run the nulls
(cognition, adhd200, etc) later.

```bash
job_name=AD_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz;
for nn in '' nonew_; do
    for target in nvVSper nvVSrem perVSrem; do
        for i in {1..150}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
        done; 
    done;
    for sx in inatt HI total; do
        for target in nvVSimp nvVSnonimp impVSnonimp; do
            for i in {1..150}; do
                echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
            done; 
        done;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait
    while grep -q ERROR swarm_wait; do
        echo "Trying $f"
        swarm -f $f -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait;
        if grep -q ERROR swarm_wait; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
job_name=ALL_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/dti_ALL_voxelwise_n223_09212018.RData.gz;
for nn in '' nonew_; do
    for target in nvVSper nvVSrem perVSrem; do
        for i in {1..150}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
        done; 
    done;
    for sx in inatt HI total; do
        for target in nvVSimp nvVSnonimp impVSnonimp; do
            for i in {1..150}; do
                echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
            done; 
        done;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait
    while grep -q ERROR swarm_wait; do
        echo "Trying $f"
        swarm -f $f -g 40 -t 16 --time 3:00:00 --partition quick --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait;
        if grep -q ERROR swarm_wait; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
job_name=thickness_DL;
swarm_file=swarm.automl_${job_name};
f=/data/NCR_SBRB/baseline_prediction/struct_thickness_09192018_260timeDiff12mo.RData.gz;
for nn in '' nonew_; do
    for target in nvVSper nvVSrem perVSrem; do
        for i in {1..150}; do
            echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
        done; 
    done;
    for sx in inatt HI total; do
        for target in nvVSimp nvVSnonimp impVSnonimp; do
            for i in {1..150}; do
                echo "Rscript --vanilla ~/research_code/automl/uni_test_autoValidation_DL.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv groupOLS_${sx}_slope_${target} /data/NCR_SBRB/baseline_prediction/models_test_DL/${target} $RANDOM" >> $swarm_file;
            done; 
        done;
    done;
done
sed -i -e "s/^/unset http_proxy; /g" $swarm_file;
split -l 1000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait
    while grep -q ERROR swarm_wait; do
        echo "Trying $f"
        swarm -f $f -g 40 -t 16 --time 3:00:00 --logdir trash_${job_name} --job-name ${job_name} -m R --gres=lscratch:10 2> swarm_wait;
        if grep -q ERROR swarm_wait; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Then I'm running AD as philip and ALL as Jen. I'm running thickness as myself on
norm.