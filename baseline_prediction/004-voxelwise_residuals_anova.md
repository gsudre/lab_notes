# 2019-01-25 15:21:06

Another approach I can try is to do the ANOVAs after residualizing the brain.
The variables to use in the residuals can be chosen in a stepwise regression.
Again, let's look at the structural results first because they run faster.

I have no idea how long they'll take to run, so before I run the randomizations
I want to check how the biggest clusters look like. I did 50 voxels in my
desktop in 2 minutes. So, all 5049 in the structural dataset should take about
3.5h.

```bash
job_name=residANOVAstruct;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls struct_*_11142018_260timeDiff12mo.RData.gz`; do
    for pp in None subjScale; do
        for target in OLS_inatt_categ OLS_HI_categ; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural_resids_anova.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
        done;
    done;
done
swarm -f $swarm_file -g 4 -t 2 --time 4:00:00 --partition quick --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2
```

For DTI it's taking less than 1min for 50 voxels, and we have 12K. So, this
should be about 4h as well. But that didn't work, so I had to increase it to 8h
just to be safe.

```bash
job_name=residANOVAdti;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls dti_??_voxelwise_n2??_09212018.RData.gz`; do
    for pp in None subjScale; do
        for target in OLS_inatt_categ OLS_HI_categ; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti_resids_anova.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
        done;
    done;
done
swarm -f $swarm_file -g 4 -t 2 --time 8:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2
```

In rsFMRI I have 44211 voxels in the intersection mask. As I know that took a
long time to begin with, let's not play around.

```bash
job_name=residANOVArsfmri;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls melodic_fancy_IC*12142018.RData.gz melodic_inter_IC*12142018.RData.gz`; do
    for pp in None subjScale; do
        for target in OLS_inatt_categ OLS_HI_categ; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/melodic_resids_anova.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
        done;
    done;
done
swarm -f $swarm_file -g 4 -t 2 --time 30:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2
```

Just because it's right before the weekend, I should leave the permutations running:

```bash
job_name=residANOVAstructRnd;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls struct_*_11142018_260timeDiff12mo.RData.gz`; do
    for pp in None subjScale; do
        for target in OLS_inatt_categ OLS_HI_categ; do
            for i in {1..250}; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural_resids_anova.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} -${RANDOM} $pp" >> $swarm_file;
            done;
        done;
    done;
done
split -l 300 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 4 -t 2 --time 8:00:00 --partition quick --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
job_name=residANOVAdtiRnd;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls dti_??_voxelwise_n2??_09212018.RData.gz`; do
    for pp in None subjScale; do
        for target in OLS_inatt_categ OLS_HI_categ; do
            for i in {1..250}; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti_resids_anova.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} -${RANDOM} $pp" >> $swarm_file;
            done;
        done;
    done;
done
split -l 300 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 4 -t 2 --time 8:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Running DTI RND as Philip.

```bash
job_name=residANOVArsfmriRnd;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls melodic_fancy_IC*12142018.RData.gz melodic_inter_IC*12142018.RData.gz`; do
    for pp in None subjScale; do
        for target in OLS_inatt_categ OLS_HI_categ; do
            for i in {1..250}; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/melodic_resids_anova.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} -${RANDOM} $pp" >> $swarm_file;
            done;
        done;
    done;
done
split -l 300 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 4 -t 2 --time 30:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```



PHILIP SUGGESTED I SHOULD TRY NOT USING THE NVS AS THE REFERENCE GROUP. JUST FOR PLOTTING! USE ONE OF THE EXTREME GROUPS TO MAKE RISK RATIO BIGGER

THEN, MAYBE DOING QUADRATIC AND CUBIC FITS WITH ALL 4 GROUPS.
