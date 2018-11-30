# 2018-11-29 16:03:41

Let's take a turn to run a simply descriptives analysis on these data. The
approach will be to look for regions or variables that are good for DX, inatt
and HI using the entire dataset. Let's employ the traditional ways for multiple
comparisons correction. Then, we take the good results and test them for
associations with outcome, slopes, and latent classes.

Regardless of yay or nay with that analysis, we can look for the last 3
assocations in the entire brain as well, just to see if variables not associated
with baseline status can be changing in outcome.

Finally, we can show our current ML results to test if non-linear relationships
on the data make life better. 

The final step, if any of these tests work on outcome, will be to use it to
predict a test set. That test set could be just the people Wendy is compiling,
or we could also reduce our initial set until right before the results
disappear, and use the remaining ones as testing. Not really an order there. We
could go chronologically, if it makes life simpler. Kinda like mimicing what
we're trying to do now with the phone interviews.

## DTI

Let's start with DTI, which is the most straight-forward one. But for all
neuroimaging analysis, we should probably decide if we're using 3dClustSim, or
going with a permutation analysis... let's do 3dClustSim for now, as we want to
stay open to new results and permutations might end up quite strict especially
in the unbalanced groups.

http://blog.cogneurostats.com/2013/03/08/correlating-brain-and-behavior-in-afni/

http://blog.cogneurostats.com/2016/10/26/performing-brain-behavior-correlations-with-3dttest/

https://afni.nimh.nih.gov/pub/dist/edu/latest/afni_handouts/Clusters_2017.pdf

https://afni.nimh.nih.gov/pub/dist/doc/program_help/3dttest++.html

Note from links above that it looks like 3dtest++ still has issues with
computing significant clusters for covariates (i.e. brain-behavior
correlations)...

Humpf, OK, let's just code it in R, and have it generate the results using
permuted data. This way we know exactly what it's doing... 

Also, might be a good idea to try the data transforms here as well, just to see
what we get.

# 2018-11-30 12:30:12

Now that we have functions to run at least the DTI, let's swarm it in the
cluster.

```bash
job_name=dti;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls dti_??_voxelwise_n2??_09212018.RData.gz`; do
    for target in nvVSadhd SX_inatt_baseline SX_HI_baseline \
        ADHDonly_SX_inatt_baseline ADHDonly_SX_HI_baseline; do
        for pp in None subjScale log subjScale-log; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
            for i in {1..250}; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} -${RANDOM} $pp" >> $swarm_file;
            done;
        done;
    done;
done
split -l 3000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 4 -t 2 --time 30:00 --partition quick --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

The structural analysis follows the same idea as the DTI analysis, except that
we call a different script because it needs to use SurfClust:

```bash
job_name=struct;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls struct_*_11142018_260timeDiff12mo.RData.gz`; do
    for target in nvVSadhd SX_inatt_baseline SX_HI_baseline \
        ADHDonly_SX_inatt_baseline ADHDonly_SX_HI_baseline; do
        for pp in None subjScale log subjScale-log; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
            for i in {1..250}; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} -${RANDOM} $pp" >> $swarm_file;
            done;
        done;
    done;
done
split -l 3000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 4 -t 2 --time 30:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Note that we need lscratch here because read.xls uses it to convert the file to
csv!

For fMRI it's very similar, except that we can play with the MELODIC outputs
here, since we're doing a straight-up neuroimaging clustering analysis. Note
that we'll need bigger machines too, and more time, as the fmri voxelwise data
is much bigger.

```bash
job_name=melodic;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls melodic_*_IC*_09212018.RData.gz`; do
    for target in nvVSadhd SX_inatt_baseline SX_HI_baseline \
        ADHDonly_SX_inatt_baseline ADHDonly_SX_HI_baseline; do
        for pp in None subjScale log subjScale-log; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
            for i in {1..250}; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} -${RANDOM} $pp" >> $swarm_file;
            done;
        done;
    done;
done
split -l 3000 $swarm_file ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 20 -t 2 --time 5:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

(waiting to finish generating the RData files, then need to see how much memory
we'l need.)

Finally, I created a generic function to run all the crappy domains, and also
some of the fMRI transformations. I might need to expand that generic script for
fMRI though, if I want to include movement in the covariates.

# TODO
 * crappy domains descriptives
 * rsfmri connectivity matrices
 * rsfmri icasso









```bash
for i in {1..250}; do
    Rscript --vanilla ~/research_code/automl/generate_random_uni_spatial.R ~/data/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz ~/data/baseline_prediction/long_clin_0918.csv nvVSadhd $RANDOM;
done
```

And change the target to run it other cores as well...

Now we just need to compile the results similar to before:

```bash
cd ~/data/baseline_prediction/neuroimage
for f in `ls ~/data/tmp/nvVSper_*_rnd_clusters.txt`; do
    grep -v \# $f | head -n 1 >> nvVSper_top_clusters.txt
done
```

And, in R:

```r
> x = read.table('~/data/baseline_prediction/neuroimage/nvVSper_top_clusters.txt')[,1]
> summary(x)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   6.00   11.00   15.00   18.02   20.00   75.00 
> sum(x>20)/length(x)
[1] 0.2310757
> sum(x>40)/length(x)
[1] 0.05577689
> sum(x>42)/length(x)
[1] 0.05179283
> sum(x>43)/length(x)
[1] 0.04780876
> sum(x>32)/length(x)
[1] 0.09960159
```

Well, nvVSper with non-shuffled labels has the biggest cluster at size 48, then
one at 43. Even if we go down to .1 we only get these two.

For perVSrem, I'm at:

```r
> x = read.table('~/data/baseline_prediction/neuroimage/perVSrem_top_clusters.txt')[,1]
> summary(x)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    6.0    11.0    15.5    18.2    22.0    62.0 
> sum(x>36)/length(x)
[1] 0.048
> sum(x>30)/length(x)
[1] 0.084
```

Hum... the biggest cluster we get with non-random data is 10?

# 2018-11-14 13:59:02

For nvVSrem I get 9, and for nvVSadhd I get 43. In other words, whenever I used the remission group I only got 10 or less voxels with the actual data. Maybe it's a reflexion of the number of subjects (only 38 rem)?

Let's see if we do any better looking at correlations. For OLS_SX_inatt I get 15, same for HI, and only 13 for total. If we restrict it to ADHDonly, we get 18 for inatt, but 20 for HI and 10 for total. Does RD or FA do better? Sticking to ADHDonly, inatt gets 12 for RD and 14 for FA. HI got 14 for RD and 14 for FA. Finally, total got 13 for RD and 14 for FA. So, let's generate some random numbers. IT looks like AD for ADHDonly is still our best candidate...

```bash
for i in {1..250}; do
    Rscript --vanilla ~/research_code/automl/generate_random_uni_spatial.R ~/data/baseline_prediction/dti_ad_voxelwise_n223_09212018.RData.gz ~/data/baseline_prediction/long_clin_0918.csv ADHDonly_OLS_HI_slope $RANDOM;
done
```

```r
> x = read.table('~/data/baseline_prediction/neuroimage/inatt_top_clusters.txt')[,1]
> summary(x)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
      6      12      15      18      20      71
> sum(x>17)/length(x)
[1] 0.368
> x = read.table('~/data/baseline_prediction/neuroimage/hi_top_clusters.txt')[,1]
> summary(x)
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
   7.00   12.00   16.00   19.47   23.00   77.00 
> sum(x>17)/length(x)
[1] 0.4056225
> sum(x>19)/length(x)
[1] 0.3333333
```

Yeah, this is a no-go. A vanilla neuroiamging analysis won't fly here.
 -->
