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

Finally, I created a generic function to run all the crappy domains, and also
some of the fMRI transformations. I might need to expand that generic script for
fMRI though, if I want to include movement in the covariates.

So, while we wait for some stuff to run, we can go ahead and grab the top 5
clusters of each result:

```bash
myfile=dti_descriptives.txt
rm $myfile; touch $myfile;
for f in `/bin/ls /data/NCR_SBRB/tmp/dti_??_voxelwise_n2??_09212018/*_42_clusters.txt`; do
    echo $f >> $myfile;
    grep -v \# $f | head -n 5 >> $myfile;
done
```

So, some things that are good to know from these preliminary results. There's no
difference between None and log results. And it looks like there are some good
results, at least based on number of voxels. But we need to see what's actually
significant, and for that we need to wait for the random runs to finish
running...

```bash
/bin/ls -1 /data/NCR_SBRB/tmp/dti_??_voxelwise_n2??_09212018/*_42_clusters.txt > result_files.txt;
for root_file in `cat result_files.txt | sed -e 's/_42_clusters.txt//g'`; do
    collect_name=${root_file}_top_rnd_clusters.txt;
    echo $collect_name;
    if [ -e $collect_name ]; then
        rm $collect_name;
    fi;
    for f in `ls ${root_file}*rnd*clusters.txt`; do
        grep -v \# $f | head -n 1 >> $collect_name;
    done
done
tar -zcvf dti_top_rnd_clusters.tar.gz dti_??_voxelwise_n2??_09212018/*top_rnd_clusters.txt
```

Now we can just write something in R to compute the p-value for each of the
clusters in the descriptives file from above.

# 2018-12-03 10:26:58

```r
res_fname = '~/tmp/dti_descriptives.txt'
res_lines = readLines(res_fname)
for (line in res_lines) {
  # starting new file summary
  if (grepl(pattern='clusters', line)) {
    root_fname = strsplit(line, '/')[[1]]
    dir_name = root_fname[length(root_fname)-1]
    root_fname = strsplit(root_fname[length(root_fname)], '_')[[1]]
    root_fname = paste0(root_fname[1:(length(root_fname)-2)], sep='', collapse='_')
    rnd_fname = sprintf('~/tmp/%s/%s_top_rnd_clusters.txt', dir_name, root_fname)
    rnd_results = read.table(rnd_fname)[, 1]
    nperms = length(rnd_results)
    cat(sprintf('%s: %s (%d perms)\n', dir_name, root_fname, nperms))
  } 
  else {
    parsed = strsplit(line, ' +')
    clus_size = as.numeric(parsed[[1]][2])
    pval = sum(rnd_results >= clus_size) / nperms
    cat(sprintf('Cluster size: %d, p<%.3f', clus_size, pval))
    if (pval < .05) {
      cat(' *')
    }
    if (pval < .01) {
      cat('*')
    }
    cat('\n')
  }
}
```

This spits out a whole lot of result. Unfortunately, the only significant or close to nominal ones were (cropped for better reading):

```
dti_ad_voxelwise_n223_09212018: SX_inatt_baseline_subjScale (499 perms)
Cluster size: 43, p<0.026 *
dti_ad_voxelwise_n223_09212018: SX_inatt_baseline_subjScale-log (250 perms)
Cluster size: 41, p<0.020 *
dti_fa_voxelwise_n272_09212018: SX_HI_baseline_subjScale (499 perms)
Cluster size: 40, p<0.034 *
```

Say that result is good (particularly the AD result, which seems robust), then where is it? And, going with the current idea, does it correlate with outcome as well, or is it mostly related to baseline?

## structural

And we run the same thing for structural runs, except that some file names change, and we have to deal with LH and RH:

```bash
myfile=struct_descriptives.txt
rm $myfile; touch $myfile;
for f in `/bin/ls /data/NCR_SBRB/tmp/struct_*_11142018_260timeDiff12mo/*_42_?h_ClstTable_e1_a1.0.1D`; do
    if ! grep -q 'rnd' $f; then
        echo $f >> $myfile;
        grep -v \# $f | head -n 5 >> $myfile;
    fi
done
```

```bash
/bin/ls -1 /data/NCR_SBRB/tmp/struct_*_11142018_260timeDiff12mo/*_42_?h_ClstTable_e1_a1.0.1D | grep -v rnd > result_files.txt;
sed -i -e 's/_42_lh_ClstTable_e1_a1.0.1D//g' result_files.txt;
sed -i -e 's/_42_rh_ClstTable_e1_a1.0.1D//g' result_files.txt;
for root_file in `cat result_files.txt`; do
    collect_name_lh=${root_file}_lh_top_rnd_clusters.txt;
    collect_name_rh=${root_file}_rh_top_rnd_clusters.txt;
    echo $collect_name_lh;
    echo $collect_name_rh;
    if [ -e $collect_name_lh ]; then
        rm $collect_name_lh $collect_name_rh;
    fi;
    for f in `ls ${root_file}_rnd*lh_ClstTable_e1_a1.0.1D`; do
        grep -v \# $f | head -n 1 >> $collect_name_lh;
    done
    for f in `ls ${root_file}_rnd*rh_ClstTable_e1_a1.0.1D`; do
        grep -v \# $f | head -n 1 >> $collect_name_rh;
    done
done
tar -zcvf struct_top_rnd_clusters.tar.gz struct_*_11142018_260timeDiff12mo/*top_rnd_clusters.txt
```

And finally, in R:

```r
res_fname = '~/tmp/struct_descriptives.txt'
res_lines = readLines(res_fname)
for (line in res_lines) {
  # starting new file summary
  if (grepl(pattern='data', line)) {
    root_fname = strsplit(line, '/')[[1]]
    dir_name = root_fname[length(root_fname)-1]
    root_fname = strsplit(root_fname[length(root_fname)], '_')[[1]]
    root_fname = paste0(root_fname[1:(length(root_fname)-5)], sep='', collapse='_')
    if (grepl(pattern='lh', line)) {
      rnd_fname = sprintf('~/tmp/%s/%s_lh_top_rnd_clusters.txt', dir_name, root_fname)
    } else {
      rnd_fname = sprintf('~/tmp/%s/%s_rh_top_rnd_clusters.txt', dir_name, root_fname)
    }
    rnd_results = read.table(rnd_fname)[, 3]
    nperms = length(rnd_results)
    if (grepl(pattern='lh', line)) {
      cat(sprintf('%s (LH): %s (%d perms)\n', dir_name, root_fname, nperms))
    } else {
      cat(sprintf('%s (RH): %s (%d perms)\n', dir_name, root_fname, nperms))
    }
  } 
  else {
    parsed = strsplit(line, ' +')
    clus_size = as.numeric(parsed[[1]][4])
    pval = sum(rnd_results >= clus_size) / nperms
    cat(sprintf('Cluster size: %.2f, p<%.3f', clus_size, pval))
    if (pval < .05) {
      cat(' *')
    }
    if (pval < .01) {
      cat('*')
    }
    cat('\n')
  }
}
```

We get these results:

```
struct_area_11142018_260timeDiff12mo (LH): nvVSadhd_log (248 perms)
Cluster size: 1303.41, p<0.008 **
struct_area_11142018_260timeDiff12mo (RH): nvVSadhd_log (248 perms)
Cluster size: 1632.22, p<0.012 *
struct_area_11142018_260timeDiff12mo (LH): nvVSadhd_None (250 perms)
Cluster size: 1279.53, p<0.040 *
struct_area_11142018_260timeDiff12mo (RH): nvVSadhd_None (250 perms)
Cluster size: 1727.71, p<0.024 *
struct_area_11142018_260timeDiff12mo (LH): SX_HI_baseline_log (247 perms)
Cluster size: 1939.40, p<0.004 **
Cluster size: 1215.41, p<0.020 *
Cluster size: 1091.86, p<0.036 *
struct_area_11142018_260timeDiff12mo (RH): SX_HI_baseline_log (247 perms)
Cluster size: 2534.08, p<0.000 **
Cluster size: 1458.02, p<0.012 *
Cluster size: 1068.31, p<0.045 *
struct_area_11142018_260timeDiff12mo (LH): SX_HI_baseline_None (248 perms)
Cluster size: 1827.40, p<0.000 **
Cluster size: 1249.68, p<0.012 *
Cluster size: 1074.18, p<0.036 *
struct_area_11142018_260timeDiff12mo (RH): SX_HI_baseline_None (248 perms)
Cluster size: 2754.37, p<0.004 **
Cluster size: 1510.89, p<0.012 *
struct_area_11142018_260timeDiff12mo (LH): SX_inatt_baseline_log (249 perms)
Cluster size: 1058.70, p<0.040 *
struct_area_11142018_260timeDiff12mo (LH): SX_inatt_baseline_None (249 perms)
Cluster size: 1273.39, p<0.032 *
struct_thickness_11142018_260timeDiff12mo (RH): ADHDonly_SX_inatt_baseline_log (250 perms)
Cluster size: 311.79, p<0.004 **
struct_thickness_11142018_260timeDiff12mo (RH): ADHDonly_SX_inatt_baseline_None (249 perms)
Cluster size: 311.79, p<0.016 *
struct_volume_11142018_260timeDiff12mo (LH): ADHDonly_SX_HI_baseline_log (249 perms)
Cluster size: 326.75, p<0.024 *
struct_volume_11142018_260timeDiff12mo (LH): nvVSadhd_log (250 perms)
Cluster size: 552.79, p<0.012 *
struct_volume_11142018_260timeDiff12mo (RH): nvVSadhd_log (250 perms)
Cluster size: 842.26, p<0.000 **
struct_volume_11142018_260timeDiff12mo (LH): nvVSadhd_None (250 perms)
Cluster size: 728.74, p<0.008 **
struct_volume_11142018_260timeDiff12mo (RH): nvVSadhd_None (250 perms)
Cluster size: 861.03, p<0.012 *
struct_volume_11142018_260timeDiff12mo (LH): SX_HI_baseline_log (246 perms)
Cluster size: 506.76, p<0.028 *
struct_volume_11142018_260timeDiff12mo (RH): SX_HI_baseline_log (246 perms)
Cluster size: 844.02, p<0.008 **
struct_volume_11142018_260timeDiff12mo (LH): SX_HI_baseline_None (250 perms)
Cluster size: 506.76, p<0.036 *
struct_volume_11142018_260timeDiff12mo (RH): SX_HI_baseline_None (250 perms)
Cluster size: 937.94, p<0.004 **
```

We have lots of results, but they're mostly showing that baseline and log don't differ much, as we had noticed before. It also looks like the cluster sizes with no transformation are a bit bigger, regardless of cluster significance, so let's go with that.

Like DTI, we still need to plot those results to see if they are in neat regions.

## fmri

It might be faster if I run permutations only for the results with big clusters... just to save some time. I'll do that starting with fMRI.

```bash
job_name=melodic;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls melodic_*_IC*_09212018.RData.gz`; do
    for target in nvVSadhd SX_inatt_baseline SX_HI_baseline \
        ADHDonly_SX_inatt_baseline ADHDonly_SX_HI_baseline; do
        for pp in None subjScale log subjScale-log; do
            echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/melodic.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ${target} 42 $pp" >> $swarm_file;
        done;
    done;
done
swarm -f $swarm_file -g 20 -t 2 --time 5:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni --gres=lscratch:2
```



# TODO

* compute p-values
* crappy domains descriptives
* rsfmri connectivity matrices
* rsfmri icasso