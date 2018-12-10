# 2018-12-10 10:51:50

BAsed on my latest thoughts from robust_descriptives.md last Friday, let's
include ADHD_NOS in our analysis, and re-run some of the results.

```bash
job_name=NOSstruct;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls struct_*_11142018_260timeDiff12mo.RData.gz`; do
    for nn in nonew_ ''; do
        for pp in None subjScale; do
            for target in OLS_inatt_slope OLS_HI_slope; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ADHDNOS_${nn}${target} 42 winsorize_$pp" >> $swarm_file;
                for i in {1..250}; do
                    echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/structural.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ADHDNOS_${nn}${target} -${RANDOM} winsorize_$pp" >> $swarm_file;
                done;
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

And let's see if this brings up some more DTI results as well:

```bash
job_name=NOSdti;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls dti_??_voxelwise_n2??_09212018.RData.gz`; do
    for nn in nonew_ ''; do
        for pp in None subjScale; do
            for target in OLS_inatt_slope OLS_HI_slope; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ADHDNOS_${nn}${target} 42 winsorize_$pp" >> $swarm_file;
                for i in {1..250}; do
                    echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ADHDNOS_${nn}${target} -${RANDOM} winsorize_$pp" >> $swarm_file;
                done;
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

And because we had so many results in fMRI, we should definitely play there too:

```bash
job_name=NOSmelodic;
mydir=/data/NCR_SBRB/baseline_prediction/;
swarm_file=swarm.desc_${job_name};
rm -rf $swarm_file;
for f in `/bin/ls melodic_*_IC*_09212018.RData.gz`; do
    for nn in nonew_ ''; do
        for pp in None subjScale; do
            for target in OLS_inatt_slope OLS_HI_slope; do
                echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ADHDNOS_${nn}${target} 42 winsorize_$pp" >> $swarm_file;
                for i in {1..250}; do
                    echo "Rscript --vanilla ~/research_code/baseline_prediction/descriptives/dti.R ${mydir}/${f} ${mydir}/long_clin_11302018.csv ADHDNOS_${nn}${target} -${RANDOM} winsorize_$pp" >> $swarm_file;
                done;
            done;
        done;
    done;
done
grep -v union $swarm_file > ${swarm_file}2;
split -l 3000 ${swarm_file}2 ${job_name}_split;
for f in `/bin/ls ${job_name}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -f $f -g 8 -t 2 --time 5:00:00 --partition norm --logdir trash_desc_${job_name} --job-name ${job_name} -m R,afni 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

## Compiling results

### structural

```bash
myfile=struct_NOSdescriptives.txt
rm $myfile; touch $myfile;
for f in `/bin/ls /data/NCR_SBRB/tmp/struct_*_11142018_260timeDiff12mo/ADHDNOS*_42_?h_ClstTable_e1_a1.0.1D`; do
    if ! grep -q 'rnd' $f; then
        echo $f >> $myfile;
        grep -v \# $f | head -n 5 >> $myfile;
    fi
done
```

```bash
/bin/ls -1 /data/NCR_SBRB/tmp/struct_*_11142018_260timeDiff12mo/ADHDNOS*_42_?h_ClstTable_e1_a1.0.1D | grep -v rnd > result_files.txt;
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
tar -zcvf struct_ADHDNOS_top_rnd_clusters.tar.gz struct_*_11142018_260timeDiff12mo/ADHDNOS*top_rnd_clusters.txt
```

```r
res_fname = '~/tmp/struct_NOSdescriptives.txt'
out_file = '~/tmp/pvals_NOSstruct.txt'
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
    if (file.exists(rnd_fname)) {
        rnd_results = read.table(rnd_fname)[, 3]
        nperms = length(rnd_results)
    } else {
        rnd_results = NA
        nperms = NA
    }
    if (grepl(pattern='lh', line)) {
      cat(sprintf('%s (LH): %s (%d perms)\n', dir_name, root_fname, nperms),
          file=out_file, append=T)
    } else {
      cat(sprintf('%s (RH): %s (%d perms)\n', dir_name, root_fname, nperms),
          file=out_file, append=T)
    }
  } 
  else {
    parsed = strsplit(line, ' +')
    clus_size = as.numeric(parsed[[1]][4])
    pval = sum(rnd_results >= clus_size) / nperms
    cat(sprintf('Cluster size: %.2f, p<%.3f', clus_size, pval),
        file=out_file, append=T)
    if ( !is.na(pval) && pval < .05) {
      cat(' *', file=out_file, append=T)
    }
    if ( !is.na(pval) && pval < .01) {
      cat('*', file=out_file, append=T)
    }
    cat('\n', file=out_file, append=T)
  }
}
```

```bash
sudregp@HG-02070684-DM2:~/tmp$ grep -B 1 "*" pvals_NOSstruct.txt
struct_area_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_HI_slope_winsorize_None (249 perms)
Cluster size: 1015.21, p<0.008 **
--
struct_area_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_HI_slope_winsorize_subjScale (249 perms)
Cluster size: 751.57, p<0.008 **
--
struct_area_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_None (248 perms)
Cluster size: 1098.04, p<0.016 *
Cluster size: 825.97, p<0.040 *
--
struct_thickness_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_None (249 perms)
Cluster size: 306.69, p<0.032 *
--
struct_thickness_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_subjScale (250 perms)
Cluster size: 372.06, p<0.020 *
--
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_HI_slope_winsorize_None (250 perms)
Cluster size: 885.11, p<0.000 **
--
struct_volume_11142018_260timeDiff12mo (LH): ADHDNOS_nonew_OLS_inatt_slope_winsorize_None (250 perms)
Cluster size: 574.77, p<0.000 **
--
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_inatt_slope_winsorize_None (250 perms)
Cluster size: 489.02, p<0.016 *
--
struct_volume_11142018_260timeDiff12mo (LH): ADHDNOS_nonew_OLS_inatt_slope_winsorize_subjScale (248 perms)
Cluster size: 387.28, p<0.012 *
--
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_inatt_slope_winsorize_subjScale (248 perms)
Cluster size: 246.19, p<0.044 *
--
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_OLS_HI_slope_winsorize_None (249 perms)
Cluster size: 604.57, p<0.008 **
--
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_OLS_HI_slope_winsorize_subjScale (240 perms)
Cluster size: 447.18, p<0.037 *
--
struct_volume_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_None (248 perms)
Cluster size: 911.02, p<0.000 **
--
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_OLS_inatt_slope_winsorize_None (248 perms)
Cluster size: 804.45, p<0.000 **
--
struct_volume_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_subjScale (248 perms)
Cluster size: 489.78, p<0.008 **
```

Like before, the none results were much better, but we need to be careful we're
not being fooled by outliers. So, let's filter them a bit:

```bash
sudregp@HG-02070684-DM2:~/tmp$ grep -B 1 "*" pvals_NOSstruct.txt | grep None
struct_area_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_HI_slope_winsorize_None (249 perms)
struct_area_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_None (248 perms)
struct_thickness_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_None (249 perms)
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_HI_slope_winsorize_None (250 perms)
struct_volume_11142018_260timeDiff12mo (LH): ADHDNOS_nonew_OLS_inatt_slope_winsorize_None (250 perms)
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_nonew_OLS_inatt_slope_winsorize_None (250 perms)
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_OLS_HI_slope_winsorize_None (249 perms)
struct_volume_11142018_260timeDiff12mo (LH): ADHDNOS_OLS_inatt_slope_winsorize_None (248 perms)
struct_volume_11142018_260timeDiff12mo (RH): ADHDNOS_OLS_inatt_slope_winsorize_None (248 perms)
```

So it looks like nonew is not making much difference here. We should probably
use the scatterplots and the other imaging modalities to see if we should go
with nonew or the entire ADHD dataset.

# TODO
* monitor calculation of other networks in same_space note