# 2018-09-12 14:46:38

Here's how we currently grab the results:

```bash
echo "target,pheno,var,nfeat,model,metric,val" > auto_summary.csv;
for y in diag_group2 random_HI_slope random_total_slope random_inatt_slope \
    OLS_inatt_slope OLS_HI_slope OLS_total_slope \
    group_HI3 group_total3 group_INATT3; do
    for f in `grep -l \"${y} *o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 4 | cut -d"_" -f 1,2 -`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        nfeat=`grep -e "26. " $f | tail -1 | awk '{ print $3 }'`;
        echo $y,$phen,$var,$nfeat,$model,$metric,$acc >> auto_summary.csv;
    done;
done
```

It looks like uni01 is not doing much. For the ones that it actually run (sometimes there were no variables p < .01)), the results were not better than when using uni. That could indicate overfitting, but also there could be features that interact well and hacving so few features (As in uni01) means shooting ourselves in the foot.

# 2018-09-14 12:34:50

I updated the code to use all the improvements/enhancements we have in the limited version, except that it's not limited to algorithms or time. In fact, I think the filtering to the structural data puts the number of variables in a manageable number, so we don't need to run the limited code anymore. And I have alreayd noticed that the raw code breaks because it runs out of memory... let's run it anyways, and let it die.

```bash
rm swarm.automl
for var in raw pca uni pca_uni uni_pca uni01; do
    for m in area thickness volume; do
        for sx in inatt HI total; do
            for sl in OLS random; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/struct_${m}_08312018_260timeDiff12mo.RData.gz ~/data/baseline_prediction/struct_gf_09052018_260timeDiff12mo.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
            done;
        done;
        for sx in INATT HI total; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/struct_${m}_08312018_260timeDiff12mo.RData.gz ~/data/baseline_prediction/struct_gf_09052018_260timeDiff12mo.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
        done;
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/struct_${m}_08312018_260timeDiff12mo.RData.gz ~/data/baseline_prediction/struct_gf_09052018_260timeDiff12mo.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl;
swarm -f swarm.automl -g 40 -t 32 --time 1-00:00:00 --logdir trash_bin --job-name struc -m R --gres=lscratch:10
```

We can also run the snp code. I also don't think there's much value in running the PCA code there, so it's just the raw code.

```bash
rm swarm.automl
var=raw;
for m in {4..9}; do
    for sx in inatt HI total; do
        for sl in OLS random; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_snps1e0${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
        done;
    done;
    for sx in INATT HI total; do
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_snps1e0${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
    done;
    echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_snps1e0${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
done;
m='prs'
for sx in inatt HI total; do
    for sl in OLS random; do
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_${m}_09142018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
    done;
done;
for sx in INATT HI total; do
    echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_${m}_09142018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
done;
echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_${m}_09142018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl;
swarm -f swarm.automl -g 40 -t 32 --time 1-00:00:00 --logdir trash_bin --job-name geno -m R --gres=lscratch:10
```

And since we're running all for comparison, let's attach neuropsych as well:

```bash
for var in raw pca uni pca_uni uni_pca uni01; do
    for m in all cpt wiscraw wiscstd wj iq; do
        for sx in inatt HI total; do
            for sl in OLS random; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/cog_${m}_09142018.RData.gz ~/data/baseline_prediction/cog_gf_09142018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
            done;
        done;
        for sx in INATT HI total; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/cog_${m}_09142018.RData.gz ~/data/baseline_prediction/cog_gf_09142018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
        done;
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/cog_${m}_09142018.RData.gz ~/data/baseline_prediction/cog_gf_09142018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
        done;
    done;
done;
```

# 2018-09-18 10:03:19

Then we collect all results again, but with a small change to the script so we
can correctly capture the number of features used:

```bash
echo "target,pheno,var,nfeat,model,metric,val" > auto_summary.csv;
for y in diag_group2 random_HI_slope random_total_slope random_inatt_slope \
    OLS_inatt_slope OLS_HI_slope OLS_total_slope \
    group_HI3 group_total3 group_INATT3; do
    for f in `grep -l \"${y} *o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 4 | cut -d"_" -f 1,2 -`;
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
        echo $y,$phen,$var,$nfeat,$model,$metric,$acc >> auto_summary.csv;
    done;
done
```

So I saved auto_summary_09182018.csv in ~/Dcouments/baseline_prediction.

# 2018-09-19 10:22:19

Now we run the resting state variables:

```bash
for var in raw pca uni pca_uni uni_pca uni01; do
    for m in aparc aparc.a2009s aparc_trimmed aparc.a2009s_trimmed; do
        for sx in inatt HI total; do
            for sl in OLS random; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/${m}_n215_09182018.RData.gz ~/data/baseline_prediction/rsfmri_gf_09182018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
            done;
        done;
        for sx in INATT HI total; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/${m}_n215_09182018.RData.gz ~/data/baseline_prediction/rsfmri_gf_09182018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
        done;
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/${m}_n215_09182018.RData.gz ~/data/baseline_prediction/rsfmri_gf_09182018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl;
cp swarm.automl swarm.automl_rsfmri;
swarm -f swarm.automl -g 40 -t 32 --time 4:00:00 --partition quick --logdir trash_bin --job-name rsfmri -m R --gres=lscratch:10
```

## Subgroup analysis

I then added some prep code to the functions to handle the different subgroups
we're analyzing. But before I could do that, in order to use a single gf file
across modalities, I had ot make a change to our older datasets, so that they
included the MRN. The merge takes care of only including the MRNs for which we
have that specific dataset...

```r
> library(gdata)
> clin = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'dti')
> load('~/data/baseline_prediction/dti_ad_voxelwise_n263_08152018.RData.gz')
> clin = clin[, c(1,4)]
> colnames(clin) = c('MRN', 'mask.id')
> data = merge(clin, data, by='mask.id')
> dim(data)
[1]   263 12024
> save(data, file='~/data/baseline_prediction/dti_ad_voxelwise_n263_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/dti_rd_voxelwise_n263_08152018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/dti_rd_voxelwise_n263_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/dti_fa_voxelwise_n263_08152018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/dti_fa_voxelwise_n263_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/dti_ALL_voxelwise_n263_08152018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/dti_ALL_voxelwise_n263_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/dti_tracts_voxelwise_n263_08152018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/dti_tracts_voxelwise_n263_09192018.RData.gz', compress=T)

> clin = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'mprage')
> clin = clin[, c(1,4)]
> colnames(clin) = c('MRN', 'mask.id')
> load('~/data/baseline_prediction/struct_rois_09062018_260timeDiff12mo.RData.gz')
> data = merge(clin, data, by='mask.id')
> dim(data)
[1] 260 273
> save(data, file='~/data/baseline_prediction/struct_rois_09192018_260timeDiff12mo.RData.gz', compress=T)
> load('~/data/baseline_prediction/struct_area_08312018_260timeDiff12mo.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/struct_area_09192018_260timeDiff12mo.RData.gz', compress=T)
> load('~/data/baseline_prediction/struct_thickness_08312018_260timeDiff12mo.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/struct_thickness_09192018_260timeDiff12mo.RData.gz', compress=T)
> load('~/data/baseline_prediction/struct_volume_08312018_260timeDiff12mo.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/struct_volume_09192018_260timeDiff12mo.RData.gz', compress=T)

> clin = read.csv('/Volumes/Shaw/prs/FINAL_FILES_08242018/REGRESSION/PRS2017_geno3_1KG_noFlip_genop05MAFbtp01rsbtp9.csv')
> clin = clin[,1:2]
> load('~/data/baseline_prediction/geno3_prs_09142018.RData.gz')
> colnames(clin) = c('mask.id', 'MRN')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_prs_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/geno3_snps1e04_09132018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_snps1e04_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/geno3_snps1e05_09132018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_snps1e05_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/geno3_snps1e06_09132018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_snps1e06_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/geno3_snps1e07_09132018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_snps1e07_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/geno3_snps1e08_09132018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_snps1e08_09192018.RData.gz', compress=T)
> load('~/data/baseline_prediction/geno3_snps1e09_09132018.RData.gz')
> data = merge(clin, data, by='mask.id')
> save(data, file='~/data/baseline_prediction/geno3_snps1e09_09192018.RData.gz', compress=T)

```

Now, we need to construct a super swarm file. I'll also called the job
differently so we can compile results just within these new log files:

```bash
for var in raw pca uni pca_uni uni_pca uni01; do
    for f in `/bin/ls -1 ~/data/baseline_prediction/*0919*gz ~/data/baseline_prediction/cog*gz ~/data/baseline_prediction/aparc*gz`;do
        for nn in nonew_ ''; do
            for g in ADHDonly_ ''; do
                for t in diag_group2 OLS_inatt_slope OLS_HI_slope OLS_total_slope; do
                    echo "Rscript --vanilla ~/research_code/automl/${var}.R $f ~/data/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} ~/tmp/${nn}${g}${t}" >> swarm.automl_subgroup;
                done; 
            done;
            for g in ADHDNOS_ ADHDNOS_group; do
                for t in OLS_inatt_slope OLS_HI_slope OLS_total_slope; do
                    echo "Rscript --vanilla ~/research_code/automl/${var}.R $f ~/data/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} ~/tmp/${nn}${g}${t}" >> swarm.automl_subgroup;
                done; 
            done;
        done;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl_subgroup;
swarm -f swarm.automl_subgroup -g 40 -t 32 --time 4:00:00 --logdir trash_bin --job-name subgroups -m R --gres=lscratch:10
```

# 2018-09-20 09:18:45

SLURM was bundling my swarms, which screwed up my output files. So, I had to
split the runs:

```bash
split -l 500 swarm.automl_subgroup 
for f in `/bin/ls -1 xa?`; do swarm -f $f -g 40 -t 32 --time 4:00:00 --partition quick --logdir trash_bin --job-name subgroups2 -m R --gres=lscratch:10; done
```

While we wait on results, let's collect the rsfmri results to add to the
previous summary:

```bash
echo "target,pheno,var,nfeat,model,metric,val" > auto_summary.csv;
for y in diag_group2 random_HI_slope random_total_slope random_inatt_slope \
    OLS_inatt_slope OLS_HI_slope OLS_total_slope \
    group_HI3 group_total3 group_INATT3; do
    for f in `grep -l \"${y} rsfmri*o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 4 | cut -d"_" -f 1,2 -`;
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
        echo $y,$phen,$var,$nfeat,$model,$metric,$acc >> auto_summary.csv;
    done;
done
```
