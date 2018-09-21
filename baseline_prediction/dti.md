# 2018-09-21 13:07:33

Unfortunately the 263 scans I was using before showed some correlation with
norm.rot and norm.trans, when using the correct movement estimates. So, in
baseline_dti_corrs_redux.Rmd I came up with two new groups: 223 and 272. The 223
set uses the same strategy as before, but it's a bit more stringent and now
shows no correlation between the mean variables and movement. The 272 only gets
rid of scans based 2STD of the DTI distributions and transformed movement
variables, but then the idea is to run them as
is and just check the important variables to make sure they're not correlated
with movement.

In any case, neither set shows correlation between the target variables we're
running and movement.

So, let's create the DTI data to re-run the ML experiments for 223 and 272:

```bash
n=223;
rm ~/data/baseline_prediction/subjs_diffeo_${n}.txt;
# ignore header
for m in `tail -n +2 ~/data/baseline_prediction/dti_gf_09212018_${n}timeDiff12mo.csv | cut -d"," -f 2 -`; do
     maskid=`printf %04d $m`;
     echo ${maskid}_tensor_diffeo.nii.gz >> ~/data/baseline_prediction/subjs_diffeo_${n}.txt;
done
cd /Volumes/Shaw/dti_robust_tsa/analysis_may2017/
/Applications/dtitk-2.3.3-Darwin-x86_64/bin/TVMean -in ~/data/baseline_prediction/subjs_diffeo_${n}.txt -out ~/data/baseline_prediction/mean_subjects_${n}.nii.gz
cd ~/data/baseline_prediction/
/Applications/dtitk-2.3.3-Darwin-x86_64/bin/TVtool -in mean_subjects_${n}.nii.gz -fa
tbss_skeleton -i mean_subjects_${n}_fa -o mean_${n}_fa_skeleton
3dcalc -a mean_${n}_fa_skeleton.nii.gz -prefix mean_${n}_fa_skeleton_mask.nii.gz -expr 'step(a-.2)'
dti_dir=/Volumes/Shaw/dti_robust_tsa/analysis_may2017/
for m in `tail -n +2 ~/data/baseline_prediction/dti_gf_09212018_${n}timeDiff12mo.csv | cut -d"," -f 2 -`; do
     maskid=`printf %04d $m`;
     3dmaskdump -mask mean_${n}_fa_skeleton_mask.nii.gz -o ${maskid}_fa.txt ${dti_dir}/${maskid}_tensor_diffeo_fa.nii.gz;
     3dmaskdump -mask mean_${n}_fa_skeleton_mask.nii.gz -o ${maskid}_ad.txt ${dti_dir}/${maskid}_tensor_diffeo_ad.nii.gz;
     3dmaskdump -mask mean_${n}_fa_skeleton_mask.nii.gz -o ${maskid}_rd.txt ${dti_dir}/${maskid}_tensor_diffeo_rd.nii.gz;
done
mkdir dti_voxels_${n}
mv *fa.txt *ad.txt *rd.txt dti_voxels_${n}/ 
```

```r
library(gdata)
clin = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'dti')
clin = clin[, c(1,4)]
colnames(clin) = c('MRN', 'mask.id')
for (n in c(272, 223)) {
    fname = sprintf('~/data/baseline_prediction/dti_gf_09212018_%dtimeDiff12mo.csv', n)
    pheno = read.csv(fname)
    junk = read.table(sprintf('~/data/baseline_prediction/dti_voxels_%d/0282_fa.txt', n))
    nvox = nrow(junk)
    for (m in c('fa', 'ad', 'rd')) {
        print(m)
        dti_data = matrix(nrow=nrow(pheno), ncol=nvox)
        for (s in 1:nrow(dti_data)) {
            a = read.table(sprintf('~/data/baseline_prediction/dti_voxels_%d/%04d_%s.txt',
                                   n, pheno[s,]$mask.id, m))
            dti_data[s, ] = a[, 4]
        }
        data = cbind(pheno$mask.id, dti_data)
        cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%05d', d)))
        colnames(data) = cnames
        data = merge(clin, data, by='mask.id')
        save(data, file=sprintf('~/data/baseline_prediction/dti_%s_voxelwise_n%d_09212018.RData.gz', m, n),
             compress=T)
    }
}
```

Just to check for higher order interactions between the DTI modalities, I also created an uber file:

```r
for (n in c(272, 223)) {
    load(sprintf('~/data/baseline_prediction/dti_fa_voxelwise_n%d_09212018.RData.gz', n))
    all_data = data
    load(sprintf('~/data/baseline_prediction/dti_ad_voxelwise_n%d_09212018.RData.gz', n))
    all_data = merge(all_data, data, by='MRN', suffixes=c('.fa', '.ad'))
    load(sprintf('~/data/baseline_prediction/dti_rd_voxelwise_n%d_09212018.RData.gz', n))
    data = merge(all_data, data, by='MRN', suffixes=c('', '.rd'))
    save(data, file=sprintf('~/data/baseline_prediction/dti_ALL_voxelwise_n%d_09212018.RData.gz', n),
         compress=T)
}
```

And we need to generate the tracts as well, even though they are quite
useless...

```bash
cd /Volumes/Shaw/dti_robust_tsa/analysis_may2017/
export DTITK_ROOT=/Applications/dtitk-2.3.3-Darwin-x86_64/
/Applications/dtitk-2.3.3-Darwin-x86_64/scripts/tsa_sampling ~/data/baseline_prediction/subjs_diffeo_223.txt ../ixi_aging_template_v3.0/tsa/ mean
python ~/research_code/lab_mgmt/convert_dti_sampling.py
Rscript ~/research_code/dti/compile_tract_table.R
```

```r
library(gdata)
clin = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'dti')
clin = clin[, c(1,4)]
colnames(clin) = c('MRN', 'mask.id')
for (n in c(272, 223)) {
    fname = sprintf('~/data/baseline_prediction/dti_mean_phenotype_%d.csv', n)
    data = read.csv(fname)
    data = merge(clin, data, by.x='mask.id', by.y='file')
    save(data, file=sprintf('~/data/baseline_prediction/dti_tracts_n%d_09212018.RData.gz', n),
             compress=T)
}
```

Now it's time to re-run it.

```bash
rm swarm.automl
for var in raw pca uni pca_uni uni_pca uni01; do
    for n in 272 223; do
        for m in fa ad rd ALL tracts; do
            for sx in inatt HI total; do
                for sl in OLS random; do
                    echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/dti_${m}_voxelwise_n${n}_09212018.RData.gz ~/data/baseline_prediction/long_clin_0918.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
                done;
            done;
            for sx in INATT HI total; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/dti_${m}_voxelwise_n${n}_09212018.RData.gz ~/data/baseline_prediction/long_clin_0918.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
            done;
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/dti_${m}_voxelwise_n${n}_09212018.RData.gz ~/data/baseline_prediction/long_clin_0918.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
        done;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl;
cp swarm.automl swarm.automl_dti;
swarm -f swarm.automl_dti -g 40 -t 32 --time 4:00:00 --partition quick --logdir trash_bin --job-name dti_redo -m R --gres=lscratch:10
```

And we also re-run the subgroups:

```bash
rm swarm.automl_subgroupDTI
for var in raw pca uni pca_uni uni_pca uni01; do
    for f in `/bin/ls -1 ~/data/baseline_prediction/dti*0921*gz`;do
        for nn in nonew_ ''; do
            for g in ADHDonly_ ''; do
                for t in diag_group2 OLS_inatt_slope OLS_HI_slope OLS_total_slope; do
                    echo "Rscript --vanilla ~/research_code/automl/${var}.R $f ~/data/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} ~/tmp/${nn}${g}${t}" >> swarm.automl_subgroupDTI;
                done; 
            done;
            for g in ADHDNOS_ ADHDNOS_group; do
                for t in OLS_inatt_slope OLS_HI_slope OLS_total_slope; do
                    echo "Rscript --vanilla ~/research_code/automl/${var}.R $f ~/data/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} ~/tmp/${nn}${g}${t}" >> swarm.automl_subgroupDTI;
                done; 
            done;
        done;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl_subgroupDTI;
split -l 900 swarm.automl_subgroupDTI;
swarm -f xaa -g 40 -t 32 --time 4:00:00 --partition quick --logdir trash_bin --job-name subgroupDTI -m R --gres=lscratch:10;
swarm -f xab -g 40 -t 32 --time 4:00:00 --partition quick --logdir trash_bin --job-name subgroupDTI -m R --gres=lscratch:10;
```

RUNNING NEW DTI ON FULL DATA AND SUBGROUPS!

Results too good to be true? getting very close to 100% accuracy for subgroups, but that's
true even when using leaderboard scoring, and forced 5-fold CV... look at
uni01.R for tests... Also, nice that loading saved models bring up the
cross-validation scores. Maybe plot it in the brain to see if it makes sense,
implement some dummy scoring, and calculate univariate filters only in training
set? 