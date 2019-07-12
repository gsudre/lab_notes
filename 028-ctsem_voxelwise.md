# 2019-07-08 16:08:42

We first need to replace the variables in the file Philip sent by the voxels. In
other words, all we need to do is replace Y* by voxel data. The data needs to be
residualized after motion and Z scored.

I couldn't find the original .nii.gz diffeos that MArine used, so I'll just use
the dataframes she was using for voxelwise data. I'll later need to find the
templates used to extract the voxels, so I can put everything back as nifti.

```R
long<-read.csv('/Volumes/Labs/cross_lag/ctsem/for_gustavo/LONG_file_for_gustavo_3obs.csv', T)
load('/Volumes/Labs/marine/dti_crosslag/databases/raw/dti_fa_voxelwise_17JUNE2019.RData.gz')
nvox = 10814
header = sapply(1:nvox, function(x) {sprintf('Y%d', x)})
colnames(data) = c('maskid', header)
d2 = merge(long, data, by='maskid', all.x=F, all.y=F)
# renaming some variables with incorrect name
colnames(d2)[colnames(d2)=='sex'] = 'TI1'
colnames(d2)[colnames(d2)=='psychostim_group_3'] = 'TI2'
colnames(d2)[colnames(d2)=='scan_age'] = 'time'
colnames(d2)[colnames(d2)=='mrn'] = 'id'
# residualizing and z-scaling
for (v in header) {
    fm_str = sprintf('%s ~ motion', v)
    d2[, v] = scale(residuals(lm(fm_str, data=d2)))
}
wide_resid_motion<-ctLongToWide(datalong=d2, id="id", time="time",
                                manifestNames=c("sx_inatt", "sx_hi", header),
                                TIpredNames=c("TI1", "TI2"))
wider_resids_motion<-ctIntervalise(datawide=wide_resid_motion, Tpoints=3,
                                   n.manifest=(length(header)+2), n.TIpred=2,
                                   manifestNames=c("sx_inatt", "sx_hi", header),
                                   TIpredNames=c("TI1", "TI2"),
                                   individualRelativeTime=FALSE)
save(wider_resids_motion,
     file='~/data/ctsem_voxelwise/FA_wider_3obs_developmental_time.RData.gz',
     compress=T)
```

Now we need to make a few changes to the ctsem_loop script to run many voxels,
ideally more efficiently.

With a future goal to eventually run permutations, I set it up so that each R
script can run a list of voxels (to avoid overloading the filesystem with
multiple R calls). Still, each R call only uses one thread, so we can use
gnu-parallel to run multiple R calls (multiple lists of voxels).

If we have 10814 voxels, then we can start with 10 voxels per list, which is 320
per node if we request 32-core nodes. Also note that a single run of a single
phenotype takes about 14 min, and it should be safe to just multiple that by the
length of the phenotype list, as after the first one no package or data needs to
be loaded again. Also, each thread is only taking .5Gb.

So, in the quick partition I can run a job for 4h. Say, 12 voxels per thread to
be in the safe side, or increments of 384 per job.

```bash
cd ~/data/ctsem_voxelwise;
nvox=10814;
jname=ad_hi;
fname=swarm.ctsem;
bundle=384;
cur_vox=1;
while [ $cur_vox -lt $nvox ]; do
    let last_vox=${cur_vox}+${bundle}-1;
    # gets the min
    last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
    echo "bash ~/research_code/run_ctsem_voxel_parallel.sh ~/data/ctsem_voxelwise/AD_wider_3obs_developmental_time.RData.gz sx_hi ${cur_vox} ${last_vox}" >> $fname;
    let cur_vox=${last_vox}+1;
done;
swarm --gres=lscratch:1 -f ${fname} --module R -g 20 -t 32 \
    --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
    --partition quick;
```

This current setup needs 29 jobs, which can be all done in less than 4h if
executed in parallel.

Note that I've been able to run 171 of such jobs in parallel, which is quite
nice... it will make things better for future permutations, or it gives me some
room to play with the length of each list. Less than 12 will run faster, but
also execute more R instances.

Because of that, if we want to go nuts, we can do:

```bash
cd ~/data/ctsem_voxelwise;
nvox=10814;
jname=all_ctsem;
fname=swarm.ctsem;
bundle=384;
rm $fname;
for sx in inatt hi; do
    for p in FA AD RD; do
        cur_vox=1;
        while [ $cur_vox -lt $nvox ]; do
            let last_vox=${cur_vox}+${bundle}-1;
            # gets the min
            last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
            echo "bash ~/research_code/run_ctsem_voxel_parallel.sh ~/data/ctsem_voxelwise/${p}_wider_3obs_developmental_time.RData.gz sx_${sx} ${cur_vox} ${last_vox}" >> $fname;
            let cur_vox=${last_vox}+1;
        done;
    done;
done
swarm --gres=lscratch:1 -f ${fname} --module R -g 30 -t 32 \
    --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
    --partition quick;
```

## compiling the results

When it's time to compile the results, I think it makes sense to create one .nii
file per drift, and then the bricks in the nifti file will be estimate, p-value,
and everything else we use to calculate it: std, ub, lb. I should also mark the
voxels that didn't converge as NA, and interpolate using nearest neighbors.

The output can be split with:

grep sx_inatt_sx_inatt ~/tmp/oi.csv
grep _sx_inatt_Y ~/tmp/oi.csv
grep "Y[0-9]\+_sx_inatt" ~/tmp/oi.csv
grep "Y[0-9]\+_Y[0-9]\+" ~/tmp/oi.csv
grep AIC ~/tmp/oi.csv
grep BIC ~/tmp/oi.csv
grep msg ~/tmp/oi.csv

# 2019-07-11 11:00:27

Apparently many of the jobs died because of no time... should I switch it to
norm?

I tested 2 voxels per core, and the whole thing took 1.5h, so it's more like
45min per voxel. Let's do 1h to be safe. That means we can do 4 voxels per core,
128 per 32-core node.

```bash
cd ~/data/ctsem_voxelwise;
nvox=10814;
jname=all_ctsem;
fname=swarm.ctsem;
bundle=128;
rm $fname;
for sx in inatt hi; do
    for p in FA AD RD; do
        cur_vox=1;
        while [ $cur_vox -lt $nvox ]; do
            let last_vox=${cur_vox}+${bundle}-1;
            # gets the min
            last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
            echo "bash ~/research_code/run_ctsem_voxel_parallel.sh ~/data/ctsem_voxelwise/${p}_wider_3obs_developmental_time.RData.gz sx_${sx} ${cur_vox} ${last_vox} ~/data/ctsem_voxelwise/TI1" >> $fname;
            let cur_vox=${last_vox}+1;
        done;
    done;
done
swarm --gres=lscratch:1 -f ${fname} --module R -g 30 -t 32 \
    --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
    --partition quick;
```

Then, we compile everything. First, generate the ijk:

```bash
3dmaskdump -mask /Volumes/Labs/marine/dti_crosslag/data/mean_FA_skeleton_mask.nii.gz \
    /Volumes/Labs/marine/dti_crosslag/data/mean_FA_skeleton_mask.nii.gz > ijk.txt;
cut -d " " -f 1,2,3 ijk.txt > ctsem_ijk.txt;
```

```bash
module load afni
rm -rf csv;
mkdir csv;
sx=hi
for f in `ls ~/data/ctsem_voxelwise/TI1/FA_*_${sx}*tgz`; do tar -zxf $f -C csv/; done
# don't output file name
grep -h sx_${sx}_sx_${sx} csv/*.csv > sx_sx.csv;
grep -h _sx_${sx}_Y csv/*.csv > sx_voxels.csv;
grep -h "Y[0-9]\+_sx_${sx}" csv/*.csv > voxels_sx.csv;
grep -h "Y[0-9]\+_Y[0-9]\+" csv/*.csv > voxel_voxel.csv;
grep -h AIC csv/*.csv > aic.csv;
grep -h BIC csv/*.csv > bic.csv;
grep -h msg csv/*.csv > msgs.csv;
python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
    sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
    ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

# and if it completes without reruns
for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
done
```

# 2019-07-12 09:44:39

When there are holes, we can re-run the list of voxels like this:

```bash
split -da 2 -l $((`wc -l < vlist_rerun.sx_sx` /$SLURM_CPUS_PER_TASK)) vlist_rerun.sx_sx vlist --additional-suffix=".txt";
ls -1 vlist*txt > file_list.txt;
cp ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R ./;
cat file_list.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ctsem_voxel_developmental_time_3_timepoints.R ~/data/ctsem_voxelwise/FA_wider_3obs_developmental_time.RData.gz sx_inatt {} output_{}.csv
tar -czf FA_wider_3obs_developmental_time_SX_inatt_redo.tgz output_*.csv;
cp *.tgz ~/data/ctsem_voxelwise/TI1/
```

# TODO
 * make sure we're setting a seed in the scripts
