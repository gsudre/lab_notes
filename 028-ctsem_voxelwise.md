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

So, since I'll be leaving this to run over the weekend...

```bash
module load afni
module load R
sx=hi
TI=TI1_TI2
cd /lscratch/$SLURM_JOBID
for m in FA AD RD; do
    rm -rf *;
    mkdir csv;
    for f in `ls ~/data/ctsem_voxelwise/${TI}/${m}_*_${sx}*tgz`; do tar -zxf $f -C csv/; done
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

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        split -da 2 -l $((`wc -l < vlist_rerun.sx_sx` /$SLURM_CPUS_PER_TASK)) vlist_rerun.sx_sx vlist --additional-suffix=".txt";
        ls -1 vlist*txt > file_list.txt;
        # made sure that the script is setup for TI1 now!
        cat file_list.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
            Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R \
                ~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz \
                sx_${sx} {} output_{}.csv;
        tar -czf ${m}_wider_3obs_developmental_time_SX_${sx}_redo.tgz output_*.csv;
        cp *.tgz ~/data/ctsem_voxelwise/${TI}/;
    fi;

    # and assuming no other reruns are needed
    rm csv/*;
    for f in `ls ~/data/ctsem_voxelwise/${TI}/${m}_*_${sx}*tgz`; do tar -zxf $f -C csv/; done
    grep -h sx_${sx}_sx_${sx} csv/*.csv > sx_sx.csv;
    grep -h _sx_${sx}_Y csv/*.csv > sx_voxels.csv;
    grep -h "Y[0-9]\+_sx_${sx}" csv/*.csv > voxels_sx.csv;
    grep -h "Y[0-9]\+_Y[0-9]\+" csv/*.csv > voxel_voxel.csv;
    grep -h AIC csv/*.csv > aic.csv;
    grep -h BIC csv/*.csv > bic.csv;
    grep -h msg csv/*.csv > msgs.csv;
    for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
        python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
            ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
            ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
    done
    mkdir ~/data/ctsem_voxelwise/${TI}/${m}_${sx}/;
    cp *nii.gz ~/data/ctsem_voxelwise/${TI}/${m}_${sx}/;
done
```

# 2019-07-16 09:37:58

We finally have voxelwise results for most of the combinations. So, let's see
what we get:

```bash
3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
```

For inatt, the cluster FA and RD were bigger for sx_voxels than voxels_sx (twice
the size). But for AD, it was bigger for voxels_sx.

```
HG-02113362-DM4:FA_inatt sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file sx_voxels.nii.gz
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities:
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     58  -106.7  -94.2   87.2  -120.0  -98.0  -102.0  -84.0   84.0   90.0   0.1207   0.0025    -0.17  -102.0  -92.0   88.0
     36  -123.7  -107.5   60.7  -128.0  -120.0  -114.0  -102.0   54.0   68.0    0.128   0.0031  -0.1681  -124.0  -106.0   60.0
     25  -131.0  -78.8   94.7  -132.0  -130.0  -88.0  -72.0   92.0   98.0    0.091   0.0055  -0.1356  -130.0  -80.0   94.0
     22  -134.1  -132.9   81.4  -136.0  -134.0  -138.0  -126.0   76.0   86.0   0.1296   0.0062  -0.2003  -134.0  -136.0   78.0
     21  -135.7  -127.6   61.8  -138.0  -134.0  -134.0  -122.0   60.0   66.0   0.1012   0.0077  -0.1882  -134.0  -124.0   60.0
     21  -85.4  -86.6   77.8  -86.0  -82.0  -92.0  -82.0   72.0   84.0   0.1134   0.0057  -0.2013  -86.0  -90.0   76.0
     16  -101.8  -107.4   57.7  -104.0  -100.0  -112.0  -104.0   54.0   64.0   0.1259    0.004  -0.1681  -104.0  -106.0   54.0
HG-02113362-DM4:FA_inatt sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 voxels_sx.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file voxels_sx.nii.gz
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities:
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     25  -135.8  -131.4   68.2  -138.0  -134.0  -136.0  -126.0   64.0   76.0   0.0077  4.1e-04   0.0114  -134.0  -136.0   70.0
     21  -128.0  -137.8   97.0  -128.0  -128.0  -144.0  -132.0   94.0  102.0   0.0053  2.3e-04   0.0074  -128.0  -138.0   96.0
     19  -142.0  -110.7   71.8  -142.0  -142.0  -116.0  -104.0   68.0   76.0   0.0111  5.7e-04   0.0164  -142.0  -114.0   68.0
     18  -134.0  -134.3   80.2  -134.0  -134.0  -138.0  -130.0   74.0   86.0     0.01  3.9e-04   0.0137  -134.0  -136.0   80.0
     16  -126.0  -141.8  101.1  -126.0  -126.0  -150.0  -132.0   98.0  104.0   0.0062  2.5e-04   0.0076  -126.0  -144.0  102.0
     16  -128.0  -113.1  107.1  -128.0  -128.0  -124.0  -104.0  106.0  108.0   0.0045  3.1e-04   0.0069  -128.0  -114.0  108.0
     16  -95.8  -89.0  120.6  -100.0  -90.0  -90.0  -88.0  118.0  124.0   0.0062  3.2e-04   0.0082  -98.0  -88.0  122.0
     15  -123.7  -85.4  117.8  -126.0  -122.0  -88.0  -82.0  114.0  122.0   0.0069  4.0e-04   0.0109  -122.0  -84.0  122.0

HG-02113362-DM4:AD_inatt sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file sx_voxels.nii.gz
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities:
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     24  -144.9  -88.0   65.5  -148.0  -142.0  -94.0  -80.0   62.0   70.0   0.1162   0.0042  -0.1678  -144.0  -92.0   64.0
     19  -112.1  -93.0   54.4  -118.0  -106.0  -94.0  -90.0   52.0   58.0   0.1108   0.0045  -0.1472  -112.0  -92.0   54.0
     15  -87.6  -95.2   70.0  -90.0  -82.0  -98.0  -92.0   66.0   74.0   0.1066   0.0071  -0.1616  -88.0  -94.0   72.0
HG-02113362-DM4:AD_inatt sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 voxels_sx.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file voxels_sx.nii.gz
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities:
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     44  -90.0  -93.7   96.6  -92.0  -88.0  -102.0  -86.0   90.0  106.0   0.0083  3.0e-04   0.0131  -88.0  -92.0   96.0
     27  -87.7  -93.0   79.0  -88.0  -86.0  -98.0  -88.0   72.0   88.0   0.0089  2.8e-04    0.011  -88.0  -92.0   76.0
     27  -84.7  -81.6   81.1  -90.0  -82.0  -84.0  -76.0   74.0   88.0   0.0084  3.0e-04   0.0123  -84.0  -82.0   80.0
     19  -133.6  -94.4   95.7  -136.0  -130.0  -98.0  -88.0   94.0   98.0   0.0084  2.8e-04   0.0104  -134.0  -98.0   98.0
     16  -132.0  -89.9   99.9  -132.0  -132.0  -94.0  -84.0   96.0  104.0   0.0078  4.4e-04   0.0112  -132.0  -90.0   96.0
HG-02113362-DM4:AD_inatt sudregp$ cd ../RD_inatt
HG-02113362-DM4:RD_inatt sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file sx_voxels.nii.gz
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities:
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     50  -105.1  -95.2   87.8  -114.0  -98.0  -106.0  -84.0   84.0   92.0   0.1182   0.0026   0.1576  -102.0  -88.0   86.0
     18  -130.9  -79.7   94.2  -132.0  -130.0  -84.0  -74.0   92.0   96.0   0.1167   0.0041   0.1515  -130.0  -78.0   94.0
     15  -102.0  -108.3   58.8  -104.0  -100.0  -112.0  -104.0   54.0   66.0   0.1253   0.0039   0.1592  -102.0  -108.0   60.0
     15  -134.0  -132.8   81.1  -134.0  -134.0  -136.0  -128.0   76.0   86.0   0.1297   0.0041    0.157  -134.0  -132.0   80.0
     15  -88.6  -125.5   92.3  -90.0  -88.0  -130.0  -120.0   88.0   96.0   0.1026   0.0024   0.1201  -90.0  -122.0   96.0
HG-02113362-DM4:RD_inatt sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 voxels_sx.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file voxels_sx.nii.gz
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities:
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     27  -141.9  -61.8   30.0  -150.0  -134.0  -72.0  -56.0   30.0   30.0   0.0125  5.1e-04   0.0183  -140.0  -56.0   30.0
     25  -155.2  -76.7   76.6  -162.0  -152.0  -80.0  -74.0   74.0   80.0   0.0072  4.3e-04   0.0123  -154.0  -76.0   78.0
     19  -69.9  -78.8   90.2  -72.0  -66.0  -82.0  -78.0   86.0   96.0   0.0066  3.1e-04   0.0095  -68.0  -78.0   88.0
     19  -151.3  -88.3   92.0  -158.0  -146.0  -92.0  -84.0   92.0   92.0   0.0096  5.5e-04   0.0133  -148.0  -88.0   92.0
     18  -145.9  -88.4   64.6  -148.0  -144.0  -96.0  -78.0   62.0   68.0   0.0082  4.1e-04    0.012  -146.0  -94.0   62.0
```

In HI, it looked like cluster for voxels_sx was bigger in FA and AD, but
sx_voxels was bigger in RD. 

```
HG-02113362-DM4:RD_inatt sudregp$ cd ../FA_hi
HG-02113362-DM4:FA_hi sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file sx_voxels.nii.gz 
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities: 
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     33  -75.6  -91.3   59.3  -78.0  -72.0  -102.0  -80.0   54.0   66.0   0.1021   0.0028  -0.1447  -74.0  -96.0   56.0 
     22  -127.8  -149.0   77.6  -128.0  -126.0  -152.0  -146.0   70.0   86.0   0.0942   0.0028  -0.1231  -126.0  -146.0   72.0 
     20  -88.0  -89.6   87.3  -88.0  -88.0  -98.0  -84.0   82.0   90.0   0.0947   0.0025  -0.1194  -88.0  -90.0   88.0 
     18  -132.7  -51.7   84.3  -134.0  -132.0  -60.0  -44.0   82.0   86.0   0.0886   0.0018     -0.1  -132.0  -46.0   86.0 
     17  -81.5  -70.7   77.0  -82.0  -80.0  -78.0  -66.0   72.0   80.0   0.1006   0.0039  -0.1431  -80.0  -68.0   74.0 
     15  -97.9  -101.4   56.3  -100.0  -94.0  -104.0  -98.0   52.0   60.0    0.106   0.0045  -0.1375  -98.0  -102.0   58.0 
HG-02113362-DM4:FA_hi sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 voxels_sx.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file voxels_sx.nii.gz 
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities: 
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     50  -81.7  -98.9  107.8  -88.0  -74.0  -104.0  -94.0  102.0  116.0   0.0105  3.2e-04   0.0162  -84.0  -100.0  104.0 
     44  -125.8  -86.7  116.5  -132.0  -120.0  -92.0  -82.0  108.0  124.0   0.0094  3.2e-04   0.0147  -122.0  -84.0  122.0 
     41  -147.1  -85.3  108.2  -152.0  -140.0  -90.0  -82.0  102.0  114.0   0.0079  2.0e-04   0.0112  -144.0  -84.0  106.0 
     40  -133.9  -138.5   78.3  -136.0  -132.0  -148.0  -128.0   70.0   86.0   0.0131  4.1e-04    0.018  -134.0  -136.0   80.0 
     35  -147.6  -100.5   57.6  -152.0  -144.0  -106.0  -94.0   52.0   62.0   0.0133  6.0e-04  -0.0204  -148.0  -102.0   56.0 
     28  -87.2  -131.2   67.8  -92.0  -84.0  -136.0  -126.0   62.0   74.0   0.0096  3.4e-04   0.0126  -88.0  -132.0   72.0 
     28  -88.6  -94.9   90.3  -90.0  -88.0  -98.0  -88.0   82.0   98.0    0.011  3.6e-04   0.0143  -88.0  -96.0   88.0 
     25  -131.2  -77.6   94.6  -132.0  -130.0  -84.0  -72.0   90.0   98.0   0.0075  2.9e-04   0.0106  -132.0  -74.0   96.0 
HG-02113362-DM4:FA_hi sudregp$ cd ../AD_hi
HG-02113362-DM4:AD_hi sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file sx_voxels.nii.gz 
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities: 
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     15  -94.0  -153.7   77.8  -94.0  -94.0  -156.0  -150.0   74.0   82.0    0.122   0.0044   0.1442  -94.0  -152.0   78.0 
HG-02113362-DM4:AD_hi sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 voxels_sx.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file voxels_sx.nii.gz 
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities: 
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
    111  -133.0  -97.9   95.6  -136.0  -128.0  -116.0  -82.0   80.0  108.0   0.0113  2.4e-04   0.0198  -134.0  -112.0   96.0 
     72  -90.3  -92.2   95.6  -94.0  -88.0  -98.0  -84.0   80.0  108.0   0.0106  2.9e-04   0.0148  -92.0  -90.0  104.0 
     40  -126.8  -149.3   69.6  -128.0  -124.0  -158.0  -144.0   62.0   78.0   0.0151  4.1e-04   0.0206  -126.0  -148.0   68.0 
     30  -79.9  -91.3   66.8  -84.0  -76.0  -96.0  -86.0   60.0   74.0   0.0123  6.7e-04   0.0201  -80.0  -92.0   68.0 
     28  -86.7  -128.9   63.9  -92.0  -84.0  -136.0  -122.0   60.0   70.0   0.0118  5.2e-04   0.0174  -92.0  -134.0   64.0 
     26  -142.5  -89.4   67.0  -148.0  -138.0  -94.0  -82.0   64.0   72.0   0.0126  5.7e-04   0.0186  -140.0  -90.0   68.0 
HG-02113362-DM4:AD_hi sudregp$ cd ../RD_hi
HG-02113362-DM4:RD_hi sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 sx_voxels.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file sx_voxels.nii.gz 
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities: 
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     89  -88.2  -96.3   90.7  -92.0  -84.0  -114.0  -78.0   78.0   98.0   0.1254   0.0022   0.1713  -88.0  -100.0   92.0 
     51  -127.3  -149.2   73.6  -128.0  -126.0  -154.0  -146.0   62.0   88.0   0.1183   0.0022   0.1465  -126.0  -148.0   62.0 
     48  -131.3  -87.2  103.5  -136.0  -128.0  -98.0  -82.0   94.0  110.0   0.1032   0.0019   0.1393  -132.0  -86.0  104.0 
     34  -134.2  -106.9   95.5  -136.0  -132.0  -116.0  -98.0   88.0   98.0   0.1349   0.0031   0.1671  -134.0  -106.0   94.0 
     31  -135.4  -90.3   89.8  -136.0  -130.0  -100.0  -84.0   86.0   92.0   0.1207   0.0038   0.1641  -136.0  -92.0   90.0 
     30  -75.0  -93.0   58.0  -78.0  -72.0  -102.0  -84.0   54.0   62.0   0.1143   0.0037   0.1595  -74.0  -94.0   56.0 
     28  -86.5  -130.6   67.1  -88.0  -84.0  -134.0  -126.0   62.0   72.0   0.0929   0.0026   0.1277  -86.0  -130.0   66.0 
     28  -93.0  -88.6  110.6  -94.0  -92.0  -92.0  -84.0  102.0  122.0    0.111   0.0036    0.153  -92.0  -84.0  106.0 
     26  -141.1  -148.6   76.4  -148.0  -136.0  -150.0  -146.0   70.0   82.0   0.1213   0.0047   0.1767  -138.0  -148.0   78.0 
     23  -148.3  -96.1   92.1  -154.0  -146.0  -104.0  -90.0   92.0   94.0   0.1209   0.0033   0.1535  -148.0  -96.0   92.0 
     22  -142.0  -108.6   71.3  -142.0  -142.0  -114.0  -102.0   68.0   76.0   0.1573   0.0064   0.2118  -142.0  -112.0   70.0 
     16  -84.1  -97.6  107.3  -86.0  -82.0  -100.0  -96.0  104.0  112.0   0.1077   0.0035   0.1349  -86.0  -98.0  104.0 
     15  -131.3  -50.3   85.8  -132.0  -130.0  -56.0  -46.0   84.0   88.0   0.0992   0.0032   0.1224  -132.0  -48.0   86.0 
HG-02113362-DM4:RD_hi sudregp$ 3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 voxels_sx.nii.gz
++ 3dclust: AFNI version=AFNI_19.1.09 (May 17 2019) [64-bit]
++ Authored by: RW Cox et alii
#
#Cluster report for file voxels_sx.nii.gz 
#[Connectivity radius = 1.11 mm  Volume threshold = 120.00 ]
#[Single voxel volume = 8.0 (microliters) ]
#[Voxel datum type    = float ]
#[Voxel dimensions    = 2.000 mm X 2.000 mm X 2.000 mm ]
#[Coordinates Order   = RAI ]
#[Fake voxel dimen    = 1.000 mm X 1.000 mm X 1.000 mm ]
#Mean and SEM based on Absolute Value of voxel intensities: 
#
#Volume  CM RL  CM AP  CM IS  minRL  maxRL  minAP  maxAP  minIS  maxIS    Mean     SEM    Max Int  MI RL  MI AP  MI IS
#------  -----  -----  -----  -----  -----  -----  -----  -----  -----  -------  -------  -------  -----  -----  -----
     44  -80.4  -99.8  107.6  -86.0  -72.0  -104.0  -94.0  100.0  116.0   0.0114  3.4e-04  -0.0171  -76.0  -104.0  104.0 
     44  -147.1  -85.1  108.7  -152.0  -138.0  -90.0  -82.0  102.0  114.0   0.0106  3.2e-04  -0.0152  -146.0  -84.0  106.0 
     44  -75.5  -86.7  110.8  -84.0  -70.0  -94.0  -82.0  104.0  116.0   0.0106  3.7e-04   -0.017  -80.0  -84.0  112.0 
     40  -146.8  -98.9   59.1  -152.0  -142.0  -106.0  -90.0   54.0   64.0   0.0162  6.1e-04   0.0237  -146.0  -98.0   60.0 
     27  -123.3  -84.4  118.6  -128.0  -120.0  -88.0  -78.0  114.0  124.0   0.0112  3.8e-04  -0.0161  -120.0  -82.0  122.0 

```

Those were all results using TI1, and I'm still waiting on TI1_TI2. Also, note
that it's all TORTOISE preprocessing, and it's not doing anything about the
masks. I need to check whether masked results have the same pattern, and also
create the script for the imputed version. But let's see where those bigger
clusters are.

Of course, the size of the cluster just impacts its significance, so it could be
that all those results are good.

To visualize where they are, I'll need to put them in MNI space.

```bash
flirt -in /Volumes/Labs/marine/dti_crosslag/data/mean_FA_skeleton.nii.gz \
    -ref /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii \
    -out ./mean_FA_skeleton_IN_HCP1065.nii.gz \
    -omat group_skeleton_to_HCP1065.mat -bins 256 -cost corratio \
    -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear

3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 \
    -overwrite -savemask mymask.nii sx_voxels.nii.gz

flirt -in mymask.nii \
    -ref /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii \
    -out mymask_inHCP1065.nii.gz -applyxfm \
    -init ../../group_skeleton_to_HCP1065.mat -interp nearestneighbour

# just to get the COM coordinates in the new space
3dclust -1Dformat -nosum -orient LPI -NN1 5 mymask_inHCP1065.nii.gz
```

OK, so this is working. Before we spend much more time here, let's work on the
permutation scripts.

There was always a series of voxels for which the whole thing didn't complete
before. So, let's make sure we give it enough time now:

```bash
cd ~/data/ctsem_voxelwise;
nvox=10814;
bundle=128;
p=FA;
sx=inatt;
for i in `seq 1 10`; do
    jname=${p}_${sx}_p${i};
    fname=swarm.${jname};
    rm $fname;
    seed=${RANDOM};
    cur_vox=1;
    while [ $cur_vox -lt $nvox ]; do
        let last_vox=${cur_vox}+${bundle}-1;
        # gets the min
        last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
        echo "bash ~/research_code/run_ctsem_voxel_parallel.sh" \
            "~/data/ctsem_voxelwise/${p}_wider_3obs_developmental_time.RData.gz" \
            "sx_${sx} ${cur_vox} ${last_vox}" \
            "~/data/ctsem_voxelwise/TI1_perms $seed" >> $fname;
        let cur_vox=${last_vox}+1;
    done;
    swarm --gres=lscratch:1 -f ${fname} --module R -g 30 -t 32 \
        --logdir=trash_${jname} --job-name ${jname} --time=8:00:00 --merge-output \
        --partition norm;
done
```

If this is taking too long to allocate, I could also try bundles of 96 (3 per
core), and run for 4h in the quick partition?

# 2019-07-17 11:42:11

For compiling permutation runs we do:

```bash
module load afni
module load R
module load python/3.7
sx=inatt
m=FA
TI=TI1_perms
cd /lscratch/$SLURM_JOBID
for p in `seq 3 20`; do
    echo "Working on perm $p...";
    rm -rf *;
    mkdir csv;
    pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.${m}_${sx}_p${p} | head -n 1`;
    ls ~/data/ctsem_voxelwise/${TI}/${m}_*_${sx}*_${pseed}.tgz > perm_files.txt;
    ls ~/data/ctsem_voxelwise/${TI}/*${pseed}*redo* >> perm_files.txt;
    for f in `cat perm_files.txt`; do
        tar -zxf $f -C csv/;
    done
    grep -h sx_${sx}_sx_${sx} csv/*.csv > sx_sx.csv;
    grep -h msg csv/*.csv > msgs.csv;
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        split -da 2 -l $((`wc -l < vlist_rerun.sx_sx` /$SLURM_CPUS_PER_TASK)) vlist_rerun.sx_sx vlist --additional-suffix=".txt";
        ls -1 vlist*txt > file_list.txt;
        # made sure that the script is setup for TI1 now!
        cat file_list.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
            Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R \
                ~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz \
                sx_${sx} {} output_{}.csv;
        tar -czf ${m}_wider_3obs_developmental_time_SX_${sx}_${pseed}_redo.tgz output_*.csv;
        cp *.tgz ~/data/ctsem_voxelwise/${TI}/;
    fi;

    # and assuming no other reruns are needed
    rm csv/*;
    rm -f sx_sx*.nii.gz;
    ls ~/data/ctsem_voxelwise/${TI}/${m}_*_${sx}*_${pseed}.tgz > perm_files.txt;
    ls ~/data/ctsem_voxelwise/${TI}/*${pseed}*redo* >> perm_files.txt;
    for f in `cat perm_files.txt`; do
        tar -zxf $f -C csv/;
    done
    grep -h sx_${sx}_sx_${sx} csv/*.csv > sx_sx.csv;
    grep -h _sx_${sx}_Y csv/*.csv > sx_voxels.csv;
    grep -h "Y[0-9]\+_sx_${sx}" csv/*.csv > voxels_sx.csv;
    grep -h "Y[0-9]\+_Y[0-9]\+" csv/*.csv > voxel_voxel.csv;
    grep -h AIC csv/*.csv > aic.csv;
    grep -h BIC csv/*.csv > bic.csv;
    grep -h msg csv/*.csv > msgs.csv;
    for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
        echo $f perm $p;
        python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
            ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
            ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
    done
    # if compiled successfully
    if [ -e sx_voxels.nii.gz ]; then
        echo "Cleaning up...";
        mkdir -p ~/data/ctsem_voxelwise/${TI}/${m}_${sx}_p${p}/;
        # reducing the number of files a bit
        tar -zcf ${m}_${sx}_p$p.tgz csv/*;
        cp -v *nii.gz ${m}_${sx}_p$p.tgz ~/data/ctsem_voxelwise/${TI}/${m}_${sx}_p${p}/;
        for f in `cat perm_files.txt`; do
            rm -f $f;
        done;
    fi;
done
```

# 2019-07-18 12:18:22

Ben brought up a good point in that we should have smaller bundles so that if
something fails, it won't take too long to re-run. But they can't be too small,
otherwise we're starting too many Rs. Let's try to find a good balance:

```bash
cd ~/data/ctsem_voxelwise;
nvox=10814;
bundle=24;
p=FA;
sx=inatt;
for i in `seq 60 70`; do
    jname=${p}_${sx}_p${i};
    fname=swarm.${jname};
    rm $fname;
    seed=${RANDOM};
    cur_vox=1;
    while [ $cur_vox -lt $nvox ]; do
        let last_vox=${cur_vox}+${bundle}-1;
        # gets the min
        last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
        echo "bash ~/research_code/run_ctsem_voxel_parallel.sh" \
            "~/data/ctsem_voxelwise/${p}_wider_3obs_developmental_time.RData.gz" \
            "sx_${sx} ${cur_vox} ${last_vox}" \
            "~/data/ctsem_voxelwise/TI1_perms $seed" >> $fname;
        let cur_vox=${last_vox}+1;
    done;
    echo "ERROR" > swarm_wait;
    while grep -q ERROR swarm_wait; do
        echo "Trying $jname"
        swarm --gres=lscratch:1 -f ${fname} --module R -g 10 -t 8 \
            --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
            --partition quick,norm 2> swarm_wait;
        if grep -q ERROR swarm_wait; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Just be careful because we start running into scheduling issues after 6 or 7 of
these... I even tried to establish dependencies, but it didn't work. So, maybe
the usual swarm wait style of submission might work better...

# 2019-07-19 09:39:45

Because we don't write to the network drive that often, we might be able to get
away writing on csv per voxel. Let's try that:

```bash
cd ~/data/ctsem_voxelwise;
nvox=10814;
bundle=24;#96;
p=AD;
sx=inatt;
for i in `seq 19 100`; do
    jname=${p}_${sx}_p${i};
    fname=swarm.${jname};
    rm -f $fname;
    seed=${RANDOM};
    cur_vox=1;
    while [ $cur_vox -lt $nvox ]; do
        let last_vox=${cur_vox}+${bundle}-1;
        # gets the min
        last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
        echo "bash ~/research_code/run_ctsem_voxel_parallel_net.sh" \
            "~/data/ctsem_voxelwise/${p}_wider_3obs_developmental_time.RData.gz" \
            "sx_${sx} ${cur_vox} ${last_vox}" \
            "~/data/ctsem_voxelwise/TI1_perms_${p}_${sx}_p${i} $seed" >> $fname;
        let cur_vox=${last_vox}+1;
    done;
    echo "ERROR" > swarm_wait;
    while grep -q ERROR swarm_wait; do
        echo "Trying $jname"
        # swarm --gres=lscratch:1 -f ${fname} --module R -g 30 -t 32 \
        #     --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
        #     --partition quick,norm 2> swarm_wait;
        swarm --gres=lscratch:1 -f ${fname} --module R -g 10 -t 8 \
            --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
            --partition quick,norm 2> swarm_wait;
        if grep -q ERROR swarm_wait; then
            echo -e "\tError, sleeping..."
            sleep 30m;
        fi;
    done;
done
```

# 2019-07-23 14:48:58

So, we have the entire result, masked impute, and masked zero. Ideally they
would all point in the same direction, but there are good chances that won't
happen. So, let's find the COM coordinates for all those cases, and leave the
figures for when I'm back in the office, with two monitors.

```bash
res_file=~/tmp/all_res.txt;
for m in FA AD RD; do
    for sx in inatt hi; do
        mydir=~/data/ctsem_voxelwise/results/TI1/${m}_${sx}/
        for f in voxels_sx sx_voxels; do
            for suf in '' '_maskImpute' '_maskZero'; do
                fname=${mydir}/${f}${suf}.nii.gz;
                echo "======= $m $sx $fname ======" >> $res_file;
                3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 \
                    -1thresh 0.95 -NN1 15 $fname >> $res_file;
            done;
        done;
    done;
done;
```

I copied the result file to ~/data/ctsem_voxelwise/results/TI1. The maskImpute
version is looking alright, and it's much more similar to the regular results
than the maskZero ones. In any case, it's not too hard to test them all. So far,
I only have compiled a few of the FA results, so let's start with them.

| mode | sx  | drift  |  mask | top cluster NN1 |
|---|---|---|---|---|
| FA  | inatt | voxels_sx  | none  | 25 | 
| FA  | inatt | voxels_sx  | impute  | 28 |
| FA  | inatt | voxels_sx  | zero | 18 | 
| FA  | inatt | sx_voxels  | none  | 58 | 
| FA  | inatt | sx_voxels  | impute  | 59 |
| FA  | inatt | sx_voxels  | zero | 20 | 
| FA  | hi | voxels_sx  | none  | 50 | 
| FA  | hi | voxels_sx  | impute  | 51 |
| FA  | hi | voxels_sx  | zero | 46 | 
| FA  | hi | sx_voxels  | none  | 33 | 
| FA  | hi | sx_voxels  | impute  | 33 |
| FA  | hi | sx_voxels  | zero | 32 | 
| AD  | inatt | voxels_sx  | none  | 44 | 
| AD  | inatt | voxels_sx  | impute  | 92 |
| AD  | inatt | voxels_sx  | zero | 15 | 
| AD  | inatt | sx_voxels  | none  | 24 | 
| AD  | inatt | sx_voxels  | impute  | 20 |
| AD  | inatt | sx_voxels  | zero | 0 | 
| AD  | hi | voxels_sx  | none  | 111 | 
| AD  | hi | voxels_sx  | impute  | 88 |
| AD  | hi | voxels_sx  | zero | 24 | 
| AD  | hi | sx_voxels  | none  | 15 | 
| AD  | hi | sx_voxels  | impute  | 26 |
| AD  | hi | sx_voxels  | zero | 0 | 
| RD  | inatt | voxels_sx  | none  | 27 | 
| RD  | inatt | voxels_sx  | impute  | 28 |
| RD  | inatt | voxels_sx  | zero | 0 | 
| RD  | inatt | sx_voxels  | none  | 50 | 
| RD  | inatt | sx_voxels  | impute  | 54 |
| RD  | inatt | sx_voxels  | zero | 25 | 
| RD  | hi | voxels_sx  | none  | 44 | 
| RD  | hi | voxels_sx  | impute  | 45 |
| RD  | hi | voxels_sx  | zero | 40 | 
| RD  | hi | sx_voxels  | none  | 89 | 
| RD  | hi | sx_voxels  | impute  | 89 |
| RD  | hi | sx_voxels  | zero | 89 | 

```bash
cd ~/data/ctsem_voxelwise/results/TI1/FA_inatt_perms
3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 25 \
    -quiet */voxels_sx.nii.gz | grep CLUSTERS | wc -l
```

Just careful because the code above gives us the number of times we DID NOT have
a cluster of 28 or more voxels. Then, we still need to do (nperms-res)/nperms.

So, a better way to get the p-value is:

```bash
cd ~/data/ctsem_voxelwise/results/TI1/FA_inatt_perms
res=`3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 25 \
    -quiet */voxels_sx.nii.gz | grep CLUSTERS | wc -l`
nperms=`ls | grep _p | wc -l`;
p=$(bc <<<"scale=3;($nperms - $res)/$nperms")
echo negatives=${res}, perms=${nperms}, pval=$p
```

We need more perms, but it looks like the result with 58 voxels is significant,
but the 25-28 result is not even close.

# 2019-07-24 10:08:26

Because the cluster went down, I'm compiling some of the results locally now. We
can do it in caterpie and ncrshell01. So, let's make sure it all works:

```bash
# desktop or philip
sx=inatt
m=AD
TI=TI1_perms
cd ~/data/ctsem_voxelwise/results/TI1/
ncpu=6;
for p in `seq 1 100`; do
    echo "Working on perm $p...";
    rm -f *.csv *.txt *.nii.gz *tgz vlist*;
    csv_dir=${TI}_${m}_${sx}_p${p};
    # create big file just once, breaking it up because initial list of arguments is too long
    for i in {1..9}; do
        echo $i;
        cat ${csv_dir}/*Y${i}*.csv >> mega_file.csv
    done
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h msg mega_file.csv > msgs.csv;
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        for v in `cat vlist_rerun.sx_sx`; do
            echo ${v} > ${v}.txt;
            echo ${v} >> var_list.txt;
        done
        pseed=`cut -d" " -f 8 swarm.${m}_${sx}_p${p} | head -n 1`;
        echo "Random seed: $pseed";
        # made sure that the script is setup for TI1 now!
        cat var_list.txt | parallel -j $ncpu --max-args=1 \
            Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R \
                ~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz \
                sx_${sx} {}.txt ${csv_dir}/redo_{}.csv $pseed;
    fi;

    # and assuming no other reruns are needed
    rm -f sx_sx*.nii.gz *.csv *.txt;
    for i in {1..9}; do
        echo $i;
        cat ${csv_dir}/*Y${i}*.csv >> mega_file.csv
    done
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h _sx_${sx}_Y mega_file.csv > sx_voxels.csv;
    grep -h "Y[0-9]\+_sx_${sx}" mega_file.csv > voxels_sx.csv;
    grep -h "Y[0-9]\+_Y[0-9]\+" mega_file.csv > voxel_voxel.csv;
    grep -h AIC mega_file.csv > aic.csv;
    grep -h BIC mega_file.csv > bic.csv;
    grep -h msg mega_file.csv > msgs.csv;
    for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
        echo $f perm $p;
        python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
            ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
            ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
    done
    # if compiled successfully
    if [ -e sx_voxels.nii.gz ]; then
        echo "Cleaning up...";
        mkdir -p ./clean/${csv_dir};
        # reducing the number of files a bit
        find ${csv_dir}/. -name '*.csv' -print > all_files.txt;
        tar -zcf ${m}_${sx}_p$p.tgz --files-from all_files.txt;
        cp -v *nii.gz ${m}_${sx}_p$p.tgz ./clean/${csv_dir}/;
    fi;
done
```

And the code above actually worked quite well, splitting the voxels in the
parallel input. I had also done something wrong incompiling the FA results, so
let's go ahead and re-run those, and then re-clean the old ones.

But I had to hack the old ones to make sure it had a .tgz with all the voxels in
it, like this:

```bash
for p in {1..100}; do
    if [ ! -e ~/data/ctsem_voxelwise/TI1_perms/FA_inatt_p${p}/FA_inatt_p${p}.tgz ]; then
        echo $p;
        pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.FA_inatt_p${p} | head -n 1`;
        mkdir csv;
        /bin/ls ~/data/ctsem_voxelwise/TI1_perms/FA_*_inatt*_${pseed}.tgz > perm_files.txt;
        for f in `cat perm_files.txt`; do
            tar -zxf $f -C csv/;
        done
        mkdir -p ~/data/ctsem_voxelwise/TI1_perms/FA_inatt_p${p}/;
        tar -zcf FA_inatt_p${p}.tgz csv/*;
        cp -v FA_inatt_p${p}.tgz ~/data/ctsem_voxelwise/TI1_perms/FA_inatt_p${p}/;
        for f in `cat perm_files.txt`; do
            rm -f $f;
        done;
        rm -rf csv FA_inatt_p${p}.tgz;
    fi;
done
```

Then, it's just a matter of running:

```bash
#bw
module load afni
module load R
module load python/3.7
sx=inatt
m=FA
TI=TI1_perms
cd /lscratch/$SLURM_JOBID
ncpu=$SLURM_CPUS_PER_TASK;
for p in `seq 1 100`; do
    clean_fname=~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/${m}_inatt_p${p}.tgz
    echo "Working on perm $p...";
    rm -rf csv/;
    tar -zxvf $clean_fname;
    csv_dir=csv/;
    rm -f *.csv *.txt *.nii.gz *tgz vlist* $csv_dir/output*.csv;
    # create big file just once, breaking it up because initial list of arguments is too long
    cat ${csv_dir}/*.csv > mega_file.csv
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h msg mega_file.csv > msgs.csv;
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        for v in `cat vlist_rerun.sx_sx`; do
            echo ${v} > ${v}.txt;
            echo ${v} >> var_list.txt;
        done
        pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.${m}_${sx}_p${p} | head -n 1`;
        # made sure that the script is setup for TI1 now!
        cat var_list.txt | parallel -j $ncpu --max-args=1 \
            Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R \
                ~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz \
                sx_${sx} {}.txt ${csv_dir}/redo_{}.csv $pseed;
    fi;

    # and assuming no other reruns are needed
    rm -f sx_sx*.nii.gz *.csv *.txt;
    cat ${csv_dir}/*.csv > mega_file.csv
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h _sx_${sx}_Y mega_file.csv > sx_voxels.csv;
    grep -h "Y[0-9]\+_sx_${sx}" mega_file.csv > voxels_sx.csv;
    grep -h "Y[0-9]\+_Y[0-9]\+" mega_file.csv > voxel_voxel.csv;
    grep -h AIC mega_file.csv > aic.csv;
    grep -h BIC mega_file.csv > bic.csv;
    grep -h msg mega_file.csv > msgs.csv;
    for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
        echo $f perm $p;
        python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
            ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
            ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
    done
    # if compiled successfully
    if [ -e sx_voxels.nii.gz ]; then
        echo "Cleaning up...";
        mkdir -p ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p};
        # reducing the number of files a bit
        tar -zcf ${m}_${sx}_p$p.tgz csv/*;
        cp -v *nii.gz ${m}_${sx}_p$p.tgz ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p}/;
    fi;
done
```

While we wait for some of these results to compile locally, let's make some
pictures of the FA results.

```bash
3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 50 \
    -overwrite -savemask mymask.nii results/TI1/FA_inatt/sx_voxels.nii.gz
fname=FA_inatt_sx_voxels_inHCP1065.nii.gz
flirt -in mymask.nii \
    -ref /usr/local/fsl/data/standard/FSL_HCP1065_FA_1mm.nii \
    -out $fname -applyxfm \
    -init group_skeleton_to_HCP1065.mat -interp nearestneighbour
# just to get the COM coordinates in the new space
3dclust -1Dformat -nosum -orient LPI -NN1 5 $fname
```

# 2019-07-29 15:46:38

Now that BW is back, we should just make a big swarm file with everything that
needs to be re-run for FA perms, and just swarm that, instead of doing them
locally in interactive machines!

```bash
#bw
module load afni
module load R
module load python/3.7

# just so we have a single place to refer to
mkdir /scratch/sudregp/ctsem/
for v in `seq 1 10814`; do
    echo Y${v} > /scratch/sudregp/ctsem/Y${v}.txt;
done

sx=inatt
m=FA
TI=TI1_perms
cd /lscratch/$SLURM_JOBID
sfile=~/data/ctsem_voxelwise/swarm.rerun5
for p in `seq 81 100`; do
    clean_fname=~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/${m}_inatt_p${p}.tgz
    echo "Working on perm $p...";
    rm -rf csv/;
    tar -zxf $clean_fname;
    csv_dir=csv/;
    rm -f *.csv *.txt vlist_rerun.sx_sx *tgz vlist* $csv_dir/output*.csv;
    # create big file just once, breaking it up because initial list of arguments is too long
    cat ${csv_dir}/*.csv > mega_file.csv
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h msg mega_file.csv > msgs.csv;
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.${m}_${sx}_p${p} | head -n 1`;
        for y in `cat vlist_rerun.sx_sx`; do
            echo "Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R" \
                "~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz" \
                "sx_${sx} /scratch/sudregp/ctsem/${y}.txt " \
                "~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/redo_${y}.csv" \
                "$pseed" >> $sfile;
        done
    fi;
done
```

Then, it's just a matter of how long we think this whole think will take, and
how much we think the cluster can handle. Since it's a one-time thing, I'll try
to run it all at once with a bundling that doesn't allow more than a few hundred
accesses to the network at a time.

```bash
#bw
cd ~/data/ctsem_voxelwise/
swarm -f swarm.rerun -t 2 --module R,afni,python/3.7 --time=1:00:00 \
    --merge-output --logdir=trash_FAcomp --job-name FAcomp --partition quick
```

# 2019-07-30 09:49:05

Now we can recompile the FA results:

```bash
#bw
module load afni
module load R
module load python/3.7
sx=inatt
m=FA
TI=TI1_perms
cd /lscratch/$SLURM_JOBID
ncpu=$SLURM_CPUS_PER_TASK;
for p in `seq 1 100`; do
    clean_fname=~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/${m}_inatt_p${p}.tgz
    echo "Working on perm $p...";
    rm -rf csv/;
    tar -zxf $clean_fname;
    csv_dir=csv/;
    rm -f *.csv *.txt *.nii.gz *tgz vlist* $csv_dir/output*.csv vlist_rerun.sx_sx;
    cp ~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/redo_*.csv csv/;
    # create big file just once, breaking it up because initial list of arguments is too long
    cat ${csv_dir}/*.csv > mega_file.csv
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h msg mega_file.csv > msgs.csv;
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        for v in `cat vlist_rerun.sx_sx`; do
            echo ${v} > ${v}.txt;
            echo ${v} >> var_list.txt;
        done
        pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.${m}_${sx}_p${p} | head -n 1`;
        echo "Random seed: $pseed";
        # made sure that the script is setup for TI1 now!
        cat var_list.txt | parallel -j $ncpu --max-args=1 \
            Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R \
                ~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz \
                sx_${sx} {}.txt ${csv_dir}/redo_{}.csv $pseed;
    fi;

    # and assuming no other reruns are needed
    rm -f sx_sx*.nii.gz *.csv *.txt;
    cat ${csv_dir}/*.csv > mega_file.csv
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h _sx_${sx}_Y mega_file.csv > sx_voxels.csv;
    grep -h "Y[0-9]\+_sx_${sx}" mega_file.csv > voxels_sx.csv;
    grep -h "Y[0-9]\+_Y[0-9]\+" mega_file.csv > voxel_voxel.csv;
    grep -h AIC mega_file.csv > aic.csv;
    grep -h BIC mega_file.csv > bic.csv;
    grep -h msg mega_file.csv > msgs.csv;
    for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
        echo $f perm $p;
        python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
            ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
            ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
    done
    # if compiled successfully
    if [ -e sx_voxels.nii.gz ]; then
        echo "Cleaning up...";
        mkdir -p ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p};
        # reducing the number of files a bit
        tar -zcf ${m}_${sx}_p$p.tgz csv/*;
        cp -v *nii.gz ${m}_${sx}_p$p.tgz ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p}/;
    fi;
done
```

OK, this is getting annoying because I have lots of processors halted because of
a single voxel. Let's redo the idea of using a swarm file instead:

```bash
#bw
module load afni
module load R
module load python/3.7
sx=inatt
m=FA
TI=TI1_perms
cd /lscratch/$SLURM_JOBID
sfile=~/data/ctsem_voxelwise/swarm.FArerun
for p in `seq 51 100`; do
    clean_fname=~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/${m}_inatt_p${p}.tgz
    echo "Working on perm $p...";
    rm -rf csv/;
    tar -zxf $clean_fname;
    csv_dir=csv/;
    rm -f *.csv *.txt *.nii.gz *tgz vlist* $csv_dir/output*.csv vlist_rerun.sx_sx;
    cp ~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/redo_*.csv csv/;
    # create big file just once, breaking it up because initial list of arguments is too long
    cat ${csv_dir}/*.csv > mega_file.csv
    grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
    grep -h msg mega_file.csv > msgs.csv;
    python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
        sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
        ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

    # re-run whatever needs to be re-run
    if [ -e vlist_rerun.sx_sx ]; then
        pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.${m}_${sx}_p${p} | head -n 1`;
        echo "Random seed: $pseed";
        for y in `cat vlist_rerun.sx_sx`; do
            echo "Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R" \
                "~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz" \
                "sx_${sx} /scratch/sudregp/ctsem/${y}.txt " \
                "~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/redo_${y}.csv" \
                "$pseed" >> $sfile;
        done
    else
        # and assuming no other reruns are needed
        rm -f sx_sx*.nii.gz *.csv *.txt;
        cat ${csv_dir}/*.csv > mega_file.csv
        grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
        grep -h _sx_${sx}_Y mega_file.csv > sx_voxels.csv;
        grep -h "Y[0-9]\+_sx_${sx}" mega_file.csv > voxels_sx.csv;
        grep -h "Y[0-9]\+_Y[0-9]\+" mega_file.csv > voxel_voxel.csv;
        grep -h AIC mega_file.csv > aic.csv;
        grep -h BIC mega_file.csv > bic.csv;
        grep -h msg mega_file.csv > msgs.csv;
        for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
            echo $f perm $p;
            python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
                ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
                ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
        done
        # if compiled successfully
        if [ -e sx_voxels.nii.gz ]; then
            echo "Cleaning up...";
            mkdir -p ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p};
            # reducing the number of files a bit
            tar -zcf ${m}_${sx}_p$p.tgz csv/*;
            cp -v *nii.gz ${m}_${sx}_p$p.tgz ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p}/;
        fi;
    fi;
done
```

# 2019-08-01 15:17:58

Let's recompile the ones that were re-run. If they still aren't done, it's
usually because of a single pesky voxel that didn't finish even after running it
twice, so I removed the set.seed() from the R script to see if it would help.
Those voxels are likely crap anyways, so I'm not too worried.

```bash
#bw
module load afni
module load R
module load python/3.7
sx=inatt
m=FA
TI=TI1_perms
cd /lscratch/$SLURM_JOBID
sfile=~/data/ctsem_voxelwise/swarm.FArerun
rm $sfile
for p in `seq 1 100`; do
    if [ ! -d ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p} ]; then
        clean_fname=~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/${m}_inatt_p${p}.tgz
        echo "Working on perm $p...";
        rm -rf csv/;
        tar -zxf $clean_fname;
        csv_dir=csv/;
        rm -f *.csv *.txt *.nii.gz *tgz vlist* $csv_dir/output*.csv vlist_rerun.sx_sx;
        cp ~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/redo_*.csv csv/;
        # create big file just once, breaking it up because initial list of arguments is too long
        cat ${csv_dir}/*.csv > mega_file.csv
        grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
        grep -h msg mega_file.csv > msgs.csv;
        python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
            sx_sx.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
            ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv

        # re-run whatever needs to be re-run
        if [ -e vlist_rerun.sx_sx ]; then
            pseed=`cut -d" " -f 8 ~/data/ctsem_voxelwise/swarm.${m}_${sx}_p${p} | head -n 1`;
            echo "Random seed: $pseed";
            for y in `cat vlist_rerun.sx_sx`; do
                echo "Rscript ~/research_code/ctsem_voxel_developmental_time_3_timepoints.R" \
                    "~/data/ctsem_voxelwise/${m}_wider_3obs_developmental_time.RData.gz" \
                    "sx_${sx} /scratch/sudregp/ctsem/${y}.txt " \
                    "~/data/ctsem_voxelwise/TI1_perms/${m}_inatt_p${p}/redo_${y}.csv" \
                    "$pseed" >> $sfile;
            done
        else
            # and assuming no other reruns are needed
            rm -f sx_sx*.nii.gz *.csv *.txt;
            cat ${csv_dir}/*.csv > mega_file.csv
            grep -h sx_${sx}_sx_${sx} mega_file.csv > sx_sx.csv;
            grep -h _sx_${sx}_Y mega_file.csv > sx_voxels.csv;
            grep -h "Y[0-9]\+_sx_${sx}" mega_file.csv > voxels_sx.csv;
            grep -h "Y[0-9]\+_Y[0-9]\+" mega_file.csv > voxel_voxel.csv;
            grep -h AIC mega_file.csv > aic.csv;
            grep -h BIC mega_file.csv > bic.csv;
            grep -h msg mega_file.csv > msgs.csv;
            for f in sx_sx sx_voxels voxels_sx voxel_voxel aic bic; do
                echo $f perm $p;
                python3 ~/research_code/fmri/compile_ctsem_voxel_results.py \
                    ${f}.csv ~/data/ctsem_voxelwise/mean_FA_skeleton_mask.nii.gz \
                    ~/data/ctsem_voxelwise/ctsem_ijk.txt msgs.csv;
            done
            # if compiled successfully
            if [ -e sx_voxels.nii.gz ]; then
                echo "Cleaning up...";
                mkdir -p ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p};
                # reducing the number of files a bit
                tar -zcf ${m}_${sx}_p$p.tgz csv/*;
                cp -v *nii.gz ${m}_${sx}_p$p.tgz ~/data/ctsem_voxelwise/TI1_perms/clean/${m}_${sx}_p${p}/;
            fi;
        fi;
    fi;
done
```

# TODO
 * finalize running and ocmpiling perms
 * visualize results
