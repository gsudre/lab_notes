# 2019-05-14 15:55:57

OK, let's see if there is anything heritable in the MELODIC ICs. We start by
defining some masks:

```bash
# caterpie
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/
cut -d"," -f 1 fmri_corr_tables/pearson_3min_n462_power.csv | tail -n +2 > 3min_mni.txt
for maskid in `cat 3min_mni.txt`; do
    m=`printf %04d $maskid`;
    3dAutomask -prefix masks/${m}_automask.nii fmri_same_space/epi/${m}_epi_NL_inMNI.nii;
done
cd masks
3dmask_tool -input ????_automask.nii -prefix ../group_epi_mask_union.nii -frac 0
3dmask_tool -input ????_automask.nii -prefix ../group_epi_mask_inter.nii -frac 1
3dmask_tool -input ????_automask.nii -prefix ../group_epi_mask_fancy.nii \
    -dilate_input 5 -5 -frac 0.7 -fill_holes
```

# 2019-05-15 09:21:19

Let's then run melodic. But we'll need to send the data to BW first:

```bash
# caterpie
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/
for maskid in `cat 3min_mni.txt`; do
    m=`printf %04d $maskid`;
    echo ${m}_epi_NL_inMNI.nii >> fmri_same_space/epi/3min_mni_epi.txt;
done
scp fmri_same_space/epi/3min_mni_epi.txt bw:~/data/heritability_change/fmri_same_space/epi/;

# bw
module load fsl/6.0.0
cd ~/data/heritability_change/fmri_same_space/epi/;
melodic -i 3min_mni_epi.txt -o groupmelodic_union.ica -v --nobet -m ../group_epi_mask_union.nii --tr=2.5 --report --Oall -a concat;
melodic -i 3min_mni_epi.txt -o groupmelodic_fancy.ica -v --nobet -m ../group_epi_mask_fancy.nii --tr=2.5 --report --Oall -a concat;
melodic -i 3min_mni_epi.txt -o groupmelodic_inter.ica -v --nobet -m ../group_epi_mask_inter.nii --tr=2.5 --report --Oall -a concat;
```

Now we performt the dual regression to get each subject's values for the ICs:

```bash
pipe='inter';
cd ~/data/heritability_change/fmri_same_space/epi/groupmelodic_${pipe}.ica
mkdir dual
while read m; do
    s=`printf %04d $m`;
    echo ${pipe} $s;
    $FSLDIR/bin/fsl_glm -i ../${s}_epi_NL_inMNI.nii -d melodic_IC \
        -o dual/dr_stage1_${s}.txt --demean -m ../../group_epi_mask_${pipe}.nii;
    $FSLDIR/bin/fsl_glm -i ../${s}_epi_NL_inMNI.nii -d dual/dr_stage1_${s}.txt \
        -o dual/dr_stage2_$s --demean -m ../../group_epi_mask_${pipe}.nii --des_norm \
        --out_z=dual/dr_stage2_${s}_Z;
done < ../../../3min_mni.txt
```

# 2019-05-16 10:57:18

We are already in MNI space, like the Yeo networks. So, let's go with that:

```bash
#bw
cd ~/data/heritability_change/fmri_same_space/epi/;
for i in {1..7}; do
    3dcalc -prefix Yeo_liberal_inMNI_net${i}.nii \
        -a /data/NCR_SBRB/software/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz -expr "amongst(a,${i})";
done
3dTcat -prefix Yeo_liberal_inMNI_combined.nii Yeo_liberal_inMNI_net1.nii \
    Yeo_liberal_inMNI_net2.nii Yeo_liberal_inMNI_net3.nii \
    Yeo_liberal_inMNI_net4.nii Yeo_liberal_inMNI_net5.nii \
    Yeo_liberal_inMNI_net6.nii Yeo_liberal_inMNI_net7.nii
3dresample -master groupmelodic_inter.ica/melodic_IC.nii.gz \
    -prefix Yeo_nets.nii -inset Yeo_liberal_inMNI_combined.nii \
    -rmode NN -overwrite
```

So, let's figure out what are the best matching networks for each mask:

```bash
#bw
cd ~/data/heritability_change/fmri_same_space/epi/groupmelodic_inter.ica/
3dMatch -inset melodic_IC.nii.gz -refset ../Yeo_nets.nii \
    -mask ../../group_epi_mask_inter.nii -prefix matches -overwrite
cat matches_REF_coeff.vals
```

Keep in mind that the code is:

```
0: visual
1: somatomotor
2: DAN
3: VAN
4: limbic
5: cognitive (frontoparietal)
6: DMN
```

inter:

```
0               41              0.173           0.251
1               48              0.240           0.144
2               35              0.325           0.083
3               23              0.214           0.170
4               214             0.261           0.025
5               228             0.174           0.191
6               0               0.179           0.296
```

I didn't like the DMN component in inter.

fancy:

```
0               60              0.364           0.138
1               15              0.344           0.122
2               2               0.381           0.097
3               8               0.365           0.086
4               62              0.413           0.074
5               5               0.350           0.122
6               0               0.322           0.181
```

The one from fancy looks much better:

![](images/2019-05-16-11-31-29.png)

And this is cognitive:

![](images/2019-05-16-11-43-11.png)

Let's go with fancy for now.

union:

```
0               107             0.340           0.054
1               138             0.295           0.048
2               70              0.310           0.038
3               34              0.310           0.033
4               67              0.405           0.029
5               98              0.256           0.048
6               48              0.267           0.073
```

Union mask has too much crap...

```bash
pipe=fancy;
cd ~/data/heritability_change/fmri_same_space/epi/groupmelodic_${pipe}.ica/dual
mkdir dumps
for m in `cat ../../../../3min_mni.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in 60 15 2 8 62 5 0; do  # fancy
        3dmaskdump -mask ../../../group_epi_mask_${pipe}.nii \
            -o dumps/${maskid}_IC${i}_Z.txt dr_stage2_${maskid}_Z.nii.gz[${i}];
    done;
done
```

Then, we collect our results in R:

```r
maskids = read.table('~/data/heritability_change/3min_mni.txt')[, 1]
nvox=154058
for (m in c(60, 15, 2, 8, 62, 5, 0)) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/fmri_same_space/epi/groupmelodic_fancy.ica/dual/dumps/%04d_IC%d_Z.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/fmri_same_space/melodic_fancy_IC%d.RData', m)
     save(brain_data, file=fname)
}
```

Now that the data is into CSV, we need to assign MRNs, calculate slopes, and
prepare it for SOLAR voxelwise.

**Note that I won' be removing movement here! Mostly because since we're using
ICA, the movement components should have been isolated already. But we can
always check any results later against correlation to movement.
**

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
m2 = read.csv('~/data/heritability_change/rsfmri_3min_assoc_n462.csv')
clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
load('~/data/heritability_change/fmri_same_space/melodic_fancy_IC0.RData')
b = brain_data
var_names = colnames(b)[2:ncol(b)]
df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id', all.x=F)

# make sure we still have two scans for everyone
rm_subjs = names(which(table(df2$Medical.Record...MRN)<2))
rm_me = df2$Medical.Record...MRN %in% rm_subjs
df2 = df2[!rm_me, ]

mres = df2
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))

res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    y = mres[idx[2], var_names] - mres[idx[1], var_names]
    x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
    slopes = y / x
    row = c(row, slopes)
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], 'SX_inatt'])
    row = c(row, mres[idx[base_DOA], 'SX_HI'])
    # DX1 is DSMV definition, DX2 will make SX >=4 as ADHD
    if (mres[idx[base_DOA], 'age_at_scan'] < 16) {
        if ((row[length(row)] >= 6) || (row[length(row)-1] >= 6)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    } else {
        if ((row[length(row)] >= 5) || (row[length(row)-1] >= 5)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    }
    if ((row[length(row)] >= 4) || (row[length(row)-1] >= 4)) {
        DX2 = 'ADHD'
    } else {
        DX2 = 'NV'
    }
    row = c(row, DX)
    row = c(row, DX2)
    res = rbind(res, row)
    print(nrow(res))
}
colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI',
                                              'inatt_baseline',
                                              'HI_baseline', 'DX', 'DX2'))
# we only open this in R, so it's OK to be RData to load faster
save(res, file='~/data/heritability_change/fmri_same_space/melodic_fancy_slopes_IC0.RData')

# and remove outliers
res_clean = res
for (t in var_names) {
    mydata = as.numeric(res_clean[, t])
    # identifying outliers
    ul = mean(mydata) + 3 * sd(mydata)
    ll = mean(mydata) - 3 * sd(mydata)
    bad_subjs = c(which(mydata < ll), which(mydata > ul))

    # remove within-variable outliers
    res_clean[bad_subjs, t] = NA
}
save(res_clean, file='~/data/heritability_change/fmri_same_space/melodic_fancy_slopesClean_IC0.RData')

# and make sure every family has at least two people
good_nuclear = names(table(m2$Nuclear.ID...FamilyIDs))[table(m2$Nuclear.ID...FamilyIDs) >= 4]
good_extended = names(table(m2$Extended.ID...FamilyIDs))[table(m2$Extended.ID...FamilyIDs) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, m2[which(m2$Nuclear.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
for (f in good_extended) {
    keep_me = c(keep_me, m2[which(m2$Extended.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
keep_me = unique(keep_me)

fam_subjs = c()
for (s in keep_me) {
    fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
}
res2 = res[fam_subjs, ]
res2_clean = res_clean[fam_subjs, ]

write.csv(res2, file='~/data/heritability_change/fmri_same_space/melodic_fancy_slopes_n111_IC0.csv', row.names=F, na='', quote=F)
write.csv(res2_clean, file='~/data/heritability_change/fmri_same_space/melodic_fancy_slopesClean_n111_IC0.csv', row.names=F, na='', quote=F)
```

And of course, redo all of the above for all 6 ICs.

# 2019-05-17 10:14:19

Now that we have all files, let's set it up to run voxelwise SOLAR:

```bash
cd ~/data/heritability_change/fmri_same_space/
jname=ic0Clean
fname=${jname}.swarm;
for i in {1..154058}; do
    echo "bash ~/research_code/run_solar_voxel.sh melodic_fancy_slopesClean_n111_IC0 ${i}" >> $fname;
done;
swarm --gres=lscratch:1 -f ${fname} --module solar -g 1 -t 1 \
            --logdir=${jname} --job-name ${jname} -p 2 --partition quick \
            --time=1 -b 240;
```

And of course, do the same for 60, 15, 2, 8, 62, and 5.

Time to compile the voxel results. First, create an ijk file:

```bash
cd ~/data/heritability_change/fmri_same_space/
cut -d " " -f 1,2,3 \
    epi/groupmelodic_fancy.ica/dual/dumps/0901_IC0_Z.txt > ../group_mask_fancy_ijk.txt
```

We compile using:

```bash
python ~/research_code/fmri/compile_solar_voxel_results.py melodic_fancy_slopesClean_n111_IC0
3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 15 \
    ~/data/tmp/polygen_results_melodic_fancy_slopesClean_n111_IC0.nii
```

And since we have a weekend tomorrow, let's leave some of the permutations
running. At least for DMN and FP. On Monday we can check which networks actually
have clusters that are also associated with ADHD:

```bash
cd /lscratch/$SLURM_JOBID
cp ~/data/heritability_change/melodic_fancy_slopesClean_n111_IC0.csv .
```

```r
# start it from lscratch
m = 0
nperms = 250
library(data.table)
dread = fread(sprintf('~/data/heritability_change/melodic_fancy_slopesClean_n111_IC%d.csv', m), header = T, sep = ',')
d = as.data.frame(dread)  # just so we can index them a bit easier
vcols = c(which(grepl("v",colnames(d))), which(grepl("sex",colnames(d))))
d2 = d
for (p in 1:nperms) {
    d2[, vcols] = d[sample(nrow(d)), vcols]
    fname = sprintf('~/data/heritability_change/perms/melodic_fancy_slopesClean_n111_IC%d_p%04d_sexAndBrainShuffled.csv', m, p)
    print(fname)
    fwrite(d2, file=fname, row.names=F, quote=F)
}
```

Unfortunately that's taking too long. Maybe do something like this? 
https://howto.lintel.in/shuffle-lines-file-linux/ using shuf()?


And we need to be careful to only submit the swarms when we can:

```bash
cd ~/data/heritability_change/fmri_same_space/
ic=0;
nperms=250;
fname=`printf net%d_p%04d $ic 0`.swarm;
# generate voxel file for the first perm
for i in {1..154058}; do
    echo "bash ~/research_code/run_solar_voxel.sh melodic_fancy_slopesClean_n111_IC0_p0000_sexAndBrainShuffled ${i}" >> $fname;
done;
# just copy it and rename IC and perm for the other ones
for n in `seq 1 $nperms`; do
    perm=`printf %04d $n`;
    cp $fname `printf net%d_p%04d $ic $n`.swarm;
    sed -i -- "s/p0000/p${perm}/g" `printf net%d_p%04d $ic $n`.swarm;
done

# runs all swarms, but wait until we can do it
for n in `seq 1 $nperms`; do
    jname=`printf net%d_p%04d $ic $n`;
    echo "ERROR" > swarm_wait_${ic}
    while grep -q ERROR swarm_wait_${ic}; do
        echo "Trying $jname"
        swarm --gres=lscratch:1 -f ${jname}.swarm --module solar -g 1 -t 1 \
            --logdir=${jname} --job-name ${jname} -p 2 --partition quick \
            --time=1 -b 240 2> swarm_wait_${ic};
        if grep -q ERROR swarm_wait_${ic}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```








<!-- OK, let's see if any of these results have some association to ADHD. Wrote
individual masks in AFNI. Then:

```bash
cd /mnt/shaw/dti_robust_tsa/heritability
for n in 1 2 3; do
    for m in fa ad rd; do
        echo $m $n
        3dcopy ${m}_n${n}_mask_0001+orig mymask.nii -overwrite;
        echo maskid,val > ${m}_n${n}.csv;
        while read s; do
            val1=`3dmaskave -q -mask mymask.nii ../analysis_may2017/${s}_tensor_diffeo_${m}.nii.gz 2>/dev/null`;
            echo ${s},${val1} >> ${m}_n${n}.csv;
        done < maskids_566.txt;
    done;
done
``` -->


# TODO

* If results are good, make sure there is no correlation between clusters and movement!