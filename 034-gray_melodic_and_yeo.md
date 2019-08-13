# 2019-08-12 16:44:41

It's possible that my MELODIC results are not so good because my mask is taking
into considerationt oo much stuff that's not grey matter. So, let's repeat the
results using a better grey matter mask, and see how they look. I'll also redo
the yeo results just for kicks, but stick to FD1.0 for now.

And, since we're at it, why not compute everything for Z and noZ?

```bash
cd ~/data/heritability_change/xcp-36p_despike
melodic -i fd1_epi.txt -o groupmelodic_gray.ica -v --nobet -m gray_matter_mask.nii --tr=2.5 --report --Oall -a concat;
```

# 2019-08-13 09:32:14

Now we perform the dual regression to get each subject's values for the ICs:

```bash
pipe='gray';
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_${pipe}.ica
mkdir dual
while read maskid; do
    m=`printf %04d $maskid`;
    echo ${pipe} $m;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d melodic_IC \
        -o dual/dr_stage1_${m}.txt --demean -m ../gray_matter_mask.nii;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d dual/dr_stage1_${m}.txt \
        -o dual/dr_stage2_${m} --demean -m ../gray_matter_mask.nii --des_norm \
        --out_z=dual/dr_stage2_${m}_Z;
done < ../ids_1.txt
```

Now, it's time to figure out which ICs we are going to use.

```bash
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_gray.ica/
3dMatch -inset melodic_IC.nii.gz -refset ../Yeo_nets.nii \
    -mask ../gray_matter_mask.nii -prefix matches -overwrite
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

```
0               1               0.367           0.214
1               33              0.331           0.214
2               10              0.283           0.165
3               6               0.372           0.147
4               27              0.368           0.109
5               9               0.264           0.212
6               8               0.300           0.310
```

DMN doesn't look great, but let's go with it for now. Since we're doing Yeo
masks as well, I won't be too picky for now.

```bash
cd ~/data/heritability_change/xcp-36p_despike/yeo_masks_gray
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
mkdir dual
while read maskid; do
    m=`printf %04d $maskid`;
    echo yeo_masks_gray $m;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d ../Yeo_nets.nii \
        -o dual/dr_stage1_${m}.txt --demean -m ../gray_matter_mask.nii;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d dual/dr_stage1_${m}.txt \
        -o dual/dr_stage2_${m} --demean -m ../gray_matter_mask.nii --des_norm \
        --out_z=dual/dr_stage2_${m}_Z;
done < ../ids_1.txt
```

Time to dump to R:

```bash
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_gray.ica/
mkdir dumps
for m in `cat ../ids_1.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in 1 33 10 6 27 9 8; do
        3dmaskdump -mask ../gray_matter_mask.nii \
            -o dumps/${maskid}_IC${i}_Z.txt dual/dr_stage2_${maskid}_Z.nii.gz[${i}];
        3dmaskdump -mask ../gray_matter_mask.nii \
            -o dumps/${maskid}_IC${i}.txt dual/dr_stage2_${maskid}.nii.gz[${i}];
    done;
done
```

<!-- ```bash
cd ~/data/heritability_change/xcp-36p_despike/yeo_masks_gray/
mkdir dumps
for m in `cat ../ids_1.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in {0..6}; do
        3dmaskdump -mask ../gray_matter_mask.nii \
            -o dumps/${maskid}_net${i}_Z.txt dual/dr_stage2_${maskid}_Z.nii.gz[${i}];
        3dmaskdump -mask ../gray_matter_mask.nii \
            -o dumps/${maskid}_net${i}.txt dual/dr_stage2_${maskid}.nii.gz[${i}];
    done;
done
``` -->

Then, we collect our results in R:

<!-- ```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=155301
for (m in c(1, 33, 10, 6, 27, 9, 8)) {
    for (s in c('', '_Z')) {
        print(m)
        print(s)
        brain_data = matrix(nrow=length(maskids), ncol=nvox)
        for (s in 1:nrow(brain_data)) {
            fname = sprintf('~/data/heritability_change/xcp-36p_despike/groupmelodic_gray.ica/dumps/%04d_IC%d%s.txt', maskids[s], m, s)
            a = read.table(fname)
            brain_data[s, ] = a[,4]
        }
        brain_data = cbind(maskids, brain_data)
        cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
        colnames(brain_data) = cnames
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/melodic_gray_IC%d%s.rds', m, s)
        saveRDS(brain_data, file=fname)
    }
}
``` -->

Then, repeat the same for the yeo_masks dual regression.

<!-- ```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=155301
for (m in 0:6) {
    for (s in c('', '_Z')) {
        print(m)
        print(s)
        brain_data = matrix(nrow=length(maskids), ncol=nvox)
        for (s in 1:nrow(brain_data)) {
            fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks_gray/dumps/%04d_net%d%s.txt', maskids[s], m, s)
            a = read.table(fname)
            brain_data[s, ] = a[,4]
        }
        brain_data = cbind(maskids, brain_data)
        cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
        colnames(brain_data) = cnames
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks_gray_net%d%s.rds', m, s)
        saveRDS(brain_data, file=fname)
    }
}
``` -->

