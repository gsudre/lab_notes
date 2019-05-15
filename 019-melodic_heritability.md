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






```bash
#bw
cd /data/NCR_SBRB/software/Yeo_JNeurophysiol11_MNI152/
@SSwarper -input FSL_MNI152_FreeSurferConformed_1mm.nii.gz -base TT_N27_SSW.nii.gz -subid FSL_MNI152_inTLRC
3dNwarpApply -nwarp "anatQQ.FSL_MNI152_inTLRC_WARP.nii anatQQ.FSL_MNI152_inTLRC.aff12.1D" \
        -source Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz \
        -master ~/data/baseline_prediction/same_space/epi/0951_epi_NL_inTLRC.nii -inter NN \
        -overwrite -prefix ~/data/baseline_prediction/same_space/epi/Yeo2011_7Networks_LiberalMask_NL_inTLRC.nii;
cd ~/data/baseline_prediction/same_space/epi/
for i in {1..7}; do
    3dcalc -prefix Yeo_liberal_inTLRC_net${i}.nii \
        -a Yeo2011_7Networks_LiberalMask_NL_inTLRC.nii -expr "amongst(a,${i})";
done
3dTcat -prefix Yeo_liberal_inTLRC_combined.nii Yeo_liberal_inTLRC_net1.nii \
    Yeo_liberal_inTLRC_net2.nii Yeo_liberal_inTLRC_net3.nii \
    Yeo_liberal_inTLRC_net4.nii Yeo_liberal_inTLRC_net5.nii \
    Yeo_liberal_inTLRC_net6.nii Yeo_liberal_inTLRC_net7.nii
```

```bash
#bw
cd ~/data/baseline_prediction/same_space/epi/groupmelodic_inter.ica/
3dMatch -inset melodic_IC.nii.gz -refset ../Yeo_liberal_inTLRC_combined.nii \
    -mask ../group_epi_mask_inter.nii -prefix matches -overwrite
cat matches_REF_coeff.vals
```

```bash
pipe=fancy;
cd ~/data/baseline_prediction/same_space/epi/groupmelodic_${pipe}.ica/dual
for maskid in `cat ../../3min_clean.txt`; do
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in 37 24 57 4 54 1 2; do  # fancy
    # for i in 17 8 7 11 31 12 2; do  # inter
    # for i in 42 23 65 8 82 3 4; do  # union
        3dmaskdump -mask ../../group_epi_mask_${pipe}.nii \
            -o dumps/${maskid}_IC${i}_Z.txt dr_stage2_${maskid}_Z.nii.gz[${i}];
    done;
done
```

```r
library(gdata)
# assuming all fmri scans also have a mprage
clin = read.xls('~/data/baseline_prediction/long_scans_08072018.xlsx', 'mprage')
clin = clin[, c(1,4)]
colnames(clin) = c('MRN', 'mask.id')
pipes = c('inter', 'fancy', 'union')
ics = list(c(17, 8, 7, 11, 31, 12, 2), c(37, 24, 57, 4, 54, 1, 2), c(42, 23, 65, 8, 82, 3, 4)) # same order as pipes!
pheno = read.table('~/data/baseline_prediction/same_space/epi/3min_clean.txt')[, 1]
for (n in 1:length(pipes)) {
    mydir = sprintf('~/data/baseline_prediction/same_space/epi/groupmelodic_%s.ica/dual/dumps/', pipes[n])
    junk = read.table(sprintf('%s/0691_IC%d_Z.txt', mydir, ics[[n]][1]))
    nvox = nrow(junk)
    for (m in ics[[n]]) {
        fmri_data = matrix(nrow=length(pheno), ncol=nvox)
        for (s in 1:nrow(fmri_data)) {
            print(sprintf('%s %d %04d', pipes[n], m, pheno[s]))
            fname = sprintf('%s/%04d_IC%d_Z.txt', mydir, pheno[s], m)
            a = read.table(fname)
            fmri_data[s, ] = a[, 4]
        }
        data = cbind(pheno, fmri_data)
        cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
        colnames(data) = cnames
        data = merge(clin, data, by='mask.id')
        save(data, file=sprintf('~/data/baseline_prediction/melodic_%s_IC%d_12142018.RData.gz', pipes[n], m),
             compress=T)
    }
}
```
