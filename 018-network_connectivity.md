# 2019-05-09 15:00:50

We were getting hammered by multiple comparisons in rsFMRI, so let's instead
check for within and acorss network connectivity. We already have all subjects
in the common space, so now we just need to mark the different areas of each
network. We could do this by voxels or by ROI. And that's assuming we want a
within-network metric, otherwise we could just take a network to network
approach.

In any case, we'll need the data ll in the same space. First, let make sure we
have all the anatomical files we need:

```bash
net_dir=/Volumes/Shaw/MR_data_by_maskid/
cd ~/data/heritability_change/fmri_same_space/anat
export OMP_NUM_THREADS=4
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    if [ ! -e anatQQ.$m.nii ]; then
        echo $m;
        mri_convert /Volumes/Shaw/freesurfer5.3_subjects/${m}/mri/orig/001.mgz ./${m}.nii &&
        @SSwarper \
            -input ./${m}.nii \
            -base TT_N27_SSW.nii.gz -subid ${m} &&
        rm ${m}.nii;
    fi;
done;
```

Now we just do:

```bash
cd ~/data/heritability_change/fmri_same_space/epi;
for m in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m2=`printf %04d $m`;
    echo $m2;
    afni_dir=/Volumes/Shaw/MR_data_by_maskid/${m2}/afni/${m2}.rest.subjectSpace.results/;
    3dNwarpApply -nwarp "../anat/anatQQ.${m2}_WARP.nii ../anat/anatQQ.${m2}.aff12.1D" \
        -source ${afni_dir}/errts.${m2}.fanaticor+orig.HEAD \
        -master ../anat/anatQQ.${m2}.nii -dxyz 2.5\
        -overwrite -prefix ${m2}_epi_NL_inTLRC.nii;
done
```

But while that's going on I made sure all the transform QC picture looked
good...



    @SSwarper -input ${m2}.nii.gz -base TT_N27_SSW.nii.gz -subid ${m2};
    3dcalc -a /Volumes/Shaw/freesurfer5.3_subjects/${m2}/SUMA/aparc+aseg_REN_gm.nii.gz \
        -prefix ${m2}_gm_mask.nii -expr "step(a)";
    3dNwarpApply -nwarp "anatQQ.${m2}_WARP.nii anatQQ.${m2}.aff12.1D" \
        -source ${m2}_gm_mask.nii -master anatQQ.${m2}.nii -inter NN \
        -overwrite -prefix ${m2}_gm_mask_NL_inTLRC.nii;
done;

cd ~/data/baseline_prediction/same_space/epi;
for m2 in `cat ../../rsfmri_3minWithClinical.tsv`; do
    echo $m2;
    3dcopy ${mylink}/afni/${m2}.rest.subjectSpace.results/errts.${m2}.fanaticor+orig.HEAD ${m2}_epi.nii;
    3dNwarpApply -nwarp "../anat/anatQQ.${m}_WARP.nii ../anat/anatQQ.${m}.aff12.1D" \
        -source ${mylink}/afni/${m2}.rest.subjectSpace.results/errts.${m2}.fanaticor+orig.HEAD -master ../anat/anatQQ.${m}.nii -dxyz 2.5\
        -overwrite -prefix ${m}_epi_NL_inTLRC.nii;
done;
```

While the transformations were going on, I made sure to check the output of
SSwarper, which looked fine.

