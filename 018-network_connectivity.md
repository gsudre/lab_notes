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

But while that's going on I made sure all the transform QC picture looked
good...

# 2019-05-10 11:12:17

So, the anatomical warping was taking a very long time to run in the interactive
nodes. We'll need to swarm it...

```bash
cd ~/data/heritability_change/fmri_same_space/anat;
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    if [ ! -e anatQQ.$m.nii ]; then
        echo "export OMP_NUM_THREADS=16; cd ~/data/heritability_change/fmri_same_space/anat; @SSwarper -input ${m}.nii -base TT_N27_SSW.nii.gz -subid ${m}" >> swarm.SSwarper;
    fi;
done;
swarm -g 10 -t 16 --job-name SSwarper --time 2:00:00 -f swarm.SSwarper -m afni --partition quick --logdir trash
```

However, it makes more sense ot run all of this in Luke's framework, but it
needs them converted to MNI space instead. Let's compute that, shall we?

```bash
cd ~/data/heritability_change/fmri_same_space/anat_mni;
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    # mri_convert /data/NCR_SBRB/freesurfer5.3_subjects/${m}/mri/orig.mgz ${m}.nii.gz;
    echo "export OMP_NUM_THREADS=16; cd ~/data/heritability_change/fmri_same_space/anat_mni; @SSwarper -input ${m}.nii.gz -base MNI152_2009_template_SSW.nii.gz -subid ${m}" >> swarm.SSwarper;
done;
swarm -g 10 -t 16 --job-name SSwarper --time 2:00:00 -f swarm.SSwarper -m afni --partition quick --logdir trash
```

And we need to check on the alignments again. But first, transfer the actual EPI
signal to MNI space:

Now we just do:

```bash
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/fmri_same_space/epi;
for m in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m2=`printf %04d $m`;
    echo $m2;
    afni_dir=/mnt/shaw/MR_data_by_maskid/${m2}/afni/${m2}.rest.subjectSpace.results/;
    3dNwarpApply -nwarp "../anat_mni/anatQQ.${m2}_WARP.nii ../anat_mni/anatQQ.${m2}.aff12.1D" \
        -source ${afni_dir}/errts.${m2}.fanaticor+orig.HEAD \
        -master ../anat_mni/anatQQ.${m2}.nii -dxyz 2.5\
        -overwrite -prefix ${m2}_epi_NL_inMNI.nii;
done
```

Then, compute the sphere averages and correlations. But instead of having all
those intermediate files with the ROI timecourse, why not create one big mask
with all the spheres, and use 3dNetCorr like before? Here' I'll use the same
numbering from the xlsx file in the Power atlas
(https://www.jonathanpower.net/2011-neuron-bigbrain.html). 

```bash
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/fmri_same_space/epi;
3dUndump -prefix spheres.nii -master 0411_epi_NL_inMNI.nii -srad 5 \
    -orient LPI -xyz coords.1D
net_dir=/mnt/shaw/MR_data_by_maskid/
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/
for m in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m2=`printf %04d $m`;
    3dNetCorr                                       \
        -inset ${m2}_epi_NL_inMNI.nii                    \
        -in_rois  spheres.nii                       \
        -prefix  ../../fmri_corr_tables/${m2}_power                           \
        -fish_z;
done
```

# TODO

* check alignment!!!
* redo EPI warping for mask ids that hadn't finished copying for the cluster
  (just check which ones don't have warped EPI)
* run script above to compute 3dNetCorr matrices