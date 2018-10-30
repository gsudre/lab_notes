# 2018-10-29 15:56:46

I'm going to do a quick push to try to ru all images in the same space using
ICA. Let's see what that gives us. First step is to convert all our images to
TT_N27 space, and we can use that matrix later to put the EPI in the same space
(https://afni.nimh.nih.gov/pub/dist/edu/latest/afni_handouts/afni10_volreg_talairach.pdf)


# 2018-10-30 08:00:53

```bash
cd ~/data/baseline_prediction/same_space/anat;
for m in `cut -d"," -f 1 ../../struct_rois_09062018_260timeDiff12mo.csv`; do
    m2=`printf %04d $m`;
    mri_convert /Volumes/Shaw/freesurfer5.3_subjects/${m2}/mri/orig.mgz ./${m2}.nii.gz;
    @auto_tlrc -base TT_N27+tlrc -input ${m2}.nii.gz -no_pre -suffix _inTLRC;
done;
```

then, just so we can run everything in my desktop without depending on caterpie for data_by_maskID:

```bash
cd ~/data/baseline_prediction/same_space/epi;
for m2 in `cat ../../rsfmri_3minWithClinical.tsv`; do
    echo $m2;
    mylink=`readlink /Volumes/Shaw/data_by_maskID/${m2} | sed "s/\.\./\/Volumes\/Shaw/"`;
    3dcopy ${mylink}/afni/${m2}.rest.subjectSpace.results/errts.${m2}.fanaticor+orig.HEAD ${m2}_epi.nii;
    #@auto_tlrc -apar ${m2}.nii -input ${m2}_epi.nii -suffix _inTLRC -dxyz 2.5;
done;
```

Should probably check that the alignment worked well though...

