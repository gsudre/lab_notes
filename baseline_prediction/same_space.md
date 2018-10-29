# 2018-10-29 15:56:46

I'm going to do a quick push to try to ru all images in the same space using
ICA. Let's see what that gives us. First step is to convert all our images to
TT_N27 space, and we can use that matrix later to put the EPI in the same space
(https://afni.nimh.nih.gov/pub/dist/edu/latest/afni_handouts/afni10_volreg_talairach.pdf)

```bash
mri_convert /Volumes/Shaw/freesurfer5.3_subjects/0119/mri/orig.mgz ./0119.nii.gz
@auto_tlrc -base TT_N27+tlrc -input 0119.nii.gz -suffix _inTLRC
3dcopy /Volumes/Shaw/MR_data/SUVJECT_NAME/2217/afni/2217.rest.subjectSpace.results/errts.2217.fanaticor+orig test_epi.nii
@auto_tlrc -apar 0119_inTLRC.nii -input test_epi.nii -suffix _inTLRC -dxyz 2
```