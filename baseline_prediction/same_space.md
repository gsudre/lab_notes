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
    # mylink=`readlink /Volumes/Shaw/data_by_maskID/${m2} | sed "s/\.\./\/Volumes\/Shaw/"`;
    # 3dcopy ${mylink}/afni/${m2}.rest.subjectSpace.results/errts.${m2}.fanaticor+orig.HEAD ${m2}_epi.nii;
    @auto_tlrc -apar ../anat/${m2}_inTLRC.nii -input ${m2}_epi.nii -suffix _inTLRC -dxyz 2.5;
done;
```

Should probably check that the alignment worked well though... I'll go on with
converting rsFMRI as well, and **NEED TO CHECK ALIGNMENT LATER!**

# dti

Let's go ahead and do this for DTI as well, which should be just a matter of
making symlinks as we have all these data already generates for QC purposes:

```bash
for n in 223 272; do
    cd ~/data/baseline_prediction/same_space/dti/n${n}
    # ignore header
    for m in `tail -n +2 ~/data/baseline_prediction/dti_gf_09212018_${n}timeDiff12mo.csv | cut -d"," -f 2 -`; do
        maskid=`printf %04d $m`;
        for p in fa ad rd; do
            ln -s /Volumes/Shaw/dti_robust_tsa/analysis_may2017/${maskid}_tensor_diffeo_${p}.nii.gz .;
        done;
    done;
done
```

But then I gunziped and split everything by directory, because the Matlab tool
was having lots of issues listing the directory.

# 2018-11-02 15:38:15

I've been generally running one version with 30 ICs, an then another one using
the estimated number of ICs. Both, using MST500 for stability. I need to study a
bit on what that actually does, because I don't see a stability metric like
ICASSO. Of course, I could just do ICASSO if I want to. `

Another thing to look into is a better mask, because for the strucutral analysis
I'm getting a lot of white matter and ventricle action. So, having a nicer gray
matter mask would be better than just using the automask. That will be
especially important for DTI too, if the property files didn't use the mask when
being generated.

Of course, STILL NEED TO CHECK THAT TLRC ALIGNMENT WORKED WELL!
