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

# 2018-11-05 11:16:22

I figured that since I'm doing a lot of this work, we might as well usea a
non-Linear transform, and also produce all the figures automatically. That's
where the new SSwarp function comes in.

On the same step, I think we'll have better success in ICA if we
use gray matter masks (or white matter, for DTI, if the files we're using are
not restricted to that yet). So, for the structural data we'll need to use
Freesurfer masks and put them into the common space.

Ideally I'd run this in Biowulf, but until I hear from them about Xfvb and
netpbm, I'l run it locally:

```bash
cd ~/data/baseline_prediction/same_space/anat;
for m in `cut -d"," -f 1 ../../struct_rois_09062018_260timeDiff12mo.csv`; do
    @SSwarper -input ${m2}.nii.gz -base TT_N27_SSW.nii.gz -subid ${m2};
    3dcalc -a /Volumes/Shaw/freesurfer5.3_subjects/${m2}/SUMA/aparc+aseg_REN_gm.nii.gz \
        -prefix ${m2}_gm_mask.nii -expr "step(a)";
    3dNwarpApply -nwarp "anatQQ.${m2}_WARP.nii anatQQ.${m2}.aff12.1D" \
        -source ${m2}_gm_mask.nii -master anatQQ.${m2}.nii -inter NN \
        -overwrite -prefix ${m2}_gm_mask_NL_inTLRC.nii;
done;
```

For DTI, they're all in the same space already, so I'll just use the FA skeleton
mask we've been using for the voxelwise analysis.

Just needed to say that using the fa_skeleton masks causes the estimates of ICs
to be quite low. For both n223 and n272 AD, I get only 4 ICs. For n223 FA, I got
only 1, so I'm not running faEstimatedMaskedMST200. I'll check n272 and also rd,
but I might end up skipping those estimated results in the end.

Yep, only one IC for FA n272 as well. Will only run IC20, and then check on RD. 