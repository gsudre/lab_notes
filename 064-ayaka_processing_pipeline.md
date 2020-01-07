# 2020-01-07 10:43:11

## DTI
Starting with DTI, we'll use FSL FDT:

```bash
cd /scratch/sudregp/
mkdir A01
cd "dMRI_dir107_PA - 21"/
dcm2niix_afni -z y -f PA *
mv PA* ../A01/
cd "../dMRI_dir107_AP - 19"/
dcm2niix_afni -z y -f AP *
mv AP* ../A01/
```

Let's re-orient each file and rename them to make sure we're in the correct
space to begin with:

```bash
cd /scratch/sudregp/A01/
fname=AP
fslroi ${fname} b0 0 1
bet b0 b0_brain -m -f 0.2
```

Now it's mostly following directions from here:

https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide

```bash
dtifit --data=${fname}.nii.gz --mask=b0_brain_mask --bvals=${fname}.bval \
    --bvecs=${fname}.bvec --sse --out=dti
```

Checking that the data (V1) looks good... seems fine to me. I'll assume for now
the other vector field is also ccorrect, but we should always check all of this
in the end, after preprocessing, anyways.

```bash
# running TOPUP
fslroi AP nodif 0 1
fslroi PA nodif_PA 0 1
fslmerge -t AP_PA_b0 nodif nodif_PA
```

By looking at the JSON and then here:

https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/Faq

I think our parameters are .68 and 140. So, we do:

```bash
echo "0 -1 0 0.095" > acqparams.txt
echo "0 1 0 0.095" >> acqparams.txt
# this takes a while
topup --imain=AP_PA_b0 --datain=acqparams.txt --config=b02b0.cnf \
    --out=topup_AP_PA_b0 --iout=topup_AP_PA_b0_iout --fout=topup_AP_PA_b0_fout
```

To check the TOPUP results, load topup_AP_PA_b0_iout and compare its two volumes
to those in AP_PA_b0.nii.gz.

If that looks good, we will only need to send it to eddy, instead of applying it
beforehand. In other words, if using eddy, no need to use applytopup!

Here's a curious step: from then on FDT only makes use of the first
acquisition from according to the documentation. In other words, it only uses AP and PA for TOPUP, but
the rest of the analysis is only conducted in one of the datasets? That's
odd. I can see the rationale, as structural data shouldn't change, but what
about data quality and movement? Should one be making use of both datasets to
begin with?

That doesn't make sense to me. So, let's apply TOPUP first, and use that
instead:

```bash
applytopup --imain=AP,PA --topup=topup_AP_PA_b0 --datain=acqparams.txt \
    --inindex=1,2 --out=topup_corrected
```

Next step is applying eddy. Note that bval and bvec are the same in AP and PA!

```bash
fslmaths topup_AP_PA_b0_iout -Tmean topup_b0_mean
bet topup_b0_mean topup_b0_brain -m -f 0.2
idx=''; for i in {1..108}; do 
    a=$a' '1;
done;
echo $a > index.txt
eddy_openmp --imain=topup_corrected --mask=topup_b0_brain_mask \
    --index=index.txt --acqp=acqparams.txt --bvecs=PA.bvec --bvals=PA.bval \
    --fwhm=0 --flm=quadratic --out=eddy_unwarped_images --data_is_shelled
```

If you have a fancier computer, it's
worth running outlier correction in eddy. But because we put together two
different sequences here, I'm not going to use slice-to-volume correction.

```
eddy_openmp --imain=topup_corrected --acqp=acqparams.txt --index=index.txt \
    --mask=topup_b0_brain_mask --bvals=PA.bval --bvecs=PA.bvec \
    --out=eddy_ol_unwarped_images --niter=8 --fwhm=10,6,4,2,0,0,0,0 \
    --repol --ol_type=sw --flm=quadratic  --cnr_maps --data_is_shelled \
    --ol_nstd=4
```

Regardless of what is used, it's good to run some eddy QC:

https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddyqc/UsersGuide


<!-- Then, it's just a matter of checking the quality of the results and running
autoPtx and TBSS:

```bash
# copying over some files to their correct names for bedpostX
cp eddy_s2v_unwarped_images.nii.gz data.nii.gz;
cp PA.bval bvals;
cp PA.bvec old_bvecs;
cp eddy_unwarped_images.eddy_rotated_bvecs bvecs;
cp topup_b0_brain_mask.nii.gz nodif_brain_mask.nii.gz;
bedpostx_gpu -n 1 ./;


autoPTx

tbss_1_preproc dti_FA.nii.gz

cp origdata/dti_FA.nii.gz ./
cp FA/dti_FA_FA.nii.gz dti_FA_eroded.nii.gz

# make directionality encoded QC pictures, checking that eroded FA looks fine.
# Note that the dtifit results, and the std alignment were run by autoPtx!
fat_proc_decmap -in_fa dti_FA_eroded.nii.gz -in_v1 dti_V1.nii.gz \
    -mask nodif_brain_mask.nii.gz -prefix DEC

# apply the transform calculated by autoPtx to a few maps. code copied from 
# tbss_non_fa
for f in FA L1 L2 L3 MD MO FA_eroded; do
    echo Warping $f;
    applywarp -i dti_${f} -o ${f}_in_FMRIB58_FA_1mm \
        -r $FSLDIR/data/standard/FMRIB58_FA_1mm -w nat2std_warp
done

# # make transformation QC figure: warped subject B0 is the overlay!
# @chauffeur_afni                             \
#     -ulay  $FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz                     \
#     -olay  FA_in_FMRIB58_FA_1mm.nii.gz                         \
#     -ulay_range 0% 150%                     \
#     -func_range_perc 50                     \
#     -pbar_posonly                           \
#     -cbar "red_monochrome"                  \
#     -opacity 8                              \
#     -prefix   QC/FA_transform              \
#     -montx 3 -monty 3                       \
#     -set_xhairs OFF                         \
#     -label_mode 1 -label_size 3             \
#     -do_clean

# new version with FSL template as the edges
@snapshot_volreg FA_in_FMRIB58_FA_1mm.nii.gz \
    $FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz \
    QC/FA_transform;

# make QC images for standard errors. Here we set our color scale to have 95th
# percentile of all errors. Meaning, more red = bigger error.
@chauffeur_afni                             \
    -ulay  data.nii.gz                       \
    -olay  dti_sse.nii.gz                          \
    -opacity 5                              \
    -pbar_posonly   \
    -cbar Spectrum:red_to_blue              \
    -set_subbricks 0 0 0     \
    -prefix   QC/sse              \
    -montx 6 -monty 6                       \
    -set_xhairs OFF                         \
    -label_mode 1 -label_size 3             \
    -thr_olay 0 \
    -func_range_perc_nz 95 \
    -do_clean

# we can derive the skeleton later, either based on FMRIB58 or group

# note that we need to look at
# https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/TBSS/UserGuide#Using_non-FA_Images_in_TBSS
# to run the non-FA files! -->




