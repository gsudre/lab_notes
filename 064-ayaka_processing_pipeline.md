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

If you have a fancier computer, or at least anything that has GPU in it, it's
worth running outlier correction in eddy. But because we put together two
different sequences here, I'm not going to use slice-to-volume correction.

```
eddy_openmp --imain=topup_corrected --acqp=acqparams.txt --index=index.txt \
    --mask=topup_b0_brain_mask --bvals=PA.bval --bvecs=PA.bvec \
    --out=eddy_ol_unwarped_images --niter=8 --fwhm=10,6,4,2,0,0,0,0 \
    --repol --ol_type=sw --flm=quadratic  --cnr_maps --data_is_shelled \
    --ol_nstd=4
```

