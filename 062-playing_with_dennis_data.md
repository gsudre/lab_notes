# 2019-12-12 13:53:58

Started by copying it to ~/data. Then, I used dcm2nii just to get the converted
files:

```bash
cd ~/data/NICHD_PROJECT/
cd MRI_DATA/DICOM_data/
cd 34470001/
dcm2niix_afni -z y *
cd ../34470000
dcm2niix_afni -z y *
```

There's one candidate file, with 50 directions. LEt's play with that one in FSL
then:

```bash
fname=34470000_ep2d_diff_mddw_20_p2_20190628193654_13
fslreorient2std ${fname}.nii.gz dwi
fslroi dwi b0 0 1
bet b0 b0_brain -m -f 0.2
idx=''; for i in {1..50}; do 
    a=$a' '1;
done;
echo $a > index.txt
echo "0 -1 0 0.102" > acqparams.txt

eddy --imain=dwi --mask=b0_brain_mask --index=index.txt \
    --acqp=acqparams.txt --bvecs=${fname}.bvec --bvals=${fname}.bval \
    --fwhm=0 --flm=quadratic --out=eddy_unwarped_images

# ... takes a while ...

# copying over some files to their correct names for bedpostX
cp eddy_s2v_unwarped_images.nii.gz data.nii.gz;
cp ${fname}.bval bvals;
cp ${fname}.bvec old_bvecs;
cp eddy_unwarped_images.eddy_rotated_bvecs bvecs;
cp b0_brain_mask.nii.gz nodif_brain_mask.nii.gz;
```

