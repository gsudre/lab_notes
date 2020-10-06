# 2020-10-06 11:02:09

Let's see what we can find on Dennis' data.

```bash
# bw
cd ~/data/dennis/24
module load afni
dcm2niix_afni *
```

Looking at the JSONs, I saw that 24__20191126213428_13 was the DTI. Let's go
with the FDT pipeline for processing because TORTOISE takes forever.

I copied that the scan files to a raw folder to keep it tidy.

```bash
# bw
module load fsl

fslreorient2std 24__20191126213428_13.nii dwi
fslroi dwi b0 0 1
bet b0 b0_brain -m -f 0.2
cp 24__20191126213428_13.bval bvals
cp 24__20191126213428_13.bvec bvec
# making sure the directions are OK
dtifit --data=dwi --mask=b0_brain --bvecs=bvecs --bvals=bvals --out=tmp_dti
```

Looks fine I think. Let's go on with processing.

```bash
a=''; for i in {1..102}; do 
    a=$a' '1;
done;
echo $a > index.txt
echo "0 1 0 0.05" > acqparams.txt
echo "0 -1 0 0.05" >> acqparams.txt

eddy_openmp --imain=dwi --mask=b0_brain --index=index.txt \
     --acqp=acqparams.txt \
     --bvecs=bvecs --bvals=bvals --fwhm=0 --flm=quadratic \
     --out=eddy_unwarped_images
```

And I'll follow the same steps for the other subjects (at least until I figure
out whether I'm importing them properly)?

