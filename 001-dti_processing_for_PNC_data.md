# 2019-01-15 12:12:36

I'm following the directions from
https://fsl.fmrib.ox.ac.uk/fslcourse/lectures/practicals/fdt1/index.html to run
the DTI pipeline on the PNC data we just downloaded.

It doesn't seem like FSL can deal with DICOMs, so the first step is to convert
them to NIFTI. I used Paul's tool to do it, and based on his e-mail I just need:

```bash
fat_proc_convert_dcm_dwis  \
        -indir  "DTI_35dir/* DTI_36dir/*"                 \
        -prefix  OFILE
```

For each subject. That generates 71 bvals, 71 vecs, and a .nii.gz with 71 bricks
(each 128 x 128 x 70). Now, we can go on with either the FATCAT or the FSL
pipeline, knowing that there is no need for DRBUDDI.

# 2019-01-18 10:44:08

Following Paul's recommendations by e-mail, let's run one subject all the way.
And by that, I mean including TORTOISE, but assuming that no volumes are bad so
I don't have to do the visual QC.

```bash
cd /data/NCR_SBRB/pnc/dti
tar -zxvf ../600009963128_1.tar.gz
module load afni
Dimon -infile_prefix "600009963128/T1_3DAXIAL/Dicoms/*.dcm" -gert_to3d_prefix 600009963128_t1.nii.gz -gert_create_dataset
fat_proc_convert_dcm_dwis -indir  "600009963128/DTI_35dir/* 600009963128/DTI_36dir/*" -prefix 600009963128_dwi
rm -rf 600009963128
@SSwarper -input 600009963128_t1.nii.gz -base TT_N27_SSW.nii.gz -subid 600009963128
fat_proc_imit2w_from_t1w -inset 600009963128_t1_ax.nii.gz -prefix 600009963128_t2_ax_immi -mask anatSS.600009963128.nii
# I got the phase information after using ImportDICOM tool from TORTOISE and checking the .list file
DIFFPREP --dwi 600009963128_dwi.nii --bvecs 600009963128_dwi_rvec.dat --bvals 600009963128_dwi_bval.dat --structural 600009963128_t2_ax_immi.nii --phase vertical
@GradFlipTest -in_dwi 600009963128_dwi_DMC.nii -in_col_matT 600009963128_dwi_DMC.bmtxt -prefix 600009963128_GradFlipTest_rec.txt
my_flip=`cat 600009963128_GradFlipTest_rec.txt`;
fat_proc_dwi_to_dt \
    -in_dwi       600009963128_dwi_DMC.nii                    \
    -in_col_matT  600009963128_dwi_DMC.bmtxt                  \
    -in_struc_res 600009963128_dwi_DMC_structural.nii               \
    -in_ref_orig  600009963128_dwi_DMC_template.nii          \
    -prefix       600009963128_dwi                           \
    -mask_from_struc                                   \
    $my_flip
fat_proc_decmap                                     \
    -in_fa       dt_FA.nii.gz     \
    -in_v1       dt_V1.nii.gz     \
    -mask        600009963128_dwi_mask.nii.gz  \
    -prefix      DEC
```

Note that we'll need to stop in the middle to allow for the IRTAs to do the
visual QC. So, let's create a wrapper script that does some of the steps above:

```bash
cd /data/NCR_SBRB/pnc
for m in `cat have_imaging.txt`; do
    echo "bash ~/research_code/dti/tortoise_pnc_wrapper.sh ${m}" >> swarm.tortoise;
done;
swarm -g 10 -t 16 --job-name tortoise --time 4:00:00 -f swarm.tortoise \
    -m afni,TORTOISE --partition quick --logdir trash