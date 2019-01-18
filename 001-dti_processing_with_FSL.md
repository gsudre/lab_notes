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
```

Note that we'll need to stop in the middle to allow for the IRTAs to do the
visual QC. In any case, to swarm it we will do:

<!-- 
```bash
cd ~/data/baseline_prediction/same_space/anat;
for m in `cut -d"," -f 1 ../../struct_rois_09062018_260timeDiff12mo.csv`; do
    m2=`printf %04d $m`;
    echo "export OMP_NUM_THREADS=16; cd ~/data/baseline_prediction/same_space/anat; @SSwarper -input ${m2}.nii.gz -base TT_N27_SSW.nii.gz -subid ${m2}" >> swarm.SSwarper
done;
swarm -g 10 -t 16 --job-name SSwarper --time 2:00:00 -f swarm.SSwarper -m afni --partition quick --logdir trash
``` -->