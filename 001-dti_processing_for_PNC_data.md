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
```

# 2019-01-29 10:38:23

I was chatting with Ryan, and apparently their approach for GenR is to process
all volumes all the time. Then, use quantitative variables to figure out which
subjects to remove. We can try that approach here, especially considering that
we will likely run DTIPrep to get some additional QC metrics, and we rarely ever
use the subjects when we remove more than a couple volumes anyways.

So, let's run the first 100 subjects in the new TORTOISE pipeline, then we can do
the same for the FSL pipeline once I have that working.

```bash
cd /data/NCR_SBRB/pnc
for m in `cat first100.txt`; do
    echo "bash ~/research_code/dti/tortoise_pnc_wrapper.sh ${m}" >> swarm.tortoise;
done;
swarm -g 10 -t 16 --job-name tortoise --time 4:00:00 -f swarm.tortoise \
    -m afni,TORTOISE --partition quick --logdir trash
```

For the FSL pipeline, we can do something like this:

```bash
fslroi dwi b0 0 1
bet b0 b0_brain -m -f 0.2
idx=''; for i in {1..71}; do a=$a' '1; done; echo $a > index.txt
echo "0 -1 0 0.102" > acqparams.txt
eddy_openmp --imain=dwi --mask=b0_brain_mask --index=index.txt --acqp=acqparams.txt --bvecs=dwi_cvec.dat --bvals=dwi_bval.dat --fwhm=0 --flm=quadratic --out=eddy_unwarped_images --cnr_maps --repol --mporder=6

# OR

sinteractive --gres=gpu:k20x:1
module load CUDA/7.5
eddy_cuda --imain=dwi --acqp=acqparams.txt --index=index.txt --mask=b0_brain_mask --bvals=dwi_bval.dat --bvecs=dwi_cvec.dat --out=eddy_s2v_unwarped_images --niter=8 --fwhm=10,6,4,2,0,0,0,0 --repol --ol_type=both --mporder=8 --s2v_niter=8 --slspec=my_slspec.txt --cnr_maps

dtifit --data=eddy_s2v_unwarped_images --mask=b0_brain_mask --bvals=dwi_bval.dat --bvecs=dwi_cvec.dat --sse --out=dti

```

I got the parameters for acquisition from running dcm2niix_afni on both DTI
sequences (35 and 35), and then looking for Phase, PE, and Echo in the jsons. In
the end, they matched the example in
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_do_I_know_what_to_put_into_my_--acqp_file.

I also used the json to construct the myslspec.txt file, using the Matlab code
from
https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_should_my_--slspec_file_look.3F

So, I put all of that in a script, which I ran like this:

```bash
cd /data/NCR_SBRB/pnc
for m in `cat first100.txt`; do
    echo "bash ~/research_code/dti/fdt_pnc_wrapper.sh ${m}" >> swarm.fdt;
done;
swarm -g 4 --job-name fdt --time 4:00:00 -f swarm.fdt --partition gpu \
    --logdir trash_fdt --gres=gpu:k80:2
```

# 2019-02-12 17:38:14

I made a few changes to the fdt wrapper to include Qc pictures and warping to
the FMRIB58 space. The results look decent, so it's time to run it for everybody
and start looking at the QC pictures. I'll then write wrappers to grab the
motion and outlier parameters, as well as something to collect the tract
averages, brain QC images (to be more organized for easy QCing), do the rest of the TBSS analysis, and also output the SSE and CNR maps
for further QCing.

Let's go ahead and all 647 overnight.

# 2019-02-13 11:31:26

I'll ignore the CNR maps for now. We have plenty here after chaking the brain
masks, warping, and SSE. We're also looking for numbers of outliers and movement
variables, so let's grab those.

```bash
out_fname=mvmt_report.csv;
echo "id,Noutliers,norm.trans,norm.rot,RMS1stVol,RMSprevVol" > $out_fname;
for m in `cat myids.txt`; do
    echo 'Collecting metrics for' $m;
    if [ -e ${m}/eddy_s2v_unwarped_images.eddy_outlier_report ]; then
        noutliers=`cat ${m}/eddy_s2v_unwarped_images.eddy_outlier_report | wc -l`;
        1d_tool.py -infile ${m}/eddy_s2v_unwarped_images.eddy_movement_over_time \
            -select_cols '0..2' -collapse_cols euclidean_norm -overwrite \
            -write trans_norm.1D;
        trans=`1d_tool.py -infile trans_norm.1D -show_mmms | \
            tail -n -1 | awk '{ print $8 }' | sed 's/,//'`;
        1d_tool.py -infile ${m}/eddy_s2v_unwarped_images.eddy_movement_over_time \
            -select_cols '3..5' -collapse_cols euclidean_norm -overwrite \
            -write rot_norm.1D;
        rot=`1d_tool.py -infile rot_norm.1D -show_mmms | \
            tail -n -1 | awk '{ print $8 }' | sed 's/,//'`;
        1d_tool.py -infile ${m}/eddy_s2v_unwarped_images.eddy_movement_rms \
            -show_mmms > mean_rms.txt;
        vol1=`head -n +2 mean_rms.txt | awk '{ print $8 }' | sed 's/,//'`;
        pvol=`tail -n -1 mean_rms.txt | awk '{ print $8 }' | sed 's/,//'`;
    else
        echo "Could not find outlier report for $m"
        noutliers='NA';
        trans='NA';
        rot='NA';
        vol1='NA';
        pvol='NA';
    fi;
    echo $m, $noutliers, $trans, $rot, $vol1, $pvol >> $out_fname;
done
```

Note that Ryan's output usually comes from bedpostX, so I'll still need to run
that. Or I can just go for a regular average over the mask mean. For example,
see Ryan's e-mail from December 06, 2016 10:11 AM.

I'm not going to go the CAMINO way, but we might end up running autoPtx. Need to
see how it looks at Biowulf though:

https://hpc.nih.gov/apps/fsl.html

Ryan's pipeline seems to spit out both OLS and RESTORE estimates, and he says
they're highly correlated. Just so we don't have to jump between programs, let's
stic to the OLS etimates for now.

For bedpostx, I did something like this:

```bash
ln -s eddy_old/eddy_s2v_unwarped_images.nii.gz data.nii.gz
ln -s dwi_bval.dat bvals
ln -s eddy_old/eddy_s2v_unwarped_images.eddy_rotated_bvecs bvecs
ln -s b0_brain_mask.nii.gz nodif_brain_mask.nii.gz
bedpostx ./
```

and that schedules a whole bunch of swarms just for the single subject. Over 60,
I'd say, wach of 2 cores, 2h. Shouldn't take too long to run in parallel, but
that's per subject, so it'll take a while. It can run faster if I use GPU, but
then I'm limited on how many GPUs I can allocate. Let's think more about it
later.

Also note that from Joelle's e-mails, we should only do one fiber orientation for our data,
instead of 2 in bedpost.

Yay! These seem to work when converting our data:

```bash
tail -n -60 cdiflist08 | split -l 20
module load TORTOISE
ImportDICOM -i mr_0008 -o s1 -b 1100 -g xaa
ImportDICOM -i mr_0009 -o s2 -b 1100 -g xab
ImportDICOM -i mr_0010 -o s3 -b 1100 -g xac
TORTOISEBmatrixToFSLBVecs s1_proc/s1.bmtxt 
TORTOISEBmatrixToFSLBVecs s2_proc/s2.bmtxt 
TORTOISEBmatrixToFSLBVecs s3_proc/s3.bmtxt 
fat_proc_convert_dcm_dwis -no_qc_view -innii s1_proc/s1.nii s2_proc/s2.nii s3_proc/s3.nii -inbval  s1_proc/s1.bvals s2_proc/s2.bvals s3_proc/s3.bvals -inbvec  s1_proc/s1.bvecs s2_proc/s2.bvecs s3_proc/s3.bvecs -prefix dwi_comb -flip_x
```

That matched exactly what I got from opening the .list file in DIFFCALC and
exporting the raw DWI (edti.list) to FSL UNSORTED format.

After this, it's just the exact same thing as the PNC pipeline, with probably
some changed in the bedpostx to only fit one orientation. In fact, I could do
that for the PNC as well, just so we don't have two different options.
