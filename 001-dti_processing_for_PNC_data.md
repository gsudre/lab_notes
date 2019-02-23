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
echo "id,Noutliers,PCToutliers,NoutVolumes,norm.trans,norm.rot,RMS1stVol,RMSprevVol" > $out_fname;
for m in `cat myids.txt`; do
    echo 'Collecting metrics for' $m;
    if [ -e ${m}/eddy_s2v_unwarped_images.eddy_outlier_report ]; then
        noutliers=`cat ${m}/eddy_s2v_unwarped_images.eddy_outlier_report | wc -l`;
        # figuring out the percetnage of total slices the outliers represent
        nslices=`tail ${m}/eddy_s2v_unwarped_images.eddy_outlier_map | awk '{ print NF; exit } '`;
        nvol=`cat ${m}/dwi_cvec.dat | wc -l`;
        let totalSlices=$nslices*$nvol;
        let pctOutliers=$noutliers/$totalSlices;
        # figuring out how many volumes were completely removed (row of 1s)
        awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' \
            ${m}/eddy_s2v_unwarped_images.eddy_outlier_map > outlier_avg.txt;
        nOutVols=`grep -c -e "^1$" outlier_avg.txt`;
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
        pctOutliers='NA';
        nOutVols='NA';
        trans='NA';
        rot='NA';
        vol1='NA';
        pvol='NA';
    fi;
    echo $m, $noutliers, $pctOutliers, $nOutVols, $trans, $rot, $vol1, $pvol >> $out_fname;
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
ln -s eddy_s2v_unwarped_images.nii.gz data.nii.gz
ln -s dwi_bval.dat bvals
ln -s eddy_s2v_unwarped_images.eddy_rotated_bvecs bvecs
ln -s b0_brain_mask.nii.gz nodif_brain_mask.nii.gz
bedpostx ./
```

and that schedules a whole bunch of swarms just for the single subject. Over 60,
I'd say, wach of 2 cores, 2h. Shouldn't take too long to run in parallel, but
that's per subject, so it'll take a while. It can run faster if I use GPU, but
then I'm limited on how many GPUs I can allocate. Let's think more about it
later. But the command would be bedpostx_gpu.

Also note that from Joelle's e-mails, we should only do one fiber orientation for our data,
instead of 2 that is default in bedpost.

Before I go nuts running bedpost on everyone, let's collect the QC that's
already finished:

```bash
mkdir /data/NCR_SBRB/pnc/dti_fdt/summary_QC
cd /data/NCR_SBRB/pnc/dti_fdt/summary_QC/
mkdir brainmask
mkdir transform
mkdir DEC
mkdir SSE
for m in `cat ../myids.txt`; do
    cp ../${m}/QC/brain_mask.axi.png brainmask/${m}.axi.png
    cp ../${m}/QC/brain_mask.sag.png brainmask/${m}.sag.png
    cp ../${m}/QC/brain_mask.cor.png brainmask/${m}.cor.png

    cp ../${m}/QC/FA_transform.axi.png transform/${m}.axi.png
    cp ../${m}/QC/FA_transform.sag.png transform/${m}.sag.png
    cp ../${m}/QC/FA_transform.cor.png transform/${m}.cor.png

    cp ../${m}/QC/DEC_qc_dec_sca07.axi.png DEC/${m}.axi.png
    cp ../${m}/QC/DEC_qc_dec_sca07.sag.png DEC/${m}.sag.png
    cp ../${m}/QC/DEC_qc_dec_sca07.cor.png DEC/${m}.cor.png

    cp ../${m}/QC/sse.axi.png SSE/${m}.axi.png
    cp ../${m}/QC/sse.cor.png SSE/${m}.cor.png
    cp ../${m}/QC/sse.sag.png SSE/${m}.sag.png
done
```

To fill in the QC spreadsheet, we check who converted properly:

```bash
cd /data/NCR_SBRB/pnc/dti_fdt
for m in `cat ../have_imaging.txt`; do
    nvol=`cat ${m}/dwi_cvec.dat | wc -l`;
    if [ ! $nvol = 71 ]; then
        echo $m,$nvol >> ~/tmp/conversion_errors.txt;
    fi;
done
```

Then, check that all brain masks were created:

```bash
cd /data/NCR_SBRB/pnc/dti_fdt
for m in `cat converted.txt`; do
    if [[ -e ${m}/QC/brain_mask.axi.png && -e ${m}/b0_brain_mask.nii.gz ]]; then
        echo $m,y >> ~/tmp/mask_status.txt;
    else
        echo $m,n >> ~/tmp/mask_status.txt;
    fi;
done
```

And then who has eddy:

```bash
cd /data/NCR_SBRB/pnc/dti_fdt
for m in `cat converted.txt`; do
    if [[ -e ${m}/eddy_s2v_unwarped_images.eddy_rotated_bvecs && -e ${m}/eddy_s2v_unwarped_images.nii.gz ]]; then
        echo $m,y >> ~/tmp/eddy_status.txt;
    else
        echo $m,n >> ~/tmp/eddy_status.txt;
    fi;
done
```

To run bedpostx, I'll split the data into 6, so that each account can run a GPU
to its limit and then a CPU one. Things might go faster that way.

```bash
for m in `cat xaf`; do
    cd /data/NCR_SBRB/pnc/dti_fdt/${m};
    ln -s eddy_s2v_unwarped_images.nii.gz data.nii.gz;
    ln -s dwi_bval.dat bvals;
    ln -s eddy_s2v_unwarped_images.eddy_rotated_bvecs bvecs;
    ln -s b0_brain_mask.nii.gz nodif_brain_mask.nii.gz;
    bedpostx_gpu -n 1 ./;
done
```

Run autoPtx first to make sure everything is fine!!!!

I had to make some changes to how I'm approaching this. That's just because
running autoPtx makes life much easier, but it redoes some of the steps I was
doing in the wrapper before. So, let's stop the wrapper right before we
calculate dtifit (which is done by autoPtx), and we can generate the QC pictures
afterwards. 

autoPtx expects the data in the same format as what we'd send for bedpostx, so
let's do that formatting first. But keep in mind that symlinks don't work, and
the data is in the end moved to preproc. So, let's actually copy them, because I
don't want to risk losing the eddy output:

```bash
# run in helix so we don't overload BW filesystem
for m in `cat xac`; do
    echo Copying $m;
    cd /data/NCR_SBRB/pnc/dti_fdt/${m};
    cp eddy_s2v_unwarped_images.nii.gz data.nii.gz;
    cp dwi_bval.dat bvals;
    cp eddy_s2v_unwarped_images.eddy_rotated_bvecs bvecs;
    cp b0_brain_mask.nii.gz nodif_brain_mask.nii.gz;
done
```

Then, we run autoPtx. We can split it by users because it adds everything to the
same directory, and just increments the final subject list.

Run it a long interactive session, because even though it schedules bedpostx, it
still runs all kinds of registrations through FSL, so biowulf headnode won't cut
it!

```bash
data='';
for m in `cat xac`; do
    data=$data' '${m}/data.nii.gz;
done
/data/NCR_SBRB/software/autoPtx/autoPtx_1_preproc $data;
```

And of course we still need part 2 when we're done.

# 2019-02-20 10:58:43

I changed the part 2 script so run the scans split between two accounts.
Hopefully there won't be issues with permissions...

But I had to reduce the number of subjects per file because we CPU recruitment
limits in the cluster. Let's try only 100. And still, better to only fire new
ones when nothing else is queued (running is OK).


It creates one swarm per tract, with one job per subject in each tract. So,
nsubjects * ntracts. I might need to use the subject file as an argument to
split it across accounts.

But I don't need to wait for the second part to do:

```bash
for m in `cat ../xab`; do
    bash ~/research_code/dti/fdt_pnc_TBSS_and_QC.sh ${m};
done
```

Now we need to copy the QC images again:

```bash
qc_dir=/data/NCR_SBRB/pnc/dti_fdt/summary_QC/
img_dir=/data/NCR_SBRB/pnc/dti_fdt/preproc/
for m in `cat ~/tmp/pnc_qc.txt`; do
    cp $img_dir/${m}/QC/FA_transform.axi.png $qc_dir/transform/${m}.axi.png
    cp $img_dir/${m}/QC/FA_transform.sag.png $qc_dir/transform/${m}.sag.png
    cp $img_dir/${m}/QC/FA_transform.cor.png $qc_dir/transform/${m}.cor.png

    cp $img_dir/${m}/QC/DEC_qc_dec_sca07.axi.png $qc_dir/DEC/${m}.axi.png
    cp $img_dir/${m}/QC/DEC_qc_dec_sca07.sag.png $qc_dir/DEC/${m}.sag.png
    cp $img_dir/${m}/QC/DEC_qc_dec_sca07.cor.png $qc_dir/DEC/${m}.cor.png

    cp $img_dir/${m}/QC/sse.axi.png $qc_dir/SSE/${m}.axi.png
    cp $img_dir/${m}/QC/sse.cor.png $qc_dir/SSE/${m}.cor.png
    cp $img_dir/${m}/QC/sse.sag.png $qc_dir/SSE/${m}.sag.png
done
```

# 2019-02-22 10:07:06

While the IRTAs QC the resulting images, I'll go ahead and start copying data to
shaw/PNC_DTI. I'm using Globus web interface, because I'll go for all
directories for now. This is justa  way to keep the results in our servers in
case I need to make room in BW.