# 2021-05-18 07:55:10

We'll need to put it all in BIDS to run fmriprep and possibly tractorFlow and
mriqc. I'm not sure if TractorFlow will work (testing it now), and if fmriprep
already runs mriqc.

At this moment, dcm2bids seems to have worked in ncrshell01, but I'm checking if
TractorFlow will find the B0 map and use it in topup. If that works, we can use
fmriprep in the BIDS folder. I only ran one subject through dcm2bids in
ncrshell01 just to check if it'd work. 

```bash
# ncrshell01
conda activate dcm2bids
cd /mnt/shaw/sudregp/PNC_BIDS
dcm2bids -d ../tmp/600009963128/ -p 600009963128 -c ~/research_code/pnc2bids.json
```

Then, I'm currently trying TractorFlow. First, checking if topup will find the
B0 phase distortion file. If it does, we can go ahead and run dcm2bids for
everyone in PNC so we can start running fmriprep. But we'll still need to figure
out the best way to run TractorFlow, in terms of parameters, etc

```bash
# bw
module load singularity/3.7.3
module load nextflow/21.04.1

nextflow run /data/NCR_SBRB/software/tractoflow/main.nf --bids /scratch/sudregp/PNC_BIDS/ --dti_shells "0 1000" --fodf_shells "0 1000 2000" -with-singularity tractoflow_2.2.1_b9a527_2021-04-13.sif -resume -profile fully_reproducible --processes 32
```

Will need to set it up as sbatch as well, like in the TractorFlow documentation.
Likely passing subject id to it.

This was getting to cumbersome to get TractoFlow to recognize the blipup/down
input, so I'll just reformat the BIDS output to conform to TractoFlow's original
root directory:

```bash
# ncrshell
fslmerge -t rev_b0 rev_b0_e1_pha rev_b0_e2_pha
fslmerge -t dwi sub-600009963128_run-01_dwi sub-600009963128_run-02_dwi
paste sub-600009963128_run-01_dwi.bvec sub-600009963128_run-02_dwi.bvec > bvec
paste sub-600009963128_run-01_dwi.bval sub-600009963128_run-02_dwi.bval > bval
cp ../../PNC_BIDS/sub-600009963128/anat/sub-600009963128_T1w.nii.gz t1.nii.gz

nextflow run tractoflow/main.nf --input /mnt/shaw/sudregp/tmp/TF/ --participants_label "s1" --dti_shells "0 1000" --fodf_shells "0 1000 2000" -with-singularity tractoflow_2.2.1_b9a527_2021-04-13.sif -resume -profile fully_reproducible
```

# 2021-05-19 11:48:02

I ran this pipeline to the end, but the end results were not really what I
expected. The format is a binary and it's not clear what it represents. I think
I'd be much more comfortable with something like this:

https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226715

And then use our usual atlases to compute the tracts, or even the ENIGMA
pipeline. We could even go with Tracula
(https://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/Tracula), which
incorporates Freesurfer, but we might need to tweak the pipeline a bit to be in
Freesurfer space. And we'd need to run bedPostX, which is always a pain.

But maybe we could combine both pipelines? Start with 3dSlicer, run eddy_openmp
on that (with the new versions and all bells and whitles). Then, run Tracula to
use Freesurfer segmentations and avoid using T2s, and then use Camino's restore
for voxelwise work. Maybe we can use that instead of dtifit too, but we'll see.

For now, let's see if we can go all the way with this pipeline.

Actually, revising what I said before about the TractoFlow pipeline, it might
not be all that bad. Why don't we run it for 100 subjects in PNC, then 100 in
ours, and contrast that with the merged pipeline, and see how it goes? We should
probably use the same inputs though, so let's run DTIPrep first in both of them
just in case.

Let's start with 10 subjects, so it's somewhat manageable. I'll still BIDSize
them because it could be useful for fMRI later. Note that looking at the PNC
paper all images are in the AP direction, so no way to use the field maps for
blip up/down. I changed the config for dcm2bids accordingly. Also, remember that
I'll need to concatenate the DWI section before doing anything!

```bash
# ncrshell01
conda activate dcm2bids
cd /mnt/shaw/sudregp/PNC_BIDS;
base_dir=/mnt/NCR/sudregp/PNC/raw/66134/NeurodevelopmentalGenomics;
tmp_dir=/mnt/shaw/sudregp/tmp/;
for s in `cat batch1.txt`; do
    echo "Uncompressing ${s}";
    gz_fname=${s}_1.tar.gz;
    if [[ -e ${base_dir}/Study_version1/${gz_fname} ]]; then
        tar -zxf ${base_dir}/Study_version1/${gz_fname} -C $tmp_dir;
    else
        tar -zxf ${base_dir}/Study_version2/${gz_fname} -C $tmp_dir;
    fi;
    echo "BIDSizing ${s}";
    dcm2bids -d ${tmp_dir}/${s} -p ${s} -c ~/research_code/pnc2bids.json;
    rm -rf $tmp_dir/${s};
done;
```

Now, let's copy them to BW and start processing:

```bash
# ncrshell01
cd /mnt/shaw/sudregp/PNC_BIDS;
for s in `cat batch1.txt`; do
    echo "Copying ${s}";
    scp -qr sub-${s} helix:/scratch/sudregp/PNC_BIDS/;
done;
```

Before we run DTIPrep, remember to concatenate the runs:

```bash
# interactive
module load fsl/6.0.0

base_dir=/scratch/sudregp/PNC_BIDS
target_dir=/scratch/sudregp/PNC_DWI
for s in `cat ${base_dir}/batch1.txt`; do
    echo "Combining scans for ${s}";
    sdir=${target_dir}/${s};
    mkdir -p $sdir;
    cd ${base_dir}/sub-${s}/dwi;
    fslmerge -t ${sdir}/dwi sub-${s}_run-01_dwi sub-${s}_run-02_dwi;
    paste sub-${s}_run-01_dwi.bvec sub-${s}_run-02_dwi.bvec > ${sdir}/bvec;
    paste sub-${s}_run-01_dwi.bval sub-${s}_run-02_dwi.bval > ${sdir}/bval;
    cp ../anat/sub-${s}_T1w.nii.gz ${sdir}/t1.nii.gz
done;
```

That takes less than 5min per subject, so I might need to parallelize in the
future. We'll see. On the other hand, tractoflow can take hours, so we might as
well make it more SLURM friendly from the get-go.

Actually, I just noticed that DTIPrep changes the bvalues by a lot, so that
might get confusing. Let's use the original data there.

Well, PNC doesn't have T2, so that somewhat solves that issue.

This is what I'll put in my run wrapper file:

```bash
subj=$1
tractoflow_home=/data/NCR_SBRB/software/tractoflow/;
nextflow run ${tractoflow_home}/main.nf \
    --bids /scratch/sudregp/PNC_BIDS/ \
    --participants_label "${subj}" \
    --dti_shells "0 1010" --fodf_shells "0 1000 2000" \
    -with-singularity ${tractoflow_home}/tractoflow_2.2.1_b9a527_2021-04-13.sif \
    -resume -profile fully_reproducible
```

I changed to 1010 instead of 1000 (the default) because there was a 1005 bvalue.

```bash
#!/bin/sh
cd /data/NCR_SBRB/tractoflow;
for s in `cat ${base_dir}/batch1.txt`; do
    echo "bash /data/NCR_SBRB/tractoflow/run.sh $s" >> swarm.1;
done
swarm -g 10 -t 16 --job-name tractoflow --time 12:00:00 -f swarm.1 \
    -m nextflow/21.04.1,singularity/3.7.3 --partition norm --logdir trash
```

That's not working... I'd need to run it batching using MPI, and run multiple
subjects at the same time. It's possible, but not great.

But I also need to use the rootdir stick, otherwise BIDS interpret the two runs:

```bash
# make sure only subjects to be processed now are in the directory!
tractoflow_home=/data/NCR_SBRB/software/tractoflow/;
target_dir=/scratch/sudregp/PNC_DWI;
nextflow run ${tractoflow_home}/main.nf \
    --input ${target_dir}/ \
    --output_dir ${target_dir}/results \
    --dti_shells "0 1010" --fodf_shells "0 1000 2000" \
    -with-singularity ${tractoflow_home}/tractoflow_2.2.1_b9a527_2021-04-13.sif \
    -resume -profile fully_reproducible
```

# 2021-05-20 06:54:42

Today, I'm rethinking this approach. Based on the ABCD notes, they did use the
FDT pipeline, but used AtlasTrack for tract measures. The most time consuming
step of Tractoflow and the Tracula pipelines is bedpostX, so I wonder whether we
need the probabilistic tracking at all. Can we just do what ENIGMA does, and
what we were doing before with the averages?

AtlasTrack does use a probialistic approach, but that software is not available
anywhere...

Also, Tracula depends on a good Freesurfer segmentation. Sure, we'd need good
quality Freesurfer data to use it in the analysis to begin with. However, we
don't want to constraint the DTI analysis on Freesurfer analysis either. Maybe
just a good T1 would be enough to use the DTI data, because for most atlases
we'd need alignment to T1s.

Or, I can just stick with the FDT pipeline, which doesn't require any other
structural. Then we can do voxel analysis after DTI-TK coregistration, or use
the JHU atlases after registering to an FA map.

For the ENIGMA-suggested pipeline, we don't have the T2 in PNC. So, let's modify
it a bit.

# 2021-05-20 10:39:40

For TractoFlow, I ran 9 subjects in a 32machine in 14h:

```
Duration    : 14h 36m 3s
CPU hours   : 212.0
```

It took a maximum of 32Gb, and only in a few instances it took all 32 cores.
Given that we can only run this once at a time, it makes sense to do more
subjects for longer.

I'm just not sure we're going to go with probabilistic tracts for now. Let's
tackle the FDT pipeline first, and then we can decide.

Going back to DTIprep, we'll need to run it quickly first ro get the xml file,
which is based on the data (bvec and bval) and then change it to only run the
QC, and not do eddy or motion correction (which we'll do later). If the bvecs
were always the save this wouldn't be an issue (potentially), but this will not
be the case for our data, and who knows if the sequences changed in the middle
for other datasets. So, better safe then sorrow here.

I did run a text and the xml file we get by a quick run (aborted) and the entire
run is the same. So, let's modify it and get going.

```bash
timeout 5s ~/tmp/DTIPrepTools-1.2.9_rhel7-Linux/bin/DTIPrep --DWINrrdFile test.nrrd --numberOfThreads 4 -d -c -p test.xml
# let's do denoising
lnumber=`grep -n \"DENOISING_bCheck\" test.xml | cut -d":" -f 1`;
let lnumber=$lnumber+1
sed "${lnumber}s/No/Yes/" test.xml > test2.xml
# don't do eddy
lnumber=`grep -n \"EDDYMOTION_bCheck\" test2.xml | cut -d":" -f 1`;
let lnumber=$lnumber+1
sed "${lnumber}s/Yes/No/" test2.xml > test3.xml

# re-run using this new file
~/tmp/DTIPrepTools-1.2.9_rhel7-Linux/bin/DTIPrep --DWINrrdFile test.nrrd --numberOfThreads 4 -c -p test3.xml
```

But I cannot find the stupid Rician module, so let's just use ANTs version
before we convert to nrrd, and then just use DTIPrep to remove bad slices.

```bash
module load ANTs/2.3.2
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK
DenoiseImage -v 1 -d 4 -i dwi.nii.gz -o dwi_denoised.nii.gz -n Rician

/data/NCR_SBRB/software/DTIPrep-1.2.11/bin/DWIConvert \
    --conversionMode FSLToNrrd --inputVolume dwi_denoised.nii.gz \
    --inputBValues bval --inputBVectors bvec -o dwi.nrrd;
timeout 5s /data/NCR_SBRB/software/DTIPrepTools-1.2.9_rhel7-Linux/bin/DTIPrep \
    --DWINrrdFile dwi.nrrd --numberOfThreads $SLURM_CPUS_PER_TASK \
    -d -c -p ${s}.xml;
# don't do eddy
lnumber=`grep -n \"EDDYMOTION_bCheck\" ${s}.xml | cut -d":" -f 1`;
let lnumber=$lnumber+1
sed "${lnumber}s/Yes/No/" ${s}.xml > ${s}_noeddy.xml
# we also don't need to calculate the tensors or the gradient check, because
# it fails if we don't run eddy
/data/NCR_SBRB/software/DTIPrepTools-1.2.9_rhel7-Linux/bin/DTIPrep \
    --DWINrrdFile dwi.nrrd --numberOfThreads $SLURM_CPUS_PER_TASK \
    -c -p ${s}_noeddy.xml;
# save the original image (bval and bvec stay the same)
mv dwi.nii.gz dwi_orig.nii.gz;
/data/NCR_SBRB/software/DTIPrep-1.2.11/bin/DWIConvert \
    --inputVolume dwi_QCed.nrrd \
    --outputVolume dwi.nii.gz \
    --outputBVectors bvec \
    --outputBValues bval \
    --allowLossyConversion \
    --conversionMode NrrdToFSL;
```

This is not working. The bval and bvec files are very different, even when
nothing is done (only the QC). Maybe I could just grab which volumes are bad and
remove them manually afterwards?

I'll need to use one of our crappiest scans to see if some get recognized...
I'll also need to check how long ANTs DenoiseImage took. I left it running with
32 cores, and it was taking over 1h...

# 2021-05-21 06:54:46

In the end, it took almost 2h, on 32 cores. Maybe worth it? It used less than
2Gb the whole time. If we get Slicer to work, we can run it for 10 subjects and
check it.

In the meanwhile, I'll go ahead and convert the 10 subjects I've been playing
with, just so they'r ready when I' done with 3dSlicer.

```bash
# interactive
module load ANTs/2.3.2
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$SLURM_CPUS_PER_TASK

base_dir=/scratch/sudregp/PNC_BIDS
target_dir=/scratch/sudregp/PNC_DWI
for s in `cat ${base_dir}/batch1.txt`; do
    echo "Denoising ${s}";
    cd ${target_dir}/${s};
    DenoiseImage -v 1 -d 4 -i dwi.nii.gz -o dwi_denoised.nii.gz -n Rician
done;
```

Going back to our dataset, the following are missing 40 volumes: 

0741
1493
1999
2145

Let's try them through Slicer and see what we get...

It didn't detect anything. OK, quitting on DTIPrep. But I did find this:

 * https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddyqc/UsersGuide
 * https://www.sciencedirect.com/science/article/pii/S1053811918319451

So we can just run eddy instead. The main question is whether the initial noise
removal makes a difference. Given that it takes 2h per subject (using 32 cores),
it'd have to make a huge difference. 

For now, let's revise what we did for FDT PNC a couple years ago, make sure the
initial data looks good (in terms of directions), check all the brain mask QCs,
and other QC that was done. Then we can try running QUAD and seeing what we get.

Then, repeat the same for our own cohort. I'll need to check which images we
started with though.


# 2021-06-28 12:39:56

Going back to this, I'll give qsiprep a chance. I'll process the 10 subjects I
started with first, but in the meanwhile I can also BIDSize another batch, say
about 100. The only oplan is that it's taking forever to install dcm2bids in BW,
so I'll just keep BIDSizing in ncrshell01:

```bash
# ncrshell01
conda activate dcm2bids
cd /mnt/shaw/sudregp/PNC_BIDS;
base_dir=/mnt/NCR/sudregp/PNC/raw/66134/NeurodevelopmentalGenomics;
tmp_dir=/mnt/shaw/sudregp/tmp/;
for s in `cat batch2.txt`; do
    echo "Uncompressing ${s}";
    gz_fname=${s}_1.tar.gz;
    if [[ -e ${base_dir}/Study_version1/${gz_fname} ]]; then
        tar -zxf ${base_dir}/Study_version1/${gz_fname} -C $tmp_dir;
    else
        tar -zxf ${base_dir}/Study_version2/${gz_fname} -C $tmp_dir;
    fi;
    echo "BIDSizing ${s}";
    dcm2bids -d ${tmp_dir}/${s} -p ${s} -c ~/research_code/pnc2bids.json;
    rm -rf $tmp_dir/${s};
done;
```

qsiprep in BW, using their generic pipeline, took only 1.5h. It used less than
15Gb of memory, and although it had a couple spikes, most of the time it stayed
under 18 threads, even though I gave it all 32.

Let me try another subject just to check, specifying less memory and threads.

Analyzing that first output (sub-600053476714), a few things:
 * need to combine both runs after denoising (I thought it was the default?)
 * what derivatives are being computed?
 * Freesurfer?

Maybe this will do it.

qsiprep ./ ../PNC_qsiprep_outputs/ participant \
    --participant_label sub-600039665619 -w /lscratch/$SLURM_JOB_ID \
    --notrack --nthreads 16 --mem_mb 15000 --stop-on-first-crash \
    --output-resolution 1.2 --recon-input ../PNC_qsiprep_outputs/ \
    --recon_spec mrtrix_singleshell_ss3t --do-reconall \
    --denoise_before_combining --combine_all_dwis \
    --recon_input ../PNC_qsiprep_outputs/qsiprep \
    --recon_spec mrtrix_singleshell_ss3t


# 2021-06-29 06:11:16

This time it took **MUCH** longer, as it should because of freesurfer. But it
looks like the freesurfer subject inside PNC_qsiprep_outputs is just
$SUBJECTS_DIR, so maybe I could symlink it and not have to run Freesurfer? It
also creates another folder in the same level, qsirecon, which likely has some
of the goodies.

```
(base) [sudregp@cn0990 PNC_BIDS]$ cd ../PNC_qsiprep_outputs/
(base) [sudregp@cn0990 PNC_qsiprep_outputs]$ ls
freesurfer  qsiprep  qsirecon
```

However, 600053476714 has Freesurfer recon-all issues. And it looks like
600039665619 died because "Exception: No MNI to T1w transform found in
anatomical directory". So, nothing to show yet...

Let me see if I can use some pre-made freesurfer then. I'll play with
sub-600062084650 and sub-600018902293 because they look like they've finished
freesurfer from Luke's fmriprep:

```bash
# will need to use symlinks in the future!
cd /data/NCR_SBRB/PNC_BIDS
cp -r /data/NCR_SBRB/PNC_fmriprep_fs_20.2.0/sub-600062084650/freesurfer /data/NCR_SBRB/PNC_qsiprep_outputs/;
cp -r /data/NCR_SBRB/PNC_fmriprep_fs_20.2.0/sub-600018902293/freesurfer/sub-600018902293 /data/NCR_SBRB/PNC_qsiprep_outputs/freesurfer/;

s=sub-600018902293;
qsiprep ./ ../PNC_qsiprep_outputs/ participant \
    --participant_label $s -w /lscratch/$SLURM_JOB_ID \
    --notrack --nthreads 16 --mem_mb 15000 --stop-on-first-crash \
    --output-resolution 1.2 --recon-input ../PNC_qsiprep_outputs/ \
    --recon_spec mrtrix_singleshell_ss3t \
    --denoise_before_combining --combine_all_dwis \
    --skip_bids_validation --unringing-method mrdegibbs \
    --force-spatial-normalization --template MNI152NLin2009cAsym
```

Started these at 6:48am. They finished only at 4:44pm. So, we're looking at a
10-11h process for reconstruction.

Someone on NeuroStarts said qsiprep doesn't need Freesurfer. Is that correct?
Let me try some other subject.

```bash
cd /data/NCR_SBRB/PNC_BIDS
s=sub-600138617917;
qsiprep ./ ../PNC_qsiprep_outputs/ participant \
    --participant_label $s -w /lscratch/$SLURM_JOB_ID \
    --notrack --nthreads 16 --mem_mb 15000 --stop-on-first-crash \
    --output-resolution 1.2 --recon-input ../PNC_qsiprep_outputs/ \
    --recon_spec mrtrix_singleshell_ss3t \
    --denoise_before_combining --combine_all_dwis \
    --skip_bids_validation --unringing-method mrdegibbs \
    --force-spatial-normalization --template MNI152NLin2009cAsym
```

Started that at 10:37am. The difference from the two above is that it doesn't
have a freesurfer directory for this particular subject.

It doesn't looks like tensors are being computed in the default recon, but
here's a way to do it:

```bash
# bw
conda activate mrtrix3

dwi2tensor -fslgrad /data/NCR_SBRB/PNC_qsiprep_outputs/qsiprep/sub-600039665619/dwi/sub-600039665619_space-T1w_desc-preproc_dwi.bvec /data/NCR_SBRB/PNC_qsiprep_outputs/qsiprep/sub-600039665619/dwi/sub-600039665619_space-T1w_desc-preproc_dwi.bval -nthreads 4 /data/NCR_SBRB/PNC_qsiprep_outputs/qsiprep/sub-600039665619/dwi/sub-600039665619_space-T1w_desc-preproc_dwi.nii.gz test.nii.gz
tensor2metric test.nii.gz -fa test_fa.nii.gz -ad test_ad.nii.gz -rd test_rd.nii.gz
```

And then we could run TBSS in it to get skeletons, or use the transformation
computed in the anat folder:

```
(mrtrix3) [sudregp@cn0990 anat]$ pwd

/data/NCR_SBRB/PNC_qsiprep_outputs/qsiprep/sub-600018902293/anat
(mrtrix3) [sudregp@cn0990 anat]$ ls

sub-600018902293_desc-brain_mask.nii.gz
sub-600018902293_desc-preproc_T1w.nii.gz
sub-600018902293_dseg.nii.gz
sub-600018902293_from-MNI152NLin2009cAsym_to-T1w_mode-image_xfm.h5
sub-600018902293_from-orig_to-T1w_mode-image_xfm.txt
sub-600018902293_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5
sub-600018902293_label-CSF_probseg.nii.gz
sub-600018902293_label-GM_probseg.nii.gz
sub-600018902293_label-WM_probseg.nii.gz
sub-600018902293_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz
sub-600018902293_space-MNI152NLin2009cAsym_desc-preproc_T1w.nii.gz
sub-600018902293_space-MNI152NLin2009cAsym_dseg.nii.gz
sub-600018902293_space-MNI152NLin2009cAsym_label-CSF_probseg.nii.gz
sub-600018902293_space-MNI152NLin2009cAsym_label-GM_probseg.nii.gz
sub-600018902293_space-MNI152NLin2009cAsym_label-WM_probseg.nii.gz
```

to put the tensor into the template space, or use DTI-TK... the sky is the limit.

So, for now I imagine I can just fire up the pipeline without recon for all
subjects, and then try the recon-only option to create the connectivity, while
running the tensors as well?

# 2021-06-30 05:53:11

Let's fire up the swarm for the first batch, preproc only.

```bash
cd /scratch/sudregp/;
s=600009963128;
bash ~/research_code/qsiprep_preproc_wrapper.sh ${s} \
    /data/NCR_SBRB/PNC_BIDS /data/NCR_SBRB/PNC_qsiprep_outputs
```

```bash
# bw
cd /scratch/sudregp

main_bids=/data/NCR_SBRB/PNC_BIDS/;
out_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/;
bname=batch5;

rm -rf swarm.$bname; 
for m in `cat ${bname}.txt`; do
    echo "bash ~/research_code/qsiprep_wrapper.sh $m $main_bids $out_dir " >> swarm.$bname;
done
swarm -g 20 -t 16 --logdir trash_qsiprep --gres=lscratch:20 --time 4:00:00 \
    -f swarm.$bname --partition quick,norm --job-name $bname -m qsiprep/0.8.0
```

In the meanwhile, let's figure out the output from the connectivity analysis.
Following this great book:

https://andysbrainbook.readthedocs.io/en/latest/MRtrix/MRtrix_Course/MRtrix_07_Streamlines.html

tckgen generates the ifod2_tck file. Then, we run tcksift2 to conter-balance the
overfitting, and generate siftweights_ifod2.csv. As expected, the last file
that's created is the .mat, which is a binary. Probably Matlab? But we can
likely also get the individual CSVs by doing something like this, using the
other MIFs in the directory:

```bash
tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in sift_1M.txt tracks_10M.tck sub-CON02_parcels.mif sub-CON02_parcels.csv -out_assignment assignments_sub-CON02_parcels.csv
```

But that's not necessary:

```r
library(R.matlab)
df = readMat('/data/NCR_SBRB/PNC_qsiprep_outputs/qsirecon/sub-600018902293/dwi/sub-60
    0018902293_space-T1w_desc-preproc_space-T1w_dhollanderconnectome.mat')
```

That gives me a list like:

```
[...]
[49] "gordon333.region.ids"                                     
[50] "gordon333.region.labels"                                  
[51] "gordon333.sift.invnodevol.radius2.count.connectivity"     
[52] "gordon333.radius2.meanlength.connectivity"                
[53] "gordon333.radius2.count.connectivity"                     
[54] "gordon333.sift.radius2.count.connectivity"                
[55] "schaefer200x17.region.ids"                                
[56] "schaefer200x17.region.labels"                             
[57] "schaefer200x17.sift.invnodevol.radius2.count.connectivity"
[58] "schaefer200x17.radius2.meanlength.connectivity"           
[59] "schaefer200x17.radius2.count.connectivity"                
[60] "schaefer200x17.sift.radius2.count.connectivity"           
[61] "brainnetome246.region.ids"                                
[62] "brainnetome246.region.labels"                             
[63] "brainnetome246.sift.invnodevol.radius2.count.connectivity"
[64] "brainnetome246.radius2.meanlength.connectivity"           
[65] "brainnetome246.radius2.count.connectivity"                
[66] "brainnetome246.sift.radius2.count.connectivity"    
```

So, many things to play with there. We have the mean length of the tracts and
tract count. invnodevol scales each contribution to the connectome edge by the
inverse of the two node volumes. radius2 means that the search radius was 2mm.
sift means it uses the sift weights to counter overfitting.

I also checked and there isn't anyone with more than 2 DTIs in PNC.

## ABCD

Not sure I'll be able to reproduce the ABCD metrics. Not only because their
pipeline was weird, but the tabular metrics don't see like anything I can
identify. It's always one of the Destrieux ROIs against the entire hemisphere?

Let's work on the tensor script:

```bash
cd /lscratch/$SLURM_JOBID
cp -r /data/NCR_SBRB/PNC_qsiprep_outputs/sub-600009963128 .
cd sub-600009963128/
conda activate mrtrix3
dwi2tensor -fslgrad sub-600009963128_space-T1w_desc-preproc_dwi.bvec sub-600009963128_space-T1w_desc-preproc_dwi.bval -nthreads $SLURM_CPUS_PER_TASK sub-600009963128_space-T1w_desc-preproc_dwi.nii.gz tensors.nii.gz
tensor2metric tensors.nii.gz -fa dti_FA.nii.gz -ad dti_AD.nii.gz -rd dti_RD.nii.gz

module load fsl/6.0.4/fsl
mkdir tbss
cp dti_FA.nii.gz tbss
cd tbss
tbss_1_preproc dti_FA.nii.gz;
# this does the same thing as "tbss_2_reg -T" without calling SLURM
cd FA; fsl_reg dti_FA_FA target dti_FA_FA_to_target -e -FA; cd ..;
tbss_3_postreg -S;
tbss_4_prestats 0.2;
for m in AD RD; do
    mkdir $m;
    cp ../dti_${m}.nii.gz ${m}/dti_FA.nii.gz;
    tbss_non_FA ${m};
done
# the results to run the analysis are all in tbss/stats/all_??_skeletonised.nii.gz
# we should probably just copy the entire tbss directory though

# now, for the JHU tracts, the dti_FA and other property files are already in the correct space, so it's just a matter of averaging over the mask. See script in note 007 for example.
```

Let's check which ones didn't finish qsiprep:

```bash
cd /scratch/sudregp/;
qsiprep_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/;
suf='space-T1w_desc-preproc_dwi';
rm -rf qsiprep_redo.txt
for s in `cat batch?.txt`; do
    if [ ! -e ${qsiprep_dir}/sub-${s}/dwi/sub-${s}_${suf}.nii.gz ]; then
        echo $s >> qsiprep_redo.txt;
    fi;
done
```

Only 26 did not finish. Let's erase their folders just in case and re-run with a
bigger wall time:

```bash
# bw
cd /scratch/sudregp

main_bids=/data/NCR_SBRB/PNC_BIDS/;
out_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/;
bname=qsiprep_redo;

rm -rf swarm.$bname; 
for m in `cat ${bname}.txt`; do
    rm -rf ${out_dir}/sub-${m};
    echo "bash ~/research_code/qsiprep_wrapper.sh $m $main_bids $out_dir " >> swarm.$bname;
done
swarm -g 20 -t 16 --logdir trash_qsiprep --gres=lscratch:20 --time 8:00:00 \
    -f swarm.$bname --partition norm --job-name $bname -m qsiprep/0.8.0
```

We're at a point where we can run tbss. Let's run for everyone, knowing that
some will fail:

```bash
# bw
cd /scratch/sudregp

conda activate mrtrix3;
out_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/;
bname=batch5;

rm -rf swarm.$bname; 
for m in `cat ${bname}.txt`; do
    rm -rf ${out_dir}/tbss/sub-${m};
    echo "bash ~/research_code/tbss_pnc_wrapper.sh $m $out_dir" >> swarm.$bname;
done
swarm -g 8 -t 2 -b 10 --logdir trash_tbss --gres=lscratch:20 --time 20:00 \
    -f swarm.$bname --partition quick,norm --job-name $bname -m fsl/6.0.4/fsl
```

Had to run a few subjects with 2h because they were taking forever.

I then did a quick check of the DEC maps for a few random subjects, and they
look fine. To do that, we need a V1 image, which I have not produced or saved.
In fact, I haven't been saving the tensors either, but they calculate quite fast
from the cleaned data. For the DEC, I followed this:

https://community.mrtrix.org/t/fa-color-maps/992/4

```bash
cd /lscratch/$SLURM_JOBID
cp -r /data/NCR_SBRB/PNC_qsiprep_outputs/sub-600009963128 .
cd sub-600009963128/
conda activate mrtrix3
dwi2tensor -fslgrad sub-600009963128_space-T1w_desc-preproc_dwi.bvec \
    sub-600009963128_space-T1w_desc-preproc_dwi.bval \
    -nthreads $SLURM_CPUS_PER_TASK \
    sub-600009963128_space-T1w_desc-preproc_dwi.nii.gz tensors.nii.gz
cd dwi
tensor2metric tensors.nii.gz -fa dti_FA.nii.gz -ad dti_AD.nii.gz \
    -rd dti_RD.nii.gz -vector dti_V1.nii.gz -force
mrcalc dti_V1.nii.gz -abs fac.nii.gz
module load afni
fat_proc_decmap -in_fa dti_FA.nii.gz -in_v1 fac.nii.gz -fa_thr .2 -prefix DEC
# then look at the pngs in the QC folder
```

Now, we just need to calculate the JHU tracts and labels. If we want our usual
11 tracts from DTI-TK, we'll need to run that entire pipeline in the tensors.
It's doable, but it will take a while.

I talked ot Marine and we agreed that it's best to wait and just play with JHU
for now. Before I do that, I'll organize everything in
/data/NCR_SBRB/PNC_qsiprep_outputs/ inside qsiprep, on the same level as tbss. I
then started rsyncing everything back to Shaw.

## Averaging JHU tracts and labels

```bash
# sinteractive
module load afni
module load fsl

tbss_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/tbss;
atlas=$FSLDIR/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz;
weighted_tracts=pnc_jhu_tracts.csv;

# copy everything to /lscratch to make it faster (need 300Gb)
cd /lscratch/$SLURM_JOBID/;
cp -r $tbss_dir .;
wrk_dir=$RANDOM;
mkdir -p $wrk_dir; cd $wrk_dir;
imcp $atlas ./atlas.nii.gz
row="id";
for t in `seq 1 20`; do
    for m in fa ad rd; do
        row=${row}','${m}_${t};
    done
done
echo $row > $weighted_tracts;
for m in `cat /scratch/sudregp/ready.txt`; do
    row="${m}";
    froot=../tbss/sub-$m/stats/all;
    echo $m;
    for t in `seq 1 20`; do
        3dcalc -a atlas.nii.gz -expr "amongst(a, $t)" -prefix mask.nii \
            -overwrite 2>/dev/null &&
        fa=`3dmaskave -q -mask mask.nii ${froot}_FA.nii.gz 2>/dev/null`;
        ad=`3dmaskave -q -mask mask.nii ${froot}_AD.nii.gz 2>/dev/null`;
        rd=`3dmaskave -q -mask mask.nii ${froot}_RD.nii.gz 2>/dev/null`;
        row=${row}','${fa}','${ad}','${rd};
    done
    echo $row >> $weighted_tracts;
done
```

And re-using the copied files from above, we can run the same thing for JHU
labels:


```bash
atlas=$FSLDIR/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz;
weighted_tracts=pnc_jhu_labels.csv;

# copy everything to /lscratch to make it faster (need 300Gb)
cd /lscratch/$SLURM_JOBID/;
wrk_dir=$RANDOM;
mkdir -p $wrk_dir; cd $wrk_dir;
imcp $atlas ./atlas.nii.gz
row="id";
for t in `seq 1 48`; do
    for m in fa ad rd; do
        row=${row}','${m}_${t};
    done
done
echo $row > $weighted_tracts;
for m in `cat /scratch/sudregp/ready.txt`; do
    row="${m}";
    froot=../tbss/sub-$m/stats/all;
    echo $m;
    for t in `seq 1 48`; do
        3dcalc -a atlas.nii.gz -expr "amongst(a, $t)" -prefix mask.nii \
            -overwrite 2>/dev/null &&
        fa=`3dmaskave -q -mask mask.nii ${froot}_FA.nii.gz 2>/dev/null`;
        ad=`3dmaskave -q -mask mask.nii ${froot}_AD.nii.gz 2>/dev/null`;
        rd=`3dmaskave -q -mask mask.nii ${froot}_RD.nii.gz 2>/dev/null`;
        row=${row}','${fa}','${ad}','${rd};
    done
    echo $row >> $weighted_tracts;
done
```

This takes a LONG time, even grabbing from lscratch. May need to parallelize it
for ABCD, or even for this if the interactive machine times out.

## Recon

I copied the 3 subjects I ran fully in qsirecon, but from now I'll only get the
matlab file and the figures for the html report to save on size. Each directory
is 15G, so we'd be looking at lots of Terabytes if we save everything. Testing
the wrapper script now to make sure we copy the correct files in the end.

For now, let's collect the QC variables:

```bash
qsiprep_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/qsiprep;
out_fname=pnc_qsiprep_QC.csv;

cd /scratch/sudregp/;
s1=`head -n 1 ready.txt`;
head -n 1 ${qsiprep_dir}/sub-${s1}/dwi/sub-${s1}_desc-ImageQC_dwi.csv > $out_fname;
for s in `cat ready.txt`; do
    tail -n +2 ${qsiprep_dir}/sub-${s}/dwi/sub-${s}_desc-ImageQC_dwi.csv >> $out_fname;
done
```




# TODO
  * compile JHU tracts
  * run recon

# useful links
 * https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226715
 * https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24691
