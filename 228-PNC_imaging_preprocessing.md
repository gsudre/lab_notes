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

For the ENIGMA-suggested pipeline, we don't have the T2 in PNC. So, let's modify
it a bit.


Now, for the ENIGMA-suggested pipeline, using DTIPrep:










    /data/NCR_SBRB/software/DTIPrep-1.2.11/bin/DWIConvert \
        --conversionMode FSLToNrrd --inputVolume dwi.nii.gz \
        --inputBValues bval --inputBVectors bvec -o dwi.nrrd;
    /data/NCR_SBRB/software/DTIPrepTools-1.2.9_rhel7-Linux/bin/DTIPrep \
        --DWINrrdFile dwi.nrrd --numberOfThreads $SLURM_CPUS_PER_TASK \
        -d -c -p ${s}.xml;
    mv dwi.nii.gz dwi_orig.nii.gz;
    mv bval bval_orig;
    mv bvec bvec_orig;
    /data/NCR_SBRB/software/DTIPrep-1.2.11/bin/DWIConvert \
        --inputVolume dwi_QCed.nrrd \
        --outputVolume dwi.nii.gz \
        --outputBVectors bvec \
        --outputBValues bval \
        --allowLossyConversion \
        --conversionMode NrrdToFSL;





# TODO
 * Maybe add the GAD filtering from here? https://discourse.slicer.org/t/using-slicer-and-slicer-modules-from-command-line/8162/2
  

# useful links
 * https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226715
 * https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24691
