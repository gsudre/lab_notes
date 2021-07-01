# 2021-07-01 06:43:31

For our data, we'll have a few extra wrinkles when creating the BIDS directory
for DTI. Talking with Marine (Teams, 5/19/21, 1:50PM), the idea would be to use
the re-acquired volumes when possible, replacing the old ones. Second best would
be using everything. 

So, why not create up to 3 sessions for each mask id? Session 1 is the vanilla
60 or 80, which we'll process similarly to PNC. Session 2 includes the 99 run.
Session 3 has the 99 run replacing the original volumes. Session 3 would only be
created after qsiprep, by removing the re-acquired volumes from the eddy output.
Not sure how much difference this will make, or even if it's feasible, because
it'd only work if all volumes are processed independently.















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

```bash
# bw
cd /scratch/sudregp

conda activate mrtrix3;
out_dir=/data/NCR_SBRB/PNC_qsiprep_outputs/;
bname=tbss_redo2;

rm -rf swarm.$bname; 
for m in `cat ${bname}.txt`; do
    rm -rf ${out_dir}/tbss/sub-${m};
    echo "bash ~/research_code/tbss_pnc_wrapper.sh $m $out_dir" >> swarm.$bname;
done
swarm -g 8 -t 2 --logdir trash_tbss --gres=lscratch:20 --time 2:00:00 \
    -f swarm.$bname --partition quick,norm --job-name $bname -m fsl/6.0.4/fsl
```


# TODO
  * write tensor fit and JHU tracts pipeline
  * review html output of first batch


# useful links
 * https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226715
 * https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24691
