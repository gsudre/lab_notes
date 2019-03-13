# 2019-03-08 19:03:45

We might end up not even doing this in the future, but since I had to do some
more DTI-TK processing for the tentative heritability change work, let's write
the pipeline here:

```bash
dataDir=/Volumes/Shaw/dti_robust_tsa/analysis_may2017/
maskidDir=/Volumes/Shaw/MR_data_by_maskID/
maskidFile=~/tmp/all_maskids.txt
while read m; do
    dtitkFolder=${maskidDir}/${m}/edti_proc/edti_DMC_DR_R1_SAVE_DTITK/;
    if [ ! -d $dtitkFolder ]; then
        dtitkFolder=${maskidDir}/${m}/edti_proc/edti_DMC_R1_SAVE_DTITK/;
        tensorFile=edti_DMC_R1_tensor.nii;
    else
        tensorFile=edti_DMC_DR_R1_tensor.nii;
    fi;
    if [ ! -d $dtitkFolder ]; then
        echo $dtitkFolder "does not exist. Skipping...";
    else
         cp -v ${dtitkFolder}/${tensorFile} ${dataDir}/${m}_tensor.nii;
    fi;
done < $maskidFile
```

With all the data now copied, and knowing from previous analysis that there
aren't any outliers, let's copy the whole thing to the cluster and run the rest
of the analysis there. Let's not worry about the resampling for specific
tensors, and just do everybody.

```bash
for m in `cat ~/tmp/all_maskids.txt`; do
    echo $m;
    scp -q $dataDir/${m}_tensor.nii helix:/scratch/sudregp/dtitk/;
done
scp -q $dataDir/0121_tensor.nii helix:/scratch/sudregp/dtitk/;  # just as a template
```

Then, in Biowulf:

```bash
cd /scratch/sudregp/dtitk;
ls -1 *tensor.nii > subj_tensor.txt;
while read f; do TVtool -in ${f} -spd -out ${f}; done < subj_tensor.txt;
while read f; do TVResample -in $f -target 0121_tensor.nii; done < subj_tensor.txt;
sed 's/_tensor.nii//g' subj_tensor.txt > subj_ids.txt;
while read m; do echo bash ~/research_code/dti/run_tsa.sh $m >> swarm.tsa; done < subj_ids.txt
swarm -f swarm.tsa -t 1 --logdir trash_tsa --job-name dtitk --partition quick
```

Then we check who failed:

```bash
while read s; do
    if [ ! -e ${s}_tensor_diffeo.nii.gz ]; then
        echo $s;
    fi;
done < subj_ids.txt
```

I have instructions for sampling in other notes, but it goes generally like
this, after copying it back locally:

```bash
for m in `cat ~/tmp/all_maskids.txt`; do
    echo $m;
    scp -q helix:/scratch/sudregp/dtitk/${m}* ${dataDir}/;
done
```

```bash
sed "s/$/_tensor_diffeo\.nii/g" ids224.txt > tensors224.txt
cd /Volumes/Shaw/dti_robust_tsa/analysis_may2017/
export DTITK_ROOT=/Applications/dtitk-2.3.3-Darwin-x86_64/
/Applications/dtitk-2.3.3-Darwin-x86_64/scripts/tsa_sampling ~/data/heritability_change/tensors224.txt ../ixi_aging_template_v3.0/tsa/ mean
# edit them first
python ~/research_code/lab_mgmt/convert_dti_sampling.py
Rscript ~/research_code/dti/compile_tract_table.R
```

Or, in Linux because shared drive access in the MAc sometimes is painfully slow:

```bash
cd ~/tmp
for m in `cat assoc3`; do echo ${m}_tensor_diffeo.nii.gz >> ready.txt; done
cd /mnt/shaw/sudregp/dti_robust_tsa/analysis_may2017/
export DTITK_ROOT=/usr/local/neuro/dti-tk/dtitk-2.3.1-Linux-x86_64
${DTITK_ROOT}/scripts/tsa_sampling ~/tmp/ready.txt ../ixi_aging_template_v3.0/tsa/ mean
```