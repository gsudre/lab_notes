# 2019-02-08 14:36:01

If I'm converting everything to .nii.gz to run Tonya's script on our data, why
not convert all MPRAGE to BIDS and run MRIQC as well? We can even do the same
for our resting state...

First, let's do the anatomicals. We start with a file that has subj code (from
Labmatrix) in one columns, and mask id in the other:

```bash
cd /Volumes/Shaw/NCR_BIDS
# installed jo through homebrew
jo -p "Name"="Neurobehavioral Clinical Research (NCR) Section neuroimaging database" "BIDSVersion"="1.0.2" >> dataset_description.json;
while read line; do
    s=`echo $line | awk '{ print $1 }'`;
    m=`echo $line | awk '{ print $2 }'`;
    if [ ! -d sub-${s} ]; then
        mkdir sub-${s};
    fi;
    mkdir sub-${s}/ses-${m};
    mkdir sub-${s}/ses-${m}/anat;
    dcm2niix_afni -o sub-${s}/ses-${m}/anat/ -z y -f sub-${s}_ses-${m}_T1w /Volumes/Shaw/best_mprages/${m}/;
done < ~/Downloads/Results\ 1.txt
```

Tonya's scripts takes single .nii.gz files, or a directory full of them. So,
let's create symbolic links so we don't duplicate the data:

```bash
cd symlinks
while read line; do
    s=`echo $line | awk '{ print $1 }'`;
    m=`echo $line | awk '{ print $2 }'`;
    ln -s /Volumes/Shaw/NCR_BIDS/sub-${s}/ses-${m}/anat/sub-${s}_ses-${m}_T1w.nii.gz .;
done < ~/Downloads/Results\ 1.txt
```

I then just ran it like this:

```bash
conda activate mri-qa
cd ~/Downloads/mri-qa-python-master
for s in `/bin/ls -1 symlinks/`; do
    python3 tw-qa.py -i symlinks/${s};
done
```

And I copied the results for 1970 mask ids to ~/data/tonya_results.txt.

For MRIQC, we swarm it in the cluster:

```bash
cd ~/tmp
awk '{ print $1 }' Results\ 1.txt | sort | uniq | head -n -1 > subjs.txt;
for s in `cat subjs.txt`; do
    echo "mriqc /scratch/sudregp/NCR_BIDS /scratch/sudregp/mriqc_output participant --participant_label ${s} -m T1w -w /scratch/sudregp/mriqc_work" >> swarm.mriqc;
done
swarm -g 8 -f swarm.mriqc --job-name mriqc --time 4:00:00 --logdir trash_mriqc -m mriqc --partition quick --gres=lscratch:40
```

# 2019-02-21 11:49:52

Biowulf made a huge deal about my big swarms. So, I'll try to run it locally to
see what I can get:

```bash
docker run -it --rm -v /Volumes/Shaw/NCR_BIDS/:/data:ro \
    -v ~/data/mriqc_output/:/out poldracklab/mriqc:latest /data /out \
    participant --participant_label 99892 -m T1w
```

That worked, so now I'll just go ahead and make a different call for each
subject. I know that's not the best way to call it, but at least it gets split
into small bits that I can simply go back to later if stopped in the middle.

```bash
for s in `cat ~/tmp/bids_ids.txt`; do
    docker run -it --rm -v /Volumes/Shaw/NCR_BIDS/:/data:ro \
        -v ~/data/mriqc_output/:/out poldracklab/mriqc:latest /data /out \
        participant --participant_label $s -m T1w;
done
```

# TODO
* Run MRIQC on all our data
* Run Tonya's script on all our data
* Convert resting to BIDS
* Run MRIQC on resting