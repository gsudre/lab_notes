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

Then, we swarm it in the cluster:

```bash
cd ~/tmp
awk '{ print $1 }' Results\ 1.txt | sort | uniq | head -n -1 > subjs.txt;
for s in `cat subjs.txt`; do
    echo "mriqc /scratch/sudregp/NCR_BIDS /scratch/sudregp/mriqc_output participant --participant_label ${s} -m T1w -w /scratch/sudregp/mriqc_work" >> swarm.mriqc;
done
swarm -g 8 -f swarm.mriqc --job-name mriqc --time 4:00:00 --logdir trash_mriqc -m mriqc --partition quick --gres=lscratch:40
```

<!-- And then we need to collect all results:

```bash
module load singularity
export SINGULARITY_CACHEDIR=/data/sudregp/singularity/
singularity exec -B /scratch/sudregp:/mnt docker://poldracklab/mriqc:latest mriqc /mnt/BIDS /mnt/mriqc_output group --no-sub -w /mnt/mriqc_work -m T1w -->


# TODO
* Run MRIQC on all our data
* Run Tonya's script on all our data
* Convert resting to BIDS
* Run MRIQC on resting