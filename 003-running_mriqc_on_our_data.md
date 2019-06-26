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

# 2019-03-25 13:49:49

OK, let's try converting the resting data to BIDS now. We assume the anatomical
folder for all subjects we have processed in AFNI have already been created, so
now all that remains to be done is create the functional directory. Based on the
IQMs, it made sense to grab all the resting DICOMs instead of just the ones we
processed. We can always just take the values that we want later.

# 2019-03-27 13:35:20

Well, I didn't really like how I had split the data per subject, then mask id.
We should keep that information only in Labmatrix. So, let's reparse and re-run
everything using just mask ids. This way we can also parse the resting state at
the same time:

```bash
net_dir=/Volumes/Shaw
ls -1 $net_dir/MR_data_by_maskid/ | sed 's/\///' > ~/tmp/maskids.txt
cd $net_dir/NCR_BIDS
# installed jo through homebrew
jo -p "Name"="Neurobehavioral Clinical Research (NCR) Section neuroimaging database" "BIDSVersion"="1.0.2" >> dataset_description.json;
while read m; do
    echo $m;
    mkdir sub-${m};

    # find name of date folders
    ls -1 $net_dir/MR_data_by_maskid/${m}/ | grep -e ^20 > ~/tmp/date_dirs;

    # for each date folder, check for mprages
    mkdir sub-${m}/anat;
    cnt=1
    while read d; do
        grep -i -e rage_ -e mprage $net_dir/MR_data_by_maskid/${m}/${d}/*README* > ~/tmp/mprage;
        awk '{for(i=1;i<=NF;i++) {if ($i ~ /Series/) print $i}}' ~/tmp/mprage | sed "s/Series://g" > ~/tmp/mprage_clean
        while read line; do
            mr_dir=`echo $line | sed "s/,//g"`;
            dcm2niix_afni -o sub-${m}/anat/ -z y -f sub-${m}_run-${cnt}_T1w ${net_dir}/MR_data_by_maskid/${m}/${d}/${mr_dir}/;
            let cnt=$cnt+1;
        done < ~/tmp/mprage_clean;
    done < ~/tmp/date_dirs;

    # for each date folder, check for resting scans
    mkdir sub-${m}/func;
    cnt=1
    while read d; do
        grep rest $net_dir/MR_data_by_maskid/${m}/${d}/*README* > ~/tmp/restall;
        # remove the movie runs for now
        grep -v movie ~/tmp/restall > ~/tmp/rest;
        awk '{for(i=1;i<=NF;i++) {if ($i ~ /Series/) print $i}}' ~/tmp/rest | sed "s/Series://g" > ~/tmp/rest_clean
        while read line; do
            mr_dir=`echo $line | sed "s/,//g"`;
            dcm2niix_afni -o sub-${m}/func/ -z y -f sub-${m}_task-rest_run-${cnt}_bold ${net_dir}/MR_data_by_maskid/${m}/${d}/${mr_dir}/;
            let cnt=$cnt+1;
        done < ~/tmp/rest_clean;
    done < ~/tmp/date_dirs;
done < ~/tmp/maskids.txt
```

```bash
for s in `cat ~/tmp/maskids.txt`; do
    docker run -it --rm -v /Volumes/Shaw/NCR_BIDS/:/data:ro \
        -v ~/data/mriqc_output/:/out poldracklab/mriqc:latest /data /out \
        participant --participant_label $s --no-sub;
done
```

# 2019-05-02 10:57:35

It's taking a long time to do it on my laptop, and also it's slowing things down
a ton. So, let's try to run it in Biowulf again. First, create the wrapper. I
limited it there to 4 cores and 8Gb, and one mask id ran for about 30min. So, we
can do a few mask ids like this:

```bash
ls -1 /scratch/sudregp/NCR_BIDS/ | sed "s/sub\-//" > ~/tmp/bids.txt;

cd ~/data/mriqc
for m in `head ~/tmp/bids.txt`; do
    echo "bash ~/research_code/mriqc_wrapper.sh $m" >> swarm.bids;
done
swarm -g 10 -t 4 --logdir trash_mriqc --gres=lscratch:10 --time 45:00 -f swarm.bids \
    --partition quick --job-name mriqc -m mriqc
```

When running in interactive node I had to do:

```bash
module purge
unset LD_LIBRARY_PATH
```

So, I'm not sure if I'll need to do something like that for swarms.

That worked well, and I got results in ~/data/mriqc_output. So, let's go ahead
and fire it up for everybody!

That seems to have worked... I just need to carefully check it later for any
errors, but the results have been copied to shaw/mriqc_output.

# TODO
* Run MRIQC on all our data (collect it)
* Run Tonya's script on all our data
