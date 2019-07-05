# 2019-05-30 15:23:00

I'll start following Luke's past experiences and the BW recommendations.
Starting with fmriprep, we do something like this:

```bash
export TMPDIR=/lscratch/$SLURM_JOBID; \
mkdir -p $TMPDIR/out; \
mkdir -p $TMPDIR/wrk; \
cp -r /data/NCR_SBRB/freesurfer5.3_subjects/${m} $TMPDIR/out/freesurfer; \
cp /usr/local/apps/freesurfer/6.0.0/license.txt ./ \
fmriprep /data/NCR_SBRB/NCR_BIDS/ $TMPDIR/out \
    participant --participant_label sub-${m} -w $TMPDIR/wrk --use-aroma \
    --nthreads $SLURM_CPUS_PER_TASK --mem_mb 100000 --notrack \
    --fs-license-file ./license.txt; \
mv $TMPDIR/out /data/NCR_SBRB/fmriprep_output/sub-${m}
```

And we can swarm something like the bit above for a few IDs, just for
benchmarking.

```bash
cd /data/NCR_SBRB/
rm fmriprep.swarm;
for m in `cat ~/tmp/testids.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'cp /usr/local/apps/freesurfer/6.0.0/license.txt ./; ' \
        'fmriprep /data/NCR_SBRB/NCR_BIDS/ $TMPDIR/out ' \
        "participant --participant_label sub-${m}" '-w $TMPDIR/wrk --use-aroma ' \
        '--nthreads $SLURM_CPUS_PER_TASK --mem_mb 10000 --notrack ' \
        '--fs-license-file ./license.txt --fs-no-reconall; ' \
        'mv $TMPDIR/out ' "/data/NCR_SBRB/fmriprep_output/sub-${m};" >> fmriprep.swarm;
done
swarm -f fmriprep.swarm --gres=lscratch:10 -g 10 -t 16 --module fmriprep \
     --time=4:00:00 --merge-output --logdir=trash_fmriprep \
     --job-name fmriprep1 --partition quick
```

But make sure the BIDS data re where they're expected to be! For this I'll
ignore the Freesurfer output, as it adds a lot to the computations. It seems
that it's only used in refining the ANTS brainmask in fmriprep preprocessing,
and I don't think this will make or break our analysis. If anything, we'd be
more dependent on good T1s... not sure if we wan to do that.

Using 16 cores, it only went up to 6 Gb, and took about 3h. Also, it only took
about 4G in the local scratch, including working and output directories.

I actually started with 3h, but it was cutting close for a few IDs, so I went up
to 4h, which is the limit for the quick partition.

Also, make sure all scans being processed have the same TR and slice timing!

# 2019-05-31 12:03:16

I created shaw/AROMA/qc_20190531.xlsx to share QC metadata with Luke. But now we
need to test xpc:

```bash
module load xcpengine

export TMPDIR=/lscratch/$SLURM_JOBID; \
mkdir -p $TMPDIR/out; \
mkdir -p $TMPDIR/wrk; \
cp /data/NCR_SBRB/fc-aroma.dsn /data/NCR_SBRB/cohort.csv $TMPDIR/
xcpEngine \
  -c $TMPDIR/cohort.csv \
  -d $TMPDIR/fc-aroma.dsn \
  -i $TMPDIR/work \
  -o $TMPDIR/out \
  -r /data/NCR_SBRB/fmriprep_output/
mv $TMPDIR/out/sub-0977 /data/NCR_SBRB/xcpengine_output/
```

That took less than 30min, using 16 CPUs and less than 5Gb of memory. Local disk
wasn't much different. So, if we were to swarm it:

```bash
rm xcpengine.swarm;
for m in `cat ~/tmp/testids.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'cp /data/NCR_SBRB/fc-aroma.dsn $TMPDIR/; ' \
        'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
        'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
        'xcpEngine -c $TMPDIR/'${m}'.csv ' \
        '-d $TMPDIR/fc-aroma.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
        '-r /data/NCR_SBRB/fmriprep_output/;' \
        'mv $TMPDIR/out/sub-'${m}' /data/NCR_SBRB/xcpengine_output/;' >> xcpengine.swarm;
done
swarm -f xcpengine.swarm --gres=lscratch:10 -g 10 -t 16 --module xcpengine \
     --time=15:00 --merge-output --logdir=trash_xcpengine \
     --job-name xcp1 --partition quick
```

# 2019-06-11 14:52:21

Based on Antonio's e-mail, I should use /usr/local/apps/freesurfer/license.txt
for the license, and that should work...

```bash
cd /data/NCR_SBRB/
rm fmriprep.swarm;
for m in `cat ~/tmp/testids.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'fmriprep /data/NCR_SBRB/NCR_BIDS/ $TMPDIR/out ' \
        "participant --participant_label sub-${m}" '-w $TMPDIR/wrk --use-aroma ' \
        '--nthreads $SLURM_CPUS_PER_TASK --mem_mb 10000 --notrack ' \
        '--fs-license-file /usr/local/apps/freesurfer/license.txt --fs-no-reconall; ' \
        'mv $TMPDIR/out ' "/data/NCR_SBRB/fmriprep_output/sub-${m};" >> fmriprep.swarm;
done
swarm -f fmriprep.swarm --gres=lscratch:10 -g 10 -t 16 --module fmriprep \
     --time=4:00:00 --merge-output --logdir=trash_fmriprep \
     --job-name fmriprep1 --partition quick
```

# 2019-06-18 16:03:42

Let's check for proper ending of fmriprep:

```bash
# caterpie
rm ~/tmp/xcp
for m in `cat /mnt/shaw/AROMA_ICA/kids_n1210_20190618.txt`; do
    fname=/mnt/shaw/AROMA_ICA/fMRIprep_output/sub-${m}/fmriprep/sub-${m}/func/sub-${m}_task-rest_run-1_desc-confounds_regressors.tsv;
    if [ ! -e $fname ]; then
        echo $m;
    else
        echo $m >> ~/tmp/xcp
    fi
done
```

Then, we can check who finished xcp:

```bash
for m in `cat ~/tmp/xcp`; do
    fname=/mnt/shaw/AROMA_ICA/xcpengine_output_AROMA/sub-${m}/fcon/power264/sub-${m}_power264.net
    if [ ! -e $fname ]; then
        echo $m;
    fi
done
```

# 2019-06-19 18:14:27

And we can also run AROMA with GSR:

```bash
rm xcpengine.swarm;
for m in `cat ~/tmp/kids_n1210_20190618.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'cp /data/NCR_SBRB/fc-aroma-gsr.dsn $TMPDIR/; ' \
        'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
        'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
        'xcpEngine -c $TMPDIR/'${m}'.csv ' \
        '-d $TMPDIR/fc-aroma-gsr.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
        '-r /data/NCR_SBRB/fmriprep_output/;' \
        'mv $TMPDIR/out/sub-'${m}' /data/NCR_SBRB/xcpengine_output/;' >> xcpengine.swarm;
done
swarm -f xcpengine.swarm --gres=lscratch:10 -g 10 -t 16 --module xcpengine \
     --time=15:00 --merge-output --logdir=trash_xcpengine \
     --job-name xcpgsr --partition quick
```

# 2019-06-20 10:46:36

Running some scrubbing:

```bash
rm xcpengine.swarm;
for m in `cat ~/tmp/kids_n1210_20190618.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'cp /data/NCR_SBRB/fc-aroma-gsr-p5.dsn $TMPDIR/; ' \
        'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
        'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
        'xcpEngine -c $TMPDIR/'${m}'.csv ' \
        '-d $TMPDIR/fc-aroma-gsr-p5.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
        '-r /data/NCR_SBRB/fmriprep_output/;' \
        'mv $TMPDIR/out/sub-'${m}' /data/NCR_SBRB/xcpengine_output/;' >> xcpengine.swarm;
done
swarm -f xcpengine.swarm --gres=lscratch:10 -g 10 -t 16 --module xcpengine \
     --time=15:00 --merge-output --logdir=trash_xcpengine \
     --job-name xcpgsr-.5 --partition quick
```


# 2019-06-26 16:45:23

Just so we can run it locally:

```bash
TMPDIR=/Users/sudregp/data/;
mkdir -p $TMPDIR/out $TMPDIR/wrk;
for m in `cat ~/data/heritability_change/kids_n1210_20190618.txt`; do
    echo id0,img > ${TMPDIR}/${m}.csv;
    echo sub-${m},sub-${m}/fmriprep/sub-${m}/func/sub-${m}_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/${m}.csv;
    xcpengine-docker -m s -c $TMPDIR/${m}.csv -d $TMPDIR/fc-aroma-p25.dsn \
        -i $TMPDIR/work -o $TMPDIR/out \
        -r /Volumes/Labs/AROMA_ICA/fMRIprep_output/;
done
```

# 2019-06-27 10:38:28

When running XCP single subjects in BW, I do this:

```bash
m=2306
echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'cp /data/NCR_SBRB/fc-aroma-gsr-p5.dsn $TMPDIR/; ' \
        'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
        'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
        'xcpEngine -c $TMPDIR/'${m}'.csv ' \
        '-d $TMPDIR/fc-aroma-gsr-p5.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
        '-r /data/NCR_SBRB/fmriprep_output/;';
```

But if we want to make the scrubbing more automatic: (note the use of single and
double quotes, because $TMPDIR has to be for each swarm, and not the one
currently defined in the environment!)

```bash
rm xcpengine.swarm;
pipe=-gsr-p5-nc;
outdir=/data/NCR_SBRB/xcpengine_output_AROMA${pipe}/
mkdir $outdir;
for m in `cat ~/tmp/kids_n1210_20190618.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
        'cp /data/NCR_SBRB/fc-aroma'${pipe}'.dsn $TMPDIR/; ' \
        'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
        'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
        'xcpEngine -c $TMPDIR/'${m}'.csv ' \
        '-d $TMPDIR/fc-aroma'${pipe}'.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
        '-r /data/NCR_SBRB/fmriprep_output/;' \
        'mv $TMPDIR/out/sub-'"${m} $outdir;">> xcpengine.swarm;
done`
swarm -f xcpengine.swarm --gres=lscratch:10 -g 10 -t 16 --module xcpengine/1.0rc1 \
     --time=20:00 --merge-output --logdir=trash_xcpengine \
     --job-name xcp${pipe} --partition quick
```

# 2019-07-01 14:15:37

Sometimes, for one reason of another, the pipeline doesn't finish. So, let's
make sure we run everyone:

```bash
grep TRUE ~/data/heritability_change/resting_demo_07012019.csv | awk '{FS=","; if ( $9 < 18 ) { print $1 }}' > ~/tmp/kids.txt
grep redo ~/data/heritability_change/resting_demo_07012019.csv | awk '{FS=","; if ( $9 < 18 ) { print $1 }}' >> ~/tmp/kids.txt
for m in `cat ~/tmp/kids.txt`; do
    if [ ! -e /Volumes/Labs/AROMA_ICA/fMRIprep_output/sub-${m}/fmriprep/sub-${m}/func/sub-${m}_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz ]; then
        echo $m;
    fi;
done
```

# 2019-07-03 10:37:43

Now I have to run several pipelines for a few subjects. So, we can do:

```bash
rm xcpengine.swarm;
for g in '' '-gsr'; do
    for p in '' '-p5' '-p25' '-p5-nc' '-p25-nc'; do
        pipe=${g}${p};
        outdir=/data/NCR_SBRB/xcpengine_output_AROMA${pipe}/;
        for m in 0807 2478 2504; do
            echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
                'mkdir -p $TMPDIR/out $TMPDIR/wrk; ' \
                'cp /data/NCR_SBRB/fc-aroma'${pipe}'.dsn $TMPDIR/; ' \
                'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
                'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
                'xcpEngine -c $TMPDIR/'${m}'.csv ' \
                '-d $TMPDIR/fc-aroma'${pipe}'.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
                '-r /data/NCR_SBRB/fmriprep_output/;' \
                'mv $TMPDIR/out/sub-'"${m} $outdir;">> xcpengine.swarm;
        done
    done
done
swarm -f xcpengine.swarm --gres=lscratch:10 -g 10 -t 16 --module xcpengine/1.0rc1 \
     --time=20:00 --merge-output --logdir=trash_xcpengine \
     --job-name xcpredo --partition quick
```
