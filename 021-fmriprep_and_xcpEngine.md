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
     --time=3:00:00 --merge-output --logdir=trash_fmriprep \
     --job-name fmriprep1 --partition quick
```

But make sure the BIDS data re where they're expected to be! For this I'll
ignore the Freesurfer output, as it adds a lot to the computations. It seems
that it's only used in refining the ANTS brainmask in fmriprep preprocessing,
and I don't think this will make or break our analysis. If anything, we'd be
more dependent on good T1s... not sure if we wan to do that.

Using 16 cores, it only went up to 6 Gb, and took about 3h. Also, it only took
about 4G in the local scratch, including working and output directories.

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
