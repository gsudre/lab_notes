# 2019-07-25 11:24:52

Il'll take a short detour here and see if the 36P designs do a bit better as far
as removing movement. Because the cluster is down, I'll keep it to the regular
pipeline first, but we can always do the other ones later if it looks promising.

```r
scans_file = '/Volumes/Labs/AROMA_ICA/filtered_minFD_2scans.csv'
scans = read.csv(scans_file)
subjs = unique(as.character(scans$subj))
scans_file = '/Volumes/Labs/AROMA_ICA/filtered_minFD_3scans.csv'
scans = read.csv(scans_file)
subjs = unique(c(subjs, unique(as.character(scans$subj))))
write.table(subjs, file='~/data/36P/maskids_23.txt', row.names=F, col.names=F, quote=F)
```

```bash
# desktop
cd ~/data/36P
TMPDIR=`pwd`;
mkdir -p $TMPDIR/out $TMPDIR/work;
for s in `cat maskids_23.txt`; do
    echo id0,img > ${TMPDIR}/${s}.csv;
    echo ${s},${s}/fmriprep/${s}/func/${s}_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/${s}.csv;
    xcpengine-docker -c $TMPDIR/${s}.csv -d /xcpEngine/designs/fc-36p.dsn \
        -i $TMPDIR/work -o $TMPDIR/out \
        -r /Volumes/Labs/AROMA_ICA/fMRIprep_output/;
done
```

# 2019-07-26 11:25:59

The worry then is 2-fold:

1) remove crappy scans; then
2) make sure changes in movement are not driving the results in change of connectivity

So, we start with the best 2 or 3 scans per person. Here, we measured best based
on mean FD, but it could even be a ranked combination of FD and RMS. Then, it's
a matter of removing any scans.

# 2019-07-29 16:10:29

My initial tests showed the vanilla 36P pipeline to remove the bivariate
distribution we were seeing across the board with AROMA. Let's see what we get
when we run it in the cluster, now that it's back in business:

```bash
pipe='fc-36p';
cd /data/NCR_SBRB/
rm xcpengine.swarm;
outdir=/data/NCR_SBRB/xcpengine_output_${pipe}/
mkdir $outdir;
for m in `cat /data/NCR_SBRB/maskids_23.txt`; do
    echo 'export TMPDIR=/lscratch/$SLURM_JOBID; ' \
        'mkdir -p $TMPDIR/out $TMPDIR/work; ' \
        'echo id0,img > ${TMPDIR}/'${m}'.csv; ' \
        'cp /data/NCR_SBRB/'${pipe}'.dsn $TMPDIR/; ' \
        'echo sub-'${m}',sub-'${m}'/fmriprep/sub-'${m}'/func/sub-'${m}'_task-rest_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz >> ${TMPDIR}/'${m}'.csv; ' \
        'xcpEngine -c $TMPDIR/'${m}'.csv ' \
        '-d $TMPDIR/'${pipe}'.dsn -i $TMPDIR/work -o $TMPDIR/out ' \
        '-r /data/NCR_SBRB/fmriprep_output/;' \
        'mv $TMPDIR/out/sub-'"${m} $outdir;">> xcpengine.swarm;
done
swarm -f xcpengine.swarm --gres=lscratch:10 -g 10 -t 16 --module xcpengine/1.0rc1 \
     --time=30:00 --merge-output --logdir=trash_xcpengine_${pipe} \
     --job-name xcp${pipe} --partition quick,norm
```
