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