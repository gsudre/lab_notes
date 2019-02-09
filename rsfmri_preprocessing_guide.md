# Resting-state FMRI (rsFMRI) pre-processing guide
* Originally written by Ruiz, Esteban (NIH/NHGRI) [F], September 12, 2018
* Updated by Gustavo Sudre, December 2018

Most of these commands will be run in caterpie, so ensure a common environment
for all users. For them to work, make sure you have SSH keys set up to connect
to Biowulf/helix, so you don't have to type your password every time. More instructions here: https://hpc.nih.gov/docs/sshkeys.html

## Copying files to Biowulf
All files need to be copied from the shared drive to Biowulf for faster
processing. This can all be accomplished with one script. Copy it from
https://github.com/gsudre/research_code/blob/master/fmri/copy_rest_files.sh
to wherever you want in caterpie. Then, for a list of mask ids:

```bash
for m in `cat ~/tmp/maskids.txt`; do
    bash ~/research_code/fmri/copy_rest_files.sh ${m} /data/NCR_SBRB/tmp/ /mnt/shaw/;
done
```

A faster option would be to use Globus. Make sure you have your NIH account
properly setup at http://www.globus.org, and that you have it installed for your
user account in caterpie (come find me if you have problems). In a nutshell,
you'll need (using Anaconda's path):

```bash
pip2.7 install --upgrade --user gcc
pip2.7 install --upgrade --user globus-cli
```

You only need to do that once. Then, still in caterpie:

```bash
bash ~/research_code/fmri/copy_rest_files_globus.sh maskids.txt /data/NCR_SBRB/tmp/ /mnt/shaw/
```

That should start a Globus transfer called "rsfmri copy", which you can monitor
in your web browser.

## Converting from DICOMs to AFNI format
Next we convert all DICOMs we just copied to the AFNI format we'll need for
pre-processing. This takes a while, so there are two options here:

1) Slowly use the conversion script as the files are being copied above.
2) Wait until everything is copied and swarm it in the cluster.

I'm just using the swarm version here, but it's straight forward to pluck out
the command and do it manually (or in a for-loop) if needed. If anything, you
can construct the swarm file, and then use each line individually in an
interactive node.

```bash
rm swarm.convert
for m in `cat ~/tmp/maskids.txt`; do
   echo "bash ~/research_code/fmri/convert_sorted_rest_files.sh ${m} /data/NCR_SBRB/tmp/dcm_mprage /data/NCR_SBRB/tmp/dcm_rsfmri/ /scratch/${USER}/rsfmri/" >> swarm.convert
done
swarm -f swarm.convert -g 4 -t 2 --time 30:00 --partition quick --logdir trash_convert --job-name afni_convert -m afni
```

And after this much work, let's go ahead and copy back to the shared drive all
converted files, in case we need to use them again in the future. There aren't
many files, so using Globus here might be more work than it needs to be:

```bash
# caterpie
for m in `cat ~/tmp/maskids.txt`; do
    echo "Copying back converted data for ${m}";
    mkdir /mnt/shaw/MR_data_by_maskid/${m}/afni;
    scp -q helix.nih.gov:/scratch/${USER}/rsfmri/${m}/* /mnt/shaw/MR_data_by_maskid/${m}/afni/;
done
```

When everything is converted, it's time to update the QC spreadsheet with the
number of resting files converted.

```bash
echo maskid,converted > ~/tmp/converted.csv
for m in `cat ~/tmp/maskids.txt`; do
   nconvert=`/bin/ls -1 /mnt/shaw/MR_data_by_maskid/${m}/afni/rest*HEAD | wc -l`;
   echo ${m},${nconvert} >> ~/tmp/converted.csv;
done
```

And then it's time to run the actual processing script. First, we need to copy
the Freesurfer files as well, as they're needed for the pre-processing. Globus
is again the best option here, as Freesurfer outputs have lots of files:

```bash
# in caterpie
hpc=e2620047-6d04-11e5-ba46-22000b92c6ec
caterpie=74b32b98-a3aa-11e7-adbf-22000a92523b
myfile=~/tmp/freesurfer_transfers.dat
rm -rf $myfile
for m in `cat ~/tmp/maskids.txt`; do 
    echo "--recursive /mnt/shaw/freesurfer5.3_subjects/${m}/ /data/NCR_SBRB/freesurfer5.3_subjects/${m}/" >> $myfile;
done;
~/.local/bin/globus transfer $caterpie $hpc --batch --label "freesurfer copy" < $myfile
```

Finally, when everything is properly copied, we run afni_proc:

```bash
#biowulf
rm swarm.ss
while read m; do 
    echo "source /usr/local/apps/freesurfer/5.3.0/SetUpFreeSurfer.sh; bash ~/research_code/fmri/run_resting_afni_proc_subjectSpace.sh ${m}" >> swarm.ss;
done < ~/tmp/maskids.txt
swarm -f swarm.ss -g 4 --time 4:00:00 --partition quick --logdir trash_ss \
    --job-name afni_proc -m afni,freesurfer/5.3.0
```

The last stage is to figure out who finished properly and needs to be copied
back. To do that, we need to check for the final files, and also for
alignment.

```bash
net_dir=/scratch/${USER}/rsfmri
while read m; do
    if [ -e ${net_dir}/${m}/${m}.rest.subjectSpace.results/errts.${m}.fanaticor+orig.HEAD ]; then
        echo $m >> ~/tmp/need_qc.txt;
    fi;
done < ~/tmp/maskids.txt
```

We could also run the negative of the command above to figure out which mask ids
didn't run and check the output file inside the mask id folder for the reasons.
Most of the time:

* there are rest runs without enough TRs to even get started (e.g. started and
  stoped the sequence fast); OR
* after censoring outlier and motion TRs, there weren't enough TRs left for the
  deconvolution procedure
  
Check the log and update the tracker spreadsheet accordingly. You might be able
to recover from the first problem above by just renaming the bad rest run and
re-running the preprocessing. For example:

```bash
m=2441;
cd /scratch/${USER}/rsfmri/${m}/
rm -rf *subjectSpace*
3drename rest1+orig ignore_rest1
bash ~/research_code/fmri/run_resting_afni_proc_subjectSpace.sh ${m}
```

But not having enough TRs for the deconvolution can only be fixed by relaxing
our censoring limits...

Then, we compute some of the metrics that will determine 
Then we prepare the alignment QC images:

```bash
net_dir=/scratch/${USER}/rsfmri
mkdir ${net_dir}/resting_alignment/
Xvfb :88 -screen 0 1024x768x24 &
while read s; do
  cd ${net_dir}/${s}/${s}.rest.subjectSpace.results;
  \@snapshot_volreg ${s}_SurfVol_al_junk+orig vr_base_min_outlier+orig align 88;
  mv -v align.jpg ${net_dir}/resting_alignment/${s}.jpg;
done < ~/tmp/need_qc.txt
```

To do the visual QC in biowulf, use eog. Look here for guidance on how to
interpret alignment:
https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/tutorials/auto_image/auto_@snapshot_volreg.html#tut-auto-snapshot-volreg


```bash
cd /scratch/${USER}/rsfmri/resting_alignment/
eog ????.jpg
```

If there are mask ids that didn't align properly, we need to check which is the
best alignment to go with next. But in fact, we should always check if there is
a better alignment than the one we used. So, we can swarm the script to make the
other options:

```bash
#biowulf
rm swarm.uber
while read m; do 
    echo "bash ~/research_code/fmri/make_alignment_options.sh ${m} /scratch/sudregp/rsfmri" >> swarm.uber;
done < ~/tmp/maskids.txt
swarm -f swarm.uber -g 4 --time 4:00:00 --partition quick --logdir trash_uber \
    --job-name afni_uber -m afni
```

# 2019-02-08 19:14:57

Use something like this in caterpie to copy back the data:

```bash
# caterpie
for m in `cat ~/tmp/maskids.txt`; do
    echo "Copying back converted and processed data for ${m}";
    mkdir /mnt/shaw/MR_data_by_maskid/${m}/afni;
    scp -qr helix.nih.gov:/scratch/${USER}/rsfmri/${m}/* /mnt/shaw/MR_data_by_maskid/${m}/afni/;
    scp -qr helix.nih.gov:/data/NCR_SBRB/freesurfer5.3_subjects/${m}/SUMA /mnt/shaw/freesurfer5.3_subjects/${m}/;
done
```

REMEMBER TO COPY BACK ALIGNMENT QC IMAGES AS WELL!!!

TODO

* check how many scans all mask ids have and their status of true and false
* visual QC on all last scans
* re-align if necessary
* re-run all QC metrics for all IDs
* check who needs to be called for rsFMRI test set