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
for m in `cat ~/tmp/maskids.txt`; do
   echo "bash ~/research_code/fmri/convert_sorted_rest_files.sh ${m} /data/NCR_SBRB/tmp/dcm_mprage /data/NCR_SBRB/tmp/dcm_rsfmri/ /scratch/${USER}/rsfmri/" >> swarm.convert
done
swarm -f swarm.convert -g 4 -t 2 --time 30:00 --partition quick --logdir trash --job-name afni_convert -m afni
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
echo maskid,converted > ~/tmp/converted.txt
for m in `cat ~/tmp/maskids.txt`; do
   nconvert=`/bin/ls -1 /mnt/shaw/MR_data_by_maskid/${m}/afni/rest*HEAD | wc -l`;
   echo ${m},${nconvert} >> ~/tmp/converted.txt;
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
swarm -f swarm.ss -g 4 --time 4:00:00 --partition quick --logdir trash \
    --job-name afni_proc -m afni,freesurfer/5.3.0
```

The last stage is to figure out who finished properly and needs to be copied
back. To do that, we need to check for the final files, and also for
alignment...

REMEMBER TO COPY BACK SUMA FILES AS WELL!!!