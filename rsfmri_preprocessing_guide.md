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
   echo "bash ~/research_code/fmri/convert_sorted_rest_files.sh ${m} /data/NCR_SBRB/tmp/dcm_mprage /data/NCR_SBRB/tmp/dcm_rsfmri/ /scratch/sudregp/rsfmri/" >> swarm.convert
done
swarm -f swarm.convert -g 4 -t 2 --time 30:00 --partition quick --logdir trash --job-name afni_convert -m afni
```

When everything is converted, it's time to update the QC spreadsheet with the
number of resting files converted.

```bash
echo maskid,converted > ~/tmp/converted.txt
for m in `cat ~/tmp/maskids.txt`; do
   nconvert=`/bin/ls -1 /scratch/sudregp/rsfmri/${m}/rest*HEAD | wc -l`;
   echo ${m},${nconvert} >> ~/tmp/converted.txt;
done
```

And after this much work, let's go ahead and copy back to the shared drive all
converted files, in case we need to use them again in the future:



# TODO
* Populate QC spreadsheet with number of converted runs
* Copy back to shared drive converted files
* Run preprocessing script