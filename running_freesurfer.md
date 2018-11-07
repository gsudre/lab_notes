# 2018-11-05 11:47:46

As I migrate my notes to GitHub, let's create a fresh one on how to run
Freesurfer on new scans we acquire. This is in a nutshell just the highlights of
the older Evernote note.

We basically single out the best MPRAGE for each scan, and then run Freesurfer
by swarming it in Biowulf:

```bash
rm -rf ~/tmp/missing*
mkdir ~/tmp/missing
for i in {2404..2435}; do echo $i >> ~/tmp/missing.txt; done 
python ~/research_code/lab_mgmt/copy_mprages.py ~/tmp/missing.txt /Volumes/Shaw/ ~/tmp/missing/
cp -rv ~/tmp/missing/* /Volumes/Shaw/best_mprages/
scp ~/tmp/missing.txt helix.nih.gov:/scratch/sudregp/mprage/
scp -r ~/tmp/missing/* helix.nih.gov:/scratch/sudregp/mprage/
```

```bash
# biowulf
cd ~/freesurfer_logs/
# do some cleanup if needed
while read s; do f=`/bin/ls /scratch/sudregp/mprage/${s}/*0001.dcm | tr -d '\n'`; echo "source /usr/local/apps/freesurfer/5.3.0/SetUpFreeSurfer.sh; recon-all -i $f -subjid ${s} -all -openmp 8 | tee -a ${s}_freesurfer.log" >> freesurfer.swarm; done < /scratch/sudregp/mprage/missing.txt

swarm -g 8 -t 8 --job-name freesurfer --time 48:00:00 -f freesurfer.swarm -m freesurfer/5.3.0 --logdir trash
```

Then, to check who finished:

```bash
cd ~/freesurfer_logs/
rm finished.txt
while read s; do 
    if grep -q "finished without error" ~/data/MEG_structural/freesurfer/${s}/scripts/recon-all.log; then
        echo $s >> finished.txt;
    else
        echo $s;
    fi;
done < /scratch/sudregp/mprage/missing.txt
```

If we need to restart a subject:

```bash
while read s; do echo "recon-all -make all -s ${s} | tee -a ${s}_freesurfer.log" >> redo.swarm; done < redo  
```

Then, to bring everyone back:

```bash
# desktop
scp helix:~/freesurfer_logs/finished.txt ~/tmp
while read s; do
    scp -r helix:~/data/MEG_structural/freesurfer/${s} /Volumes/Shaw/freesurfer5.3_subjects/;
done < ~/tmp/finished.txt
```

When all done, run something like this in biowulf for cleanup:

```bash
while read s; do
    echo $s;
    mv ${s}_freesurfer.log old/;
    rm -rf /scratch/sudregp/mprage/${s};
    rm -rf ~/data/MEG_structural/freesurfer/${s};
done < finished.txt
```


