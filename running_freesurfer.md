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
scp ~/tmp/missing.txt biowulf.nih.gov:/scratch/sudregp/mprage/
scp -r ~/tmp/missing/* biowulf.nih.gov:/scratch/sudregp/mprage/
```

