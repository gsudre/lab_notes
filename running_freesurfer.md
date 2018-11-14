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

# 2018-11-14 13:46:02

Here are some lines on how to create a downsamples version of the Freesurfer results, vertex-wide. Here are some llustrations: 

https://brainder.org/2016/05/31/downsampling-decimating-a-brain-surface/

But the default options generate 163842 values per hemisphere. For this, I'll use fsaverage4, which gives me 2562 per hemisphere:

```bash
for meas in volume area thickness; do
    for h in lh rh; do
        mris_preproc --f subjects_list.txt --target fsaverage4 --meas $meas --hemi $h --out ${h}.${meas}.ico4.mgh;
    done;
done
```

Of course, to do the usual stuff without downsampling that we used to do, just use the regular fsaverage. But you might also want to smooth it a bit, like:

```bash
mri_surf2surf --hemi lh --s fsaverage --sval lh.thickness.00.mgh --fwhm 10 --cortex --tval lh.thickness.10.mgh
```

Finally, let's put it into a format we can read. Bring it back form Biowulf, and run this locally:

```python
nsubjs = 260
from rpy2.robjects import r
from rpy2.robjects.numpy2ri import numpy2ri
import surfer
import numpy as np
for meas in ['area', 'volume', 'thickness']:
    for h in ['lh', 'rh']:
        fname = '%s.%s.ico4' % (h, meas)
        data = surfer.io.read_scalar_data(fname + '.mgh')
        nvox = data.shape[0]/nsubjs
        print('%s, nvox = %.2f' % (fname, nvox))
        data = data.reshape([nsubjs, int(nvox)])
        array = np.array(data, dtype="float64")
        ro = numpy2ri(array)
        r.assign('data', ro)
        r("save(data, file='%s.gzip', compress=TRUE)" % fname)
```


