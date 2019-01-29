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
cd ~/freesurfer_logs
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

# 2018-12-18 10:28:32

## Generating QC pictures

We use the ENIGMA scripts to generate the Freesurfer QC pictures. There are
internal and external pictures. For external, it's simple:

```bash
for i in {2404..2435}; do echo $i >> ~/tmp/fs_ids.txt; done
cd /Volumes/Shaw/segmentation_qc/enigma_external_freesurfer5.3/3t/
bash ~/enigma/fsqc_noLowRes.sh ~/tmp/fs_ids.txt
```

For internal, we need Matlab. Then, we run something like this:

```matlab
addpath(genpath('~/enigma/QC_ENIGMA/'))
a = dlmread('~/tmp/fs_ids.txt');
FS_directory = '/Volumes/Shaw/freesurfer5.3_subjects/';
QC_output_directory='/Volumes/Shaw/segmentation_qc/enigma_internal_freesurfer5.3/3t/';
for x = 1:size(a,1)
b=sprintf('%04d', a(x));
func_make_corticalpngs_ENIGMA_QC_erasmus(QC_output_directory, b, [FS_directory,'/',b,'/mri/orig.mgz'], [FS_directory,'/',b,'/mri/aparc+aseg.mgz']);
end
```

# 2019-01-28 13:08:43

I needed to extract the Euler number for each subject. Freesurfer extracts one
per hemisphere, and more information can be seen here: 

https://www.sciencedirect.com/science/article/pii/S1053811917310832?via%3Dihub

https://vitallongevity.utdallas.edu/cnl/wp-content/uploads/2015/10/FREESURFER_GUIDE_CNL_2015Oct.pdf

The idea is that more holes in the original surface gives a more negative
number. So, we average the numbers in both hemispheres, and the higher the
number, the better the original result (before any corrections):

```bash
res_file=~/tmp/euler_numbers.csv
echo subjid,euler_LH,euler_RH,mean_euler > $res_file;
for s in `cat ~/tmp/have_imaging.txt`; do
    euler_lh=`grep -A 1 "Computing euler" /Volumes/Shaw/freesurfer_pnc/${s}/scripts/recon-all.log | tail -1 | awk '{ print $4 }' | sed "s/,//"`;
    euler_rh=`grep -A 1 "Computing euler" /Volumes/Shaw/freesurfer_pnc/${s}/scripts/recon-all.log | tail -1 | awk '{ print $7 }' | sed "s/,//"`;
    mean_euler=`echo "( $euler_lh + $euler_rh ) / 2" | bc`;
    echo $s,$euler_lh,$euler_rh,$mean_euler >> $res_file;
done
```

And then we can also run MRIQC. That's a bit tricky, as it needs a BIDS root. We
have not converted our data ot the BIDS format yet, so we need to fake it. One
of the requirements is a .json for each .nii, so I'll need to recreate
everything that was already done for Freesurfer...

```bash
module load jo
module load afni
mkdir /scratch/sudregp/BIDS
cd /scratch/sudregp/BIDS
# installed jo through homebrew
jo -p "Name"="Some bogus name" "BIDSVersion"="1.0.2" >> dataset_description.json;
for s in `cat /data/NCR_SBRB/pnc/have_imaging.txt`; do
    mkdir sub-${s};
    m=001;
    mkdir sub-${s}/ses-${m};
    mkdir sub-${s}/ses-${m}/anat;
    # assuming the .tar.gz have been decompressed already
    dcm2niix_afni -o sub-${s}/ses-${m}/anat/ -f sub-${s}_ses-${m}_T1w /data/NCR_SBRB/pnc/${s}/T1_3DAXIAL/Dicoms/;
done
```

Then, we need to run it in the cluster:

```bash
for s in `cat /data/NCR_SBRB/pnc/have_imaging.txt`; do
    echo "bash ~/research_code/lab_mgmt/run_mriqc.sh ${s}" >> swarm.mriqc;
done
swarm -g 8 -f swarm.mriqc --job-name mriqc --time 8:00:00 --logdir trash_mriqc --gres=lscratch:40
```

Actually, I didn't run the above in the cluster. I ended up just doing a massive
call to mriqc, without specifying oatient labs, which calls al possible
subjects. I ran it in lscratch in an interactive session. The problem is
that I didn't specify a maximum number of processes, so it's doubling what it
can use. If no one complains, I'll leave it like that. (nevermind, it broke...
ran out of room).

Otherwise, it looks like mriqc is a module in Biowulf now, so I can simply do
something like:

```bash
module load mriqc
mriqc /scratch/sudregp/BIDS /scratch/sudregp/out.ds001/ participant --participant_label 608720575588 -m T1w -w MRIqc.work.ds001/
```

And do the above in a swarm (i.e. I don't need my wrapper anymore).

The swarm looks like this:

```bash
for s in `cat /data/NCR_SBRB/pnc/have_imaging.txt`; do
    echo "mriqc /scratch/sudregp/BIDS /scratch/sudregp/mriqc_output participant --participant_label ${s} -m T1w -w /scratch/sudregp/mriqc_work" >> swarm.mriqc;
done
swarm -g 8 -f swarm.mriqc --job-name mriqc --time 4:00:00 --logdir trash_mriqc -m mriqc --partition quick --gres=lscratch:40
```

And then we need to collect all results:

```bash
module load singularity
export SINGULARITY_CACHEDIR=/data/sudregp/singularity/
singularity exec -B /scratch/sudregp:/mnt docker://poldracklab/mriqc:latest mriqc /mnt/BIDS /mnt/mriqc_output group --no-sub -w /mnt/mriqc_work -m T1w
```

I saved everything back in the shared drive in PNC_structural.