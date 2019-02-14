# 2019-02-14 09:12:45

Let's keep track of how we're running the FDT pipeline on our own data. First,
let's do it only for the people that have genotype data, as they'l likely be
used for gene expression project.

First, copy everything to the cluster so we can batch convert the data. We can
create a Globus script for that, getting the DICOM folders straight from the
README file to avoid any interactions with the previous DTI preprocessing.

First, start Globus on caterpie:

```bash
/usr/local/neuro/globusconnectpersonal-2.3.3/globusconnectpersonal -start \
    -restrict-paths r/mnt/shaw/MR_data_by_maskid
```

I need to run this in the same computer where globus connect personal is
running, because the same network path used to get directory names goes in the
Globus input file!

```bash
maskid_file=~/tmp/prs_dti.txt
out_dir=/scratch/sudregp/
net_dir=/mnt/shaw/

hpc=e2620047-6d04-11e5-ba46-22000b92c6ec
caterpie=74b32b98-a3aa-11e7-adbf-22000a92523b
myfile=~/tmp/dti_dcm_transfers.dat

ssh -qt helix.nih.gov "if [ ! -d ${out_dir}/dcm_dti ]; then mkdir ${out_dir}/dcm_dti/; fi";

rm -rf $myfile
for m in `cat $maskid_file`; do
    echo Copying $m
    ssh -qt helix.nih.gov "if [ ! -d ${out_dir}/dcm_dti/${m} ]; then mkdir ${out_dir}/dcm_dti/${m}; fi";
    scp -qt ${net_dir}/MR_data_by_maskid/${m}/E*/cdi* helix:${out_dir}/dcm_dti/${m}/;

    # find name of date folders
    ls -1 $net_dir/MR_data_by_maskid/${m}/ | grep -e ^20 > ~/tmp/date_dirs;

    # for each date folder, check for dti scans
    cnt=1
    while read d; do
        grep dti $net_dir/MR_data_by_maskid/${m}/${d}/*README* > ~/tmp/dti;
        awk '{for(i=1;i<=NF;i++) {if ($i ~ /Series/) print $i}}' ~/tmp/dti | sed "s/Series://g" > ~/tmp/dti_clean
        while read line; do
            mr_dir=`echo $line | sed "s/,//g"`;
            echo "--recursive ${net_dir}/MR_data_by_maskid/${m}/${d}/${mr_dir}/ ${out_dir}/dcm_dti/${m}/${mr_dir}/" >> $myfile;
            let cnt=$cnt+1;
        done < ~/tmp/dti_clean;
    done < ~/tmp/date_dirs;
done;

# assuming globus cli is installed as in:
# pip2.7 install --upgrade --user globus-cli
~/.local/bin/globus transfer $caterpie $hpc --batch --label "dti copy" < $myfile
```

Next thing is to start converting them.

```bash
module load TORTOISE
# quick hack that if there are two cdi files we should use the last one
gradient_file=`/bin/ls -1 | grep cdi | tail -1`;
echo Using $gradient_file for gradients!
/bin/ls -1 | grep -v cdi > mr_dirs.txt;
# remove first line
tail -n +2 $gradient_file | split -l 20 -a 1 -d - grads;
cnt=0
nii='-innii'
bval='-inbval'
bvec='-inbvec'
for mr_dir in `cat mr_dirs.txt`; do
    ImportDICOM -i $mr_dir -o s${cnt} -b 1100 -g grads${cnt};
    TORTOISEBmatrixToFSLBVecs s${cnt}_proc/s${cnt}.bmtxt;
    nii=$nii' 's${cnt}_proc/s${cnt}.nii
    bval=$bval' 's${cnt}_proc/s${cnt}.bvals
    bvec=$bvec' 's${cnt}_proc/s${cnt}.bvecs
    let cnt=${cnt}+1;
done
fat_proc_convert_dcm_dwis $nii $bval $bvec -prefix dwi_comb -flip_x
```

That matched exactly what I got from opening the .list file in DIFFCALC and
exporting the raw DWI (edti.list) to FSL UNSORTED format.

The idea behind doing it per session comes from an email from Irfan, who said
that TORTOISE v2 used directory structure, but ImportDICOM (and dcm2nixx) use
the DICOM series number written in the header. Somehow the scanner is writing
the series number wrong, and then the DICOMs get imported incorrectly. If we do
it by scan, like TORTOISE v2 used to do, and provide the correct gradients, it
should work fine.

After this, it's just the exact same thing as the PNC pipeline, with probably
some changed in the bedpostx to only fit one orientation. In fact, I could do
that for the PNC as well, just so we don't have two different options.
