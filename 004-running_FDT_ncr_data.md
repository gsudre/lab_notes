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
    scp -q ${net_dir}/MR_data_by_maskid/${m}/E*/cdi* helix:${out_dir}/dcm_dti/${m}/;

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
for m in `cat ~/tmp/prs_dti.txt`; do
    echo "bash ~/research_code/dti/convert_ncr_to_nii.sh /scratch/sudregp/dcm_dti/$m" >> swarm.ncr
done
swarm -t 2 --job-name ncr2nii --time 15:00 -f swarm.ncr --partition quick \
    --logdir trash_ncr2nii -m TORTOISE,afni
```

That matched exactly what I got from opening the .list file in DIFFCALC and
exporting the raw DWI (edti.list) to FSL UNSORTED format.

After this, it's just the exact same thing as the PNC pipeline, with probably
some changed in the bedpostx to only fit one orientation. In fact, I could do
that for the PNC as well, just so we don't have two different options.

# 2019-02-17 18:55:03

For the next step, we need to check who converted correctly. Maybe just check
the number of vectors in the gradient file against the number of converted
volumes in the NIFTI?

```bash
res_file=~/tmp/ncr_conversion.csv
echo mask.id,expected,converted,nii > $res_file;
cd /scratch/sudregp/dcm_dti
for m in `cat ~/tmp/prs_dti.txt`; do
    nvol=`cat ${m}/dwi_comb_cvec.dat | wc -l`;
    gradient_file=`/bin/ls -1 ${m}/ | /bin/grep cdi | tail -1`;
    nexp=`cat ${m}/${gradient_file} | wc -l`;
    let nexp=$nexp-1;  # removing first line
    if [ -e ${m}/dwi_comb.nii.gz ]; then
        nii='y';
    else
        nii='n';
    fi;
    echo $m,$nexp,$nvol,$nii >> $res_file;
done
```

Just noticed it'll be a bitch to get the sl_spec file for our data... we'll need
to run dcm2niix individually, then add to the sequence to offset the acquisiton
of 20 volumes that were later combined. Finally, we need to find a way to
combine in the 99 sequence... I sent an e-mail to Joelle to see if she can help.

# 2019-02-20 11:11:33

Joelle said that the temporal order for the slices in the diffusion protocol
will be interleaved, being a first pass through the volume to acquire the odd
slices then the next pass for the even slices.  This seems to be what is
reflected in the my-slspec.txt file.

So, I'll just go ahead and use the same slspec file that I created for the PNC
data, adjusting for the total number of slices. But that's the trick, as the
number of slices varies for each scan!!! OK, so we'll have to construct the
slices file on the fly.

## Copying data to NCR_SBRB

The problem of having this in scratch is that I'm the only one who can process
it. Also, it expires. So, let's copy the converted files to the shared space.
But we only copy the successful stuff, and the converted files. We can keep on
dealing with the rest of the DCMs later.

```bash
mkdir /data/NCR_SBRB/dti_fdt
for m in `cat ~/tmp/mylist.txt`; do
    mkdir /data/NCR_SBRB/dti_fdt/$m;
    echo $m;
    cd /scratch/sudregp/dcm_dti/${m};
    cp -r dwi_comb* QC /data/NCR_SBRB/dti_fdt/$m/;
    chgrp -R NCR_SBRB /data/NCR_SBRB/dti_fdt/$m;
    chmod -R 770 /data/NCR_SBRB/dti_fdt/$m;
done
```

Then, it's easy to do:

```bash
cd /data/NCR_SBRB/dti_fdt
for m in `cat list1`; do
    echo "bash ~/research_code/dti/fdt_ncr_eddy.sh /data/NCR_SBRB/dti_fdt/${m}" >> swarm.fdt;
done;
split -l 310 swarm.fdt
swarm -g 4 --job-name fdt --time 4:00:00 -f xaa --partition gpu \
    --logdir trash_fdt --gres=gpu:k80:2
```

I'm having to redo some of the conversions manually. Here's how I'm approaching
that:

```bash
cd ../0545; rm -rf QC __* dwi_comb* *proc mr_dirs.txt grad*; ls
bash ~/research_code/dti/convert_ncr_to_nii.sh /scratch/sudregp/dcm_dti/0545
fslinfo dwi_comb | grep -e "^dim4" && wc -l dwi_comb_cvec.dat && wc -l cdiflist0*
```

# 2019-02-22 10:04:30

Let's start running the first part of autoPtx for the IDs that finished eddy.
Note that here I'll need to change the assumptions because we have two different
bvals, as described in https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FDT/UserGuide.

I'm following the same steps I did for PNC, except that in that case I didn't
run model=2!

```bash
# run in helix so we don't overload BW filesystem
for m in `cat /data/NCR_SBRB/dti_fdt/autoptx1`; do
    if [ ! -e /data/NCR_SBRB/dti_fdt/${m}/eddy_s2v_unwarped_images.nii.gz ]; then
        echo $m;
    fi;
done
```

```bash
# run in helix so we don't overload BW filesystem
for m in `cat /data/NCR_SBRB/dti_fdt/done_eddy2.txt`; do
    echo Copying $m;
    cd /data/NCR_SBRB/dti_fdt/${m};
    cp eddy_s2v_unwarped_images.nii.gz data.nii.gz;
    cp dwi_comb_bval.dat bvals;
    cp eddy_s2v_unwarped_images.eddy_rotated_bvecs bvecs;
    cp b0_brain_mask.nii.gz nodif_brain_mask.nii.gz;
    chgrp NCR_SBRB data.nii.gz bvals bvecs nodif_brain_mask.nii.gz;
    chmod 770 data.nii.gz bvals bvecs nodif_brain_mask.nii.gz;
done
```

Then it's time to run autoPtx.

We can split it by users because it adds everything to the
same directory, and just increments the final subject list.

Run it a long interactive session, because even though it schedules bedpostx, it
still runs all kinds of registrations through FSL, so biowulf headnode won't cut
it! Also, it takes a while to load new scans because it does the processing in
between. So, it *might* be alright to do 400 or so scans at the same time. While
the scan is processing in the interactive node, it gives time for the ones
queued up to be evaluated. This way the queue doesn't get clogged up.

```bash
data='';
for m in `cat xab`; do
    data=$data' '${m}/data.nii.gz;
done
/data/NCR_SBRB/software/autoPtx/autoPtx_1_preproc $data;
```

xaa: jen
xab: philip

And of course we still need part 2 when we're done.