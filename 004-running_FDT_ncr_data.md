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
    cp -r dwi_comb* /data/NCR_SBRB/dti_fdt/$m/;
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

And of course we still need part 2 when we're done.

# 2019-02-28 08:19:32

Let's start running part 2 while we wait on some manual conversions.

And we also need to figure out which IDs had the set_slice error, which might
indicate outliers from eddy:

```bash
for f in `grep -l set_slice trash_fdt/*o`; do
    head -n 2 $f | tail -n -1 | cut -d" " -f 5 | cut -d"/" -f 5 >> ~/tmp/set_slice.txt;
done
```

# 2019-03-01 15:14:28

Now that eddy is ran for all IDs, time to copy over the brain masks. When
autoptx1 is done, we can copy over the rest of the QC.

```bash
mkdir /data/NCR_SBRB/dti_fdt/summary_QC
cd /data/NCR_SBRB/dti_fdt/summary_QC/
mkdir brainmask
mkdir transform
mkdir DEC
mkdir SSE
for m in `cat ~/tmp/myids.txt`; do
    echo ${m}
    cp ../${m}/QC/brain_mask.axi.png brainmask/${m}.axi.png
    cp ../${m}/QC/brain_mask.sag.png brainmask/${m}.sag.png
    cp ../${m}/QC/brain_mask.cor.png brainmask/${m}.cor.png
done
```

Some IDs didn't have brain mask QC ready for some reason...

```bash
for s in `cat ~/tmp/myids.txt`; do
    if [ ! -e /data/NCR_SBRB/dti_fdt/${s}/QC/brain_mask.axi.png ]; then
        cd /data/NCR_SBRB/dti_fdt/${s};
        mkdir QC;
        @chauffeur_afni                             \
            -ulay  dwi_comb.nii.gz[0]                         \
            -olay  b0_brain_mask.nii.gz                        \
            -opacity 4                              \
            -prefix   QC/brain_mask              \
            -montx 6 -monty 6                       \
            -set_xhairs OFF                         \
            -label_mode 1 -label_size 3             \
            -do_clean
    fi;
done
```

Even though autoPtx is not done for all IDs yet, we can run our mean props
estimate just to see if we are getting that weird binormal distribution again...

```bash
mydir=/lscratch/${SLURM_JOBID}/
mean_props=~/tmp/mean_props.csv;
echo "id,mean_fa,mean_ad,mean_rd,nvox" > $mean_props;
for m in `cat ~/tmp/myids.txt`; do
    echo $m;
    cd /data/NCR_SBRB/dti_fdt/preproc/$m &&
    if [ -e dti_FA.nii.gz ]; then
        3dcalc -a dti_FA.nii.gz -expr "step(a-.2)" -prefix ${mydir}/my_mask.nii 2>/dev/null &&
        fa=`3dmaskave -q -mask ${mydir}/my_mask.nii dti_FA.nii.gz 2>/dev/null` &&
        ad=`3dmaskave -q -mask ${mydir}/my_mask.nii dti_L1.nii.gz 2>/dev/null` &&
        3dcalc -a dti_L2.nii.gz -b dti_L3.nii.gz -expr "(a + b) / 2" \
            -prefix ${mydir}/RD.nii 2>/dev/null &&
        rd=`3dmaskave -q -mask ${mydir}/my_mask.nii ${mydir}/RD.nii 2>/dev/null` &&
        nvox=`3dBrickStat -count -non-zero ${mydir}/my_mask.nii 2>/dev/null` &&
        echo ${m},${fa},${ad},${rd},${nvox} >> $mean_props;
        rm ${mydir}/*nii;
    else
        echo ${m},NA,NA,NA,NA >> $mean_props;
    fi
done
```

Our data looks a bit noisier than PNC, but at least we don't have the binormal
pattern anymore!

![](images/2019-03-01-18-31-32.png)

And there are some bad masks there as well, along with other bad DTI estimates
like PNC, so it will likely get better.

Let's also grab the movement variables:

```bash
out_fname=~/tmp/mvmt_report.csv;
echo "id,Noutliers,PROPoutliers,NoutVolumes,norm.trans,norm.rot,RMS1stVol,RMSprevVol" > $out_fname;
for m in `cat ~/tmp/myids.txt`; do
    echo 'Collecting metrics for' $m;
    if [ -e ${m}/eddy_s2v_unwarped_images.eddy_outlier_report ]; then
        noutliers=`cat ${m}/eddy_s2v_unwarped_images.eddy_outlier_report | wc -l`;
        # figuring out the percetnage of total slices the outliers represent
        nslices=`tail ${m}/eddy_s2v_unwarped_images.eddy_outlier_map | awk '{ print NF; exit } '`;
        nvol=`cat ${m}/dwi_comb_cvec.dat | wc -l`;
        let totalSlices=$nslices*$nvol;
        pctOutliers=`echo "scale=4; $noutliers / $totalSlices" | bc`;
        # figuring out how many volumes were completely removed (row of 1s)
        awk '{sum=0; for(i=1; i<=NF; i++){sum+=$i}; sum/=NF; print sum}' \
            ${m}/eddy_s2v_unwarped_images.eddy_outlier_map > outlier_avg.txt;
        nOutVols=`grep -c -e "^1$" outlier_avg.txt`;
        1d_tool.py -infile ${m}/eddy_s2v_unwarped_images.eddy_movement_over_time \
            -select_cols '0..2' -collapse_cols euclidean_norm -overwrite \
            -write trans_norm.1D;
        trans=`1d_tool.py -infile trans_norm.1D -show_mmms | \
            tail -n -1 | awk '{ print $8 }' | sed 's/,//'`;
        1d_tool.py -infile ${m}/eddy_s2v_unwarped_images.eddy_movement_over_time \
            -select_cols '3..5' -collapse_cols euclidean_norm -overwrite \
            -write rot_norm.1D;
        rot=`1d_tool.py -infile rot_norm.1D -show_mmms | \
            tail -n -1 | awk '{ print $8 }' | sed 's/,//'`;
        1d_tool.py -infile ${m}/eddy_s2v_unwarped_images.eddy_movement_rms \
            -show_mmms > mean_rms.txt;
        vol1=`head -n +2 mean_rms.txt | awk '{ print $8 }' | sed 's/,//'`;
        pvol=`tail -n -1 mean_rms.txt | awk '{ print $8 }' | sed 's/,//'`;
    else
        echo "Could not find outlier report for $m"
        noutliers='NA';
        pctOutliers='NA';
        nOutVols='NA';
        trans='NA';
        rot='NA';
        vol1='NA';
        pvol='NA';
    fi;
    echo $m, $noutliers, $pctOutliers, $nOutVols, $trans, $rot, $vol1, $pvol >> $out_fname;
done
```

# 2019-03-04 10:07:45

Time fit autoPtx2... I'm guessing we'll have many fails, but it's worth running
everything and then just checking for errors:

xaa: g
xab: p
xac: j
xad: g
xae: p
xaf: j
xag: g
xah: p
xai: j
xaj: g
xak: g
xal: p
xam: j
xan: p

And while we wait on that, let's go ahead and make the QC images:

```bash
for m in `cat xab`; do
    bash ~/research_code/dti/fdt_TBSS_and_QC.sh /data/NCR_SBRB/dti_fdt/preproc/${m};
done
```

And we make sure everything ran:

```bash
for m in `cat ~/tmp/myids.txt`; do 
    if [ ! -e preproc/${m}/QC/sse.sag.png ]; then
        echo $m;
    fi;
done
```

Now, we just need to copy everything:

```bash
cd /data/NCR_SBRB/dti_fdt/summary_QC/
for m in `cat ../myids.txt`; do
    echo $m;
    cp ../preproc/${m}/QC/FA_transform.axi.png transform/${m}.axi.png
    cp ../preproc/${m}/QC/FA_transform.sag.png transform/${m}.sag.png
    cp ../preproc/${m}/QC/FA_transform.cor.png transform/${m}.cor.png

    cp ../preproc/${m}/QC/DEC_qc_dec_sca07.axi.png DEC/${m}.axi.png
    cp ../preproc/${m}/QC/DEC_qc_dec_sca07.sag.png DEC/${m}.sag.png
    cp ../preproc/${m}/QC/DEC_qc_dec_sca07.cor.png DEC/${m}.cor.png

    cp ../preproc/${m}/QC/sse.axi.png SSE/${m}.axi.png
    cp ../preproc/${m}/QC/sse.cor.png SSE/${m}.cor.png
    cp ../preproc/${m}/QC/sse.sag.png SSE/${m}.sag.png
done
```

# 2019-03-07 09:24:34

Let's make sure everything ran fine:

```bash
cd /data/NCR_SBRB/dti_fdt/tracts
for s in `cat ~/tmp/myids.txt`; do
    if [ ! -e ${s}/fmi/tracts/tractsNorm.nii.gz ]; then
        echo $s;
    fi;
done
```

Then, it should just be a matter of weighting the property maps by tractNorm:

```bash
mydir=/lscratch/${SLURM_JOBID}/
weighted_tracts=~/tmp/ncr_weighted_tracts.csv;
row="id";
for t in `cut -d" " -f 1 /data/NCR_SBRB/software/autoPtx/structureList`; do
    for m in fa ad rd; do
        row=${row}','${t}_${m};
    done
done
echo $row > $weighted_tracts;
for m in `cat ~/tmp/myids.txt`; do
    echo $m;
    row="${m}";
    cd /data/NCR_SBRB/dti_fdt/preproc/$m &&
    for t in `cut -d" " -f 1 /data/NCR_SBRB/software/autoPtx/structureList`; do
        if [ -e ../../tracts/${m}/${t}/tracts/tractsNorm.nii.gz ]; then
            # tract mask is higher dimension!
            3dresample -master dti_FA.nii.gz -prefix ${mydir}/mask.nii \
                -inset ../../tracts/${m}/${t}/tracts/tractsNorm.nii.gz \
                -rmode NN -overwrite &&
            fa=`3dmaskave -q -mask ${mydir}/mask.nii dti_FA.nii.gz 2>/dev/null` &&
            ad=`3dmaskave -q -mask ${mydir}/mask.nii dti_L1.nii.gz 2>/dev/null` &&
            3dcalc -a dti_L2.nii.gz -b dti_L3.nii.gz -expr "(a + b) / 2" \
                -prefix ${mydir}/RD.nii 2>/dev/null &&
            rd=`3dmaskave -q -mask ${mydir}/mask.nii ${mydir}/RD.nii 2>/dev/null` &&
            row=${row}','${fa}','${ad}','${rd};
            rm ${mydir}/*nii;
        else
            row=${row}',NA,NA,NA';
        fi;
    done
    echo $row >> $weighted_tracts;
done
```

# 2019-03-28 10:32:15

I discovered that our data needs not only the standard flip from FSL, but also a flip_x, based on looking at the vectors. So, let's redo it:

```bash
for s in `cat converted.txt`; do
    cd /data/NCR_SBRB/dti_fdt/${s};
    rm -rf *eddy* FA DEC* origdata dti*;
done

cd /data/NCR_SBRB/dti_fdt
rm -rf swarm.fdt;
for m in `cat xak xal xam xan xao`; do
    echo "bash ~/research_code/dti/fdt_ncr_eddy.sh /data/NCR_SBRB/dti_fdt/${m}" >> swarm.fdt;
done;
swarm -g 4 --job-name fdt --time 4:00:00 -f swarm.fdt --partition gpu \
    --logdir trash_fdt --gres=gpu:k80:2

data='';
for m in `cat xan xao`; do
    data=$data' '${m}/data.nii.gz;
done
/data/NCR_SBRB/software/autoPtx/autoPtx_1_preproc $data;
```

eddy
g: a, b, c, d, e
p: f, g, h, i, j
j: k, l, m, n, o

ptx2
g: a, b, c, d
p: f, g, h, i
j: k, l, m