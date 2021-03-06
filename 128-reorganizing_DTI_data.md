# 2020-07-31 11:34:55

Let's do some light re-organizing of the DTI data processed in DTI-TK. I saw in
a previous Evernote note "DTI robust TSA" that showed the pipeline used in
analysis_may2017 was just a regular transformation to the aging template. So,
let's use that same code to re-run the scans that are missing thay output file:

```bash
dataDir=/mnt/shaw/sudregp/dtitk_processing/tortoise_output/
maskidDir=/mnt/shaw/sudregp/MR_data_by_maskID/
maskidFile=~/tmp/dti.txt

cd $dataDir
for m in `cat $maskidFile`; do
    echo $m;
    dtitkFolder=${maskidDir}/${m}/edti_proc/edti_DMC_DR_R1_SAVE_DTITK/;
    if [ ! -d $dtitkFolder ]; then
        dtitkFolder=${maskidDir}/${m}/edti_proc/edti_DMC_R1_SAVE_DTITK;
        tensorFile=edti_DMC_R1_tensor.nii;
    else
        tensorFile=edti_DMC_DR_R1_tensor.nii;
    fi;
    if [ ! -d $dtitkFolder ]; then
        echo $dtitkFolder "does not exist. Skipping...";
    else
        cp ${dtitkFolder}/${tensorFile} ${dataDir}/${m}_tensor.nii;
    fi;
done
```

Now let's move to a more organized folder the stuff we need:

```bash
source_dir=/mnt/shaw/sudregp/dti_robust_tsa/analysis_may2017/;
cd /mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate;
for m in `cat ~/tmp/dti.txt`; do
    echo ${m};
    if [ -e ${source_dir}/${m}_tensor_diffeo.nii.gz ]; then
        mv ${source_dir}/${m}_* .;
    else
        echo $m >> ~/tmp/no_dtitk.txt;
    fi;
done
```

Let's run DTI-TK for everyone that we have tensors:

```bash
cd /mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate;
template=../ixi_aging_template_v3.0/template/ixi_aging_template

for m in `cat ~/tmp/dti.txt`; do
    if [ -e ../tortoise_output/${m}_tensor.nii ]; then
        if [ ! -e ./${m}_tensor_diffeo.nii.gz ]; then
            echo $m;
            echo ../tortoise_output/${m}_tensor.nii > ~/tmp/rm${m};
            dti_rigid_sn ${template}.nii.gz ~/tmp/rm${m} EDS;
            dti_affine_sn ${template}.nii.gz ~/tmp/rm${m} EDS 1;
            echo ../tortoise_output/${m}_tensor_aff.nii.gz > ~/tmp/rm${m}_aff;
            dti_diffeomorphic_sn ${template}.nii.gz ~/tmp/rm${m}_aff \
                ${template}_brain_mask.nii.gz 6 0.002;
            dti_warp_to_template_group ~/tmp/rm${m} ${template}.nii.gz 2 2 2;
            mv ../tortoise_output/${m}*aff* ../tortoise_output/${m}_tensor_* .
        fi;
    fi;
done
```

I copied and modified the code above from run_tsa.sh, which can be used to run
it in the cluster. Here it's just a handful of maskids, so it doesn't take very
long to run it in ncrshell01.

Finally, let's create the property files for everyone that doesn't have them:

```bash
cd /mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate;
for m in `cat ~/tmp/dti.txt`; do
    if [ -e ./${m}_tensor_diffeo.nii.gz ]; then
        if [ ! -e ./${m}_tensor_diffeo_rd.nii.gz ]; then
            TVtool -in ${m}_tensor_diffeo.nii.gz -fa;
            TVtool -in ${m}_tensor_diffeo.nii.gz -ad;
            TVtool -in ${m}_tensor_diffeo.nii.gz -rd;
        fi;
    fi;
done
```

For the reliability scans, I'm just copying the first session number as the mask
id scan. Copying the TORTOISE output tensors manually.

And some quick code to get the JHU tracts:

```bash
#ncrshell01
mydir=/mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate
weighted_tracts=jhu_tracts_183.csv;
cd $mydir
row="id";
for t in ATR CST cin_cin cin_hip CC IFO ILF SLF unc SLFtemp; do
    for h in l r; do
        for m in fa ad rd; do
            row=${row}','${m}_${t}_${h};
        done;
    done;
done
echo $row > $weighted_tracts;
for m in `cat ~/tmp/m2.txt`; do
    echo ${m}
    3dresample -master ./${m}_tensor_diffeo_fa.nii.gz -prefix ./rois.nii \
                -inset ../JHU_ICBM_tractsThr25_inAging.nii.gz \
                -rmode NN -overwrite 2>/dev/null &&
    row="${m}";
    for t in `seq 1 20`; do
        3dcalc -a rois.nii -expr "amongst(a, $t)" -prefix mask.nii \
            -overwrite 2>/dev/null &&
        fa=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_fa.nii 2>/dev/null`;
        ad=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_ad.nii 2>/dev/null`;
        rd=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_rd.nii 2>/dev/null`;
        row=${row}','${fa}','${ad}','${rd};
    done
    echo $row >> $weighted_tracts;
    rm -rf rois.nii mask.nii;
done
```

And we can pull out the DTI-TK tracts this way for the same IDs:

```bash
cd /mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate
for m in `cat ~/tmp/m2.txt`; do
    echo "${m}_tensor_diffeo.nii.gz" >> ~/tmp/m2_diffeo.txt;
done
tsa_sampling ~/tmp/m2_diffeo.txt ../ixi_aging_template_v3.0/tsa/ mean
python3 ~/research_code/lab_mgmt/convert_dti_sampling.py
Rscript ~/research_code/dti/compile_tract_table.R
```