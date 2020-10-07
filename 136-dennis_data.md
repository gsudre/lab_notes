# 2020-10-06 11:02:09

Let's see what we can find on Dennis' data.

```bash
# bw
cd ~/data/dennis/24
module load afni
dcm2niix_afni *
```

Looking at the JSONs, I saw that 24__20191126213428_13 was the DTI. Let's go
with the FDT pipeline for processing because TORTOISE takes forever.

I copied that the scan files to a raw folder to keep it tidy.

```bash
# bw
module load fsl

fslreorient2std 24__20191126213428_13.nii dwi
fslroi dwi b0 0 1
bet b0 b0_brain -m -f 0.2
cp 24__20191126213428_13.bval bvals
cp 24__20191126213428_13.bvec bvecs
# making sure the directions are OK
dtifit --data=dwi --mask=b0_brain --bvecs=bvecs --bvals=bvals --out=tmp_dti
```

Looks fine I think. Let's go on with processing.

```bash
a=''; for i in {1..102}; do 
    a=$a' '1;
done;
echo $a > index.txt
echo "0 1 0 0.05" > acqparams.txt
echo "0 -1 0 0.05" >> acqparams.txt

eddy_openmp --imain=dwi --mask=b0_brain --index=index.txt \
     --acqp=acqparams.txt \
     --bvecs=bvecs --bvals=bvals --fwhm=0 --flm=quadratic \
     --out=eddy_unwarped_images
```

And I'll follow the same steps for the other subjects (at least until I figure
out whether I'm importing them properly)?

Continuing to bedpostX:

```bash
cp eddy_unwarped_images.nii.gz data.nii.gz;
cp eddy_unwarped_images.eddy_rotated_bvecs bvecs;
cp b0_brain_mask.nii.gz nodif_brain_mask.nii.gz;
/data/NCR_SBRB/software/autoPtx/autoPtx_1_preproc data.nii.gz
```

When bedpost is done, run the second part of autoPtx:

```bash
cd ~/data/dennis
mkdir preproc

scan=27
mkdir preproc/$scan
mkdir preproc/${scan}.bedpostX
execPath=/data/NCR_SBRB/software/autoPtx
structures=$execPath/structureList
track=$execPath/trackSubjectStruct

cp -r $scan/raw/preproc/* preproc/$scan/
cp -r ${scan}/raw/preproc.bedpostX/* preproc/${scan}.bedpostX/
while read structstring; do
    struct=`echo $structstring | awk '{print $1}'`
    nseed=`echo $structstring | awk '{print $2}'`
#    echo $struct;
    $track $scan $struct $nseed 2>&1 & #> /dev/null &
done < $structures
```

# 2020-10-07 05:59:28

Time to extract the tract values. 

```bash
# bw
module load afni

mydir=/lscratch/${SLURM_JOBID}/
weighted_tracts=~/tmp/pnc_weighted_tracts.csv;
row="id";
for t in `cut -d" " -f 1 /data/NCR_SBRB/software/autoPtx/structureList`; do
    for m in fa ad rd; do
        row=${row}','${t}_${m};
    done
done
echo $row > $weighted_tracts;
for m in 24 26 27; do
    echo $m;
    row="${m}";
    cd ~/data/dennis/preproc/$m &&
    for t in `cut -d" " -f 1 /data/NCR_SBRB/software/autoPtx/structureList`; do
        if [ -e ../../tracts/${m}/${t}/tracts/tractsNorm.nii.gz ]; then
            # tract mask is higher dimension!
            3dresample -master dti_FA.nii.gz -prefix ${mydir}/mask.nii \
                -inset ../../tracts/${m}/${t}/tracts/tractsNorm.nii.gz \
                -rmode NN -overwrite &&
            nvox=`3dBrickStat -count -non-zero ${mydir}/mask.nii 2>/dev/null` &&
            if [ $nvox -gt 0 ]; then
                fa=`3dmaskave -q -mask ${mydir}/mask.nii dti_FA.nii.gz 2>/dev/null` &&
                ad=`3dmaskave -q -mask ${mydir}/mask.nii dti_L1.nii.gz 2>/dev/null` &&
                3dcalc -a dti_L2.nii.gz -b dti_L3.nii.gz -expr "(a + b) / 2" \
                    -prefix ${mydir}/RD.nii 2>/dev/null &&
                rd=`3dmaskave -q -mask ${mydir}/mask.nii ${mydir}/RD.nii 2>/dev/null` &&
                row=${row}','${fa}','${ad}','${rd};
            else
                row=${row}',NA,NA,NA';
            fi;
        else
            row=${row}',NA,NA,NA';
        fi;
    done
    echo $row >> $weighted_tracts;
done
```

And I copied the file, with pnc_, to the dennis folder. Let's generate some more QC files, just because...

```bash
m=24;
cd ~/data/dennis/preproc/$m

# make directionality encoded QC pictures, checking that eroded FA looks fine.
# Note that the dtifit results, and the std alignment were run by autoPtx!
fat_proc_decmap -in_fa dti_FA.nii.gz -in_v1 dti_V1.nii.gz \
    -mask nodif_brain_mask.nii.gz -prefix DEC

# apply the transform calculated by autoPtx to a few maps. code copied from 
# tbss_non_fa
for f in FA L1 L2 L3 MD MO; do
    echo Warping $f;
    applywarp -i dti_${f} -o ${f}_in_FMRIB58_FA_1mm \
        -r $FSLDIR/data/standard/FMRIB58_FA_1mm -w nat2std_warp
done

# make transformation QC figure: FSL template as the edges
@snapshot_volreg FA_in_FMRIB58_FA_1mm.nii.gz \
    $FSLDIR/data/standard/FMRIB58_FA_1mm.nii.gz \
    QC/FA_transform;

# make QC images for standard errors. Here we set our color scale to have 95th
# percentile of all errors. Meaning, more red = bigger error.
@chauffeur_afni                             \
    -ulay  data.nii.gz                       \
    -olay  dti_sse.nii.gz                          \
    -opacity 5                              \
    -pbar_posonly   \
    -cbar Spectrum:red_to_blue              \
    -set_subbricks 0 0 0     \
    -prefix   QC/sse              \
    -montx 6 -monty 6                       \
    -set_xhairs OFF                         \
    -label_mode 1 -label_size 3             \
    -thr_olay 0 \
    -func_range_perc_nz 95 \
    -do_clean
```

And then copied everything back to shaw:/dennis/.

QC images look fine, but forceps minor and major could not be estimated for one
of the subjects. Might need to tweak things a bit there later.