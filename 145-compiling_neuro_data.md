# 2020-11-09 20:27:15

Since we had this long break because of COVID, all Freesurfer and DTI data
finished being processed. It makes sense to create a big csv with all the data
in it.

## MPRAGE

We start from Labmatrix to get the best mprage QC, SID, and age_scan. We also
make sure to only get the scans that were processed.

Then:

```bash
# ncrshell01
bash ~/research_code/lab_mgmt/get_freesurfer_roi_data.sh ~/tmp/ids10.txt
```

Then we organize everything in Excel, and to create the summary metrics we do:

```r
rois = read.csv('~/research_Code/REGIONAL_ANALYSES_FREESURFER.csv')
brain_data = read.csv('~/data/all_freesurfer.csv')
brain_vars = colnames(brain_data)[14:287]
for (part in c('sublobar', 'lobar', 'theoryDriven')) {
    for (roi in unique(rois[, part])) {
        labels = rois[which(rois[, part]==roi), 'region']
        to_avg = c()
        for (l in labels) {
            to_avg = c(to_avg,
                    brain_vars[grepl(brain_vars, pattern=sprintf("^%s", l))])
        }
        # only use variable if it's selected initially and defined
        if (length(to_avg) > 0 && sum(is.na(brain_data[, to_avg])) == 0 &&
            nchar(roi) > 0) {
            if (length(to_avg) == 1) {
                brain_data[, sprintf('%s_%s', part, roi)] = brain_data[, to_avg]
            } else {
                brain_data[, sprintf('%s_%s', part, roi)] = rowMeans(brain_data[, to_avg])
            }
        }
    }
}
write.csv(brain_data, file='~/data/all_freesurfer_11122020.csv', row.names=F, quote=F)
```

## DTI

We proceed the same way we did for MPRAGE by getting data for DTI in Labmatrix.
Except that here I removed everything with QualityControl==4, as they were not
processed to the end. If that number is needed in the future one would have to
go to Labmatrix.

Then we need to get the movement metrics, DTITK, and JHU tracts.

```bash
# ncrshell01
cd /mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate
rm ~/tmp/ids1_diffeo.txt;
for m in `cat ~/tmp/ids1`; do
    echo "${m}_tensor_diffeo.nii.gz" >> ~/tmp/ids1_diffeo.txt;
done
tsa_sampling ~/tmp/ids1_diffeo.txt ../ixi_aging_template_v3.0/tsa/ mean
# make sure the script is set to use the correct subject file
python3 ~/research_code/lab_mgmt/convert_dti_sampling.py
Rscript ~/research_code/dti/compile_tract_table.R
```

And we also compiled the ones for JHU_tracts:

```bash
#ncrshell01
mydir=/mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate
weighted_tracts=jhu_tracts_1831.csv;
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
for m in `cat ~/tmp/ids1`; do
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

# 2021-03-22 10:50:56

Let's also add a mean_fa, mean_ad, and mean_rd column:

```bash
#ncrshell01
mydir=/mnt/shaw/sudregp/dtitk_processing/tortoise_dtitk_crossSec_agingTemplate
weighted_tracts=mean_tracts_1831.csv;
cd $mydir
row="id";
for m in fa ad rd; do
    row=${row}','mean_${m};
done
echo $row > $weighted_tracts;
for m in `cat ~/tmp/ids1`; do
    echo ${m}
    row="${m}";
    # not skeletonizing, just simple threshold
    3dcalc -a ./${m}_tensor_diffeo_fa.nii.gz -prefix mask.nii -overwrite \
        -expr "step(a-.2)" 2>/dev/null &&
    fa=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_fa.nii 2>/dev/null`;
    ad=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_ad.nii 2>/dev/null`;
    rd=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_rd.nii 2>/dev/null`;
    row=${row}','${fa}','${ad}','${rd};
    echo $row >> $weighted_tracts;
    rm -rf mask.nii;
done
```
