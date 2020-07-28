# 2020-07-23 13:10:02

Let's see if we can come up with similar numbers, starting with more recent
files.

Whatever I do, I should include the 3 columns in the same file I sent Philip:
in393, in362, hasPRS.

# 2020-07-24 20:30:21

I started from some files Philip had sent me (362_399).

Just some notes on how I'm re-creating the JHU tracts:

```bash
#bw
mydir=~/data/baseline_prediction/dti/
weighted_tracts=jhu_tracts_291.csv;
cd $mydir
row="id";
for t in ATR CST cin_cin cin_hip CC IFO ILF SLF unc SLFtemp; do
    for h in l r; do
        m=fa;
        row=${row}','${m}_${t}_${h};
    done;
done
echo $row > $weighted_tracts;
for m in `cat ~/tmp/m2.txt`; do
    echo ${m}
    3dresample -master ./${m}_tensor_diffeo_fa.nii.gz -prefix ./rois.nii \
                -inset ../../JHU_ICBM_tractsThr25_inAging.nii.gz \
                -rmode NN -overwrite 2>/dev/null &&
    row="${m}";
    for t in `seq 1 20`; do
        3dcalc -a rois.nii -expr "amongst(a, $t)" -prefix mask.nii \
            -overwrite 2>/dev/null &&
        fa=`3dmaskave -q -mask mask.nii ${m}_tensor_diffeo_fa.nii 2>/dev/null`;
        row=${row}','${fa};
    done
    echo $row >> $weighted_tracts;
done
```

# 2020-07-26 20:59:26

After recreating the data, here are some more ML results:

```bash
my_dir=~/data/baseline_prediction
cd $my_dir
my_script=~/research_code/baseline_prediction/twoClass_ROC_splitFirst.R;
sx="categ_all_lm";
imp=dti;
clf=slda;
cov=F;
res_file=res362_splitFirst_sexAgeAll_lastAge.csv;
for cs in "worsening improvers" "worsening never_affected" \
            "worsening stable" "improvers never_affected" \
            "improvers stable" "stable never_affected"; do
    Rscript $my_script ${my_dir}/FINAL_DATA_to_gs_JULY_26_2020b.csv $sx $cs $clf $imp $cov $res_file;
done;
```

# 2020-07-27 17:37:22

Let's maybe do some upsampling tests within the LOOCV framework:

```bash
my_dir=~/data/baseline_prediction
cd $my_dir
my_script=~/research_code/baseline_prediction/twoClass_ROC_experiments.R;
sx="categ_all_lm";
imp=dti;
clf=slda;
cov=F;
res_file=res362_LOOCV_up.csv;
for cs in "worsening improvers" "worsening never_affected" \
            "worsening stable" "improvers never_affected" \
            "improvers stable" "stable never_affected"; do
    Rscript $my_script ${my_dir}/FINAL_DATA_to_gs_JULY_26_2020b.csv $sx $cs $clf $imp $cov $res_file;
done;
```