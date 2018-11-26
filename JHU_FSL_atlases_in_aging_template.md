# 2018-11-26 13:58:07

Just a few notes on how I converted the JHU labels and tracts to the aging
template. I checked them visually using a mean_FA over several subjects I had
lying around
(/Volumes/Shaw/dti_robust_tsa/analysis_may2017/mean_subjects_730_fa.nii.gz) and
it they looked acceptable, so I'll go with those.

```bash
flirt -in /usr/local/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz -ref ~/data/ixi_aging_template_v3.0/template/ixi_aging_template_brain_mask.nii.gz -out MNI152_brain_mask_IN_agingTemplate.nii.gz -omat MNI152_to_aging.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12 -interp trilinear

flirt -in /usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz -ref ixi_aging_template_v3.0/template/ixi_aging_template.nii.gz -out JHU_ICBM_labels_inAging.nii.gz -applyxfm -init MNI152_to_aging.mat -interp nearestneighbour

flirt -in /usr/local/fsl/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr25-1mm.nii.gz -ref ixi_aging_template_v3.0/template/ixi_aging_template.nii.gz -out JHU_ICBM_tractsThr25_inAging.nii.gz -applyxfm -init MNI152_to_aging.mat -interp nearestneighbour
```

Then, to create the tables we need:

```bash
awk '{ FS=">"; if ( NR > 16 && NR < 65 ) { print $2 } }' /usr/local/fsl/data/atlases/JHU-labels.xml | sed -e "s/<\/label//g" > tmp.txt;
i=1; while read l; do echo $i $l >> JHU_labels.txt; let i=$i+1; done < tmp.txt

awk '{ FS=">"; if ( NR > 15 && NR < 36 ) { print $2 } }' /usr/local/fsl/data/atlases/JHU-tracts.xml | sed -e "s/<\/label//g" > tmp.txt
i=1; while read l; do echo $i $l >> JHU_tracts.txt; let i=$i+1; done < tmp.txt
```

I also checked the integer values in fsleyes and they make sense.

Finally, we can get the averages with something like this (for FA, in this case):

```bash
# caterpie
mydir=/mnt/shaw/dti_robust_tsa/analysis_may2017/
# selects one of the labels
3dcalc -a ~/data/JHU_ICBM_labels_inAging.nii.gz -prefix mymask.nii -expr 'amongst(a, 1)';
while read s; do 
    echo $s;
    myfile=${mydir}/${s}_tensor_diffeo_fa.nii.gz
    if [ ! -e $myfile ]; then
        echo "Trying to create DTI property file for $s"
        /usr/local/neuro/dti-tk/dtitk-2.3.1-Linux-x86_64/bin/TVtool \
            -in $mydir/${s}_tensor_diffeo.nii.gz -fa \
            -out $mydir/${s}_tensor_diffeo_fa.nii.gz;
    fi;
    val=`3dmaskave -q -mask mymask.nii ${mydir}/${s}_tensor_diffeo_fa.nii.gz`;
    echo ${s},${val} >> ~/tmp/mean_fa.txt;
done < ~/tmp/dmaskids3.txt
```

I'm keeping the transformed tract maps in ~/data for now.