# 2019-07-08 16:08:42

We first need to replace the variables in the file Philip sent by the voxels. In
other words, all we need to do is replace Y* by voxel data. The data needs to be
residualized after motion and Z scored.

So, let's first extract the list of mask ids and grab the voxel-level data:

```R
long<-read.csv('/Volumes/Shaw/cross_lag/ctsem/for_gustavo/LONG_file_for_gustavo_3obs.csv', T)
m = sapply(long$maskid, function(x) {sprintf('%04d', x)})
write.table(m, file='/Volumes/Shaw/cross_lag/ctsem/for_gustavo/maskids_long.txt', row.names=F, col.names=F, quote=F)
```

I couldn't find the original .nii.gz diffeos that MArine used, so I'll just use
the dataframes she was using for voxelwise data. I'll later need to find the
templates used to extract the voxels, so I can put everything back as nifti.


```bash
# caterpie
cd /mnt/shaw/cross_lag/ctsem/for_gustavo/
sed "s/$/_tensor_diffeo.nii.gz/" maskids_long.txt > subjs_tensor.txt;
sed "s/^/\/mnt\/shaw\/dti_robust_tsa\/analysis_may2017\//" subjs_tensor.txt > subjs_tensor2.txt;








TVMean -in subjs_tensor2.txt -out mean_final_high_res.nii.gz
TVtool -in mean_final_high_res.nii.gz -fa
tbss_skeleton -i mean_final_high_res -o mean_fa_skeleton
flirt -in /usr/local/neuro/fsl/data/standard/FMRIB58_FA_1mm.nii.gz \
    -ref mean_final_high_res_fa.nii.gz \
    -out FMRIB58_FA_IN_groupTemplate.nii.gz -omat FMRIB58_to_group.mat \
    -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 \
    -searchrz -90 90 -dof 12 -interp trilinear

flirt -in /usr/local/neuro/fsl/data/standard/FMRIB58_FA-skeleton_1mm.nii.gz \
    -ref mean_final_high_res_fa.nii.gz \
    -out FMRIB58_FA-skeleton_inGroup.nii.gz -applyxfm \
    -init FMRIB58_to_group.mat -interp nearestneighbour

3dcalc -a FMRIB58_FA-skeleton_inGroup.nii.gz -prefix fa_skeleton_mask.nii.gz \
    -expr 'step(a-.2)'

dti_dir=/mnt/shaw/dti_robust_tsa/analysis_may2017/
mkdir dti_voxels
for maskid in `cat maskids_566.txt`; do
     3dmaskdump -mask fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_fa.txt ${dti_dir}/${maskid}_tensor_diffeo_fa.nii.gz;
     3dmaskdump -mask fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_ad.txt ${dti_dir}/${maskid}_tensor_diffeo_ad.nii.gz;
     3dmaskdump -mask fa_skeleton_mask.nii.gz -o dti_voxels/${maskid}_rd.txt ${dti_dir}/${maskid}_tensor_diffeo_rd.nii.gz;
done
```

Then we construct the data files in R:

```r
maskids = read.table('/mnt/shaw/dti_robust_tsa/heritability/maskids_566.txt')[, 1]
nvox=14681
for (m in c('fa', 'ad', 'rd')) {
     print(m)
     dti_data = matrix(nrow=length(maskids), ncol=nvox)
     for (s in 1:nrow(dti_data)) {
          a = read.table(sprintf('/mnt/shaw/dti_robust_tsa/heritability/dti_voxels/%04d_%s.txt',
                                 maskids[s], m))
          dti_data[s, ] = a[,4]
     }
     dti_data = cbind(maskids, dti_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%05d', d)))
     colnames(dti_data) = cnames
     write.csv(dti_data, file=sprintf('/mnt/shaw/dti_robust_tsa/heritability/dti_%s_voxelwise_n566_03052019.csv', m), row.names=F)
}
```

