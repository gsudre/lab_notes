# 2019-05-09 15:00:50

We were getting hammered by multiple comparisons in rsFMRI, so let's instead
check for within and acorss network connectivity. We already have all subjects
in the common space, so now we just need to mark the different areas of each
network. We could do this by voxels or by ROI. And that's assuming we want a
within-network metric, otherwise we could just take a network to network
approach.

In any case, we'll need the data ll in the same space. First, let make sure we
have all the anatomical files we need:

```bash
net_dir=/Volumes/Shaw/MR_data_by_maskid/
cd ~/data/heritability_change/fmri_same_space/anat
export OMP_NUM_THREADS=4
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    if [ ! -e anatQQ.$m.nii ]; then
        echo $m;
        mri_convert /Volumes/Shaw/freesurfer5.3_subjects/${m}/mri/orig/001.mgz ./${m}.nii &&
        @SSwarper \
            -input ./${m}.nii \
            -base TT_N27_SSW.nii.gz -subid ${m} &&
        rm ${m}.nii;
    fi;
done;
```

But while that's going on I made sure all the transform QC picture looked
good...

# 2019-05-10 11:12:17

So, the anatomical warping was taking a very long time to run in the interactive
nodes. We'll need to swarm it...

```bash
cd ~/data/heritability_change/fmri_same_space/anat;
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    if [ ! -e anatQQ.$m.nii ]; then
        echo "export OMP_NUM_THREADS=16; cd ~/data/heritability_change/fmri_same_space/anat; @SSwarper -input ${m}.nii -base TT_N27_SSW.nii.gz -subid ${m}" >> swarm.SSwarper;
    fi;
done;
swarm -g 10 -t 16 --job-name SSwarper --time 2:00:00 -f swarm.SSwarper -m afni --partition quick --logdir trash
```

However, it makes more sense ot run all of this in Luke's framework, but it
needs them converted to MNI space instead. Let's compute that, shall we?

```bash
cd ~/data/heritability_change/fmri_same_space/anat_mni;
for maskid in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m=`printf %04d $maskid`;
    # mri_convert /data/NCR_SBRB/freesurfer5.3_subjects/${m}/mri/orig.mgz ${m}.nii.gz;
    echo "export OMP_NUM_THREADS=16; cd ~/data/heritability_change/fmri_same_space/anat_mni; @SSwarper -input ${m}.nii.gz -base MNI152_2009_template_SSW.nii.gz -subid ${m}" >> swarm.SSwarper;
done;
swarm -g 10 -t 16 --job-name SSwarper --time 2:00:00 -f swarm.SSwarper -m afni --partition quick --logdir trash
```

And we need to check on the alignments again. But first, transfer the actual EPI
signal to MNI space:

Now we just do:

```bash
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/fmri_same_space/epi;
for m in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m2=`printf %04d $m`;
    echo $m2;
    afni_dir=/mnt/shaw/MR_data_by_maskid/${m2}/afni/${m2}.rest.subjectSpace.results/;
    3dNwarpApply -nwarp "../anat_mni/anatQQ.${m2}_WARP.nii ../anat_mni/anatQQ.${m2}.aff12.1D" \
        -source ${afni_dir}/errts.${m2}.fanaticor+orig.HEAD \
        -master ../anat_mni/anatQQ.${m2}.nii -dxyz 2.5\
        -overwrite -prefix ${m2}_epi_NL_inMNI.nii;
done
```

Then, compute the sphere averages and correlations. But instead of having all
those intermediate files with the ROI timecourse, why not create one big mask
with all the spheres, and use 3dNetCorr like before? Here' I'll use the same
numbering from the xlsx file in the Power atlas
(https://www.jonathanpower.net/2011-neuron-bigbrain.html). 

```bash
cd /mnt/shaw/Gustavo/desktop_backup/data/heritability_change/fmri_same_space/epi;
3dUndump -prefix spheres.nii -master 0411_epi_NL_inMNI.nii -srad 5 \
    -orient LPI -xyz coords.1D
net_dir=/mnt/shaw/MR_data_by_maskid/
for m in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m2=`printf %04d $m`;
    3dNetCorr                                       \
        -inset ${m2}_epi_NL_inMNI.nii                    \
        -in_rois  spheres.nii                       \
        -prefix  ../../fmri_corr_tables/${m2}_power                           \
        -fish_z;
done
```

# 2019-05-13 14:55:17

I created a streamlined version of Neuro_consensus_264 as a csv, only containing
the important columns.

Some of the IDs failed to create the QC image, so I'll just use the same script
I used for FDT so that we can have uniform pictures on all scans:

```bash
module load afni
cd ~/data/heritability_change/fmri_same_space/anat_mni;
for m in `cut -d"," -f 1 ../../rsfmri_3min_assoc_n462.csv | tail -n +2`; do
    m2=`printf %04d $m`;
    @chauffeur_afni                             \
        -ulay  /usr/local/apps/afni/current/linux_centos_7_64/MNI152_2009_template.nii.gz                     \
        -olay  anatQQ.${m2}.nii                         \
        -ulay_range 0% 150%                     \
        -func_range_perc 50                     \
        -pbar_posonly                           \
        -cbar "red_monochrome"                  \
        -opacity 8                              \
        -prefix   ${m2}_QC              \
        -montx 3 -monty 3                       \
        -set_xhairs OFF                         \
        -label_mode 1 -label_size 3             \
        -do_clean
done
```

I then created a function that collectes the Power parcellation from 3dNetCorr,
but also combines the results into the specific networks:

```r
source('~/research_code/fmri/combine_3dNetCorr_grids.R')
```

# 2019-05-14 09:43:59

Note that for now I'm removing the mask ids that don't have all 264 spheres. It
could easily be that their transformation is bad, or that the prescription
didn't cover those regions. In any case, I'll proceed with that for now (we have
458 mask ids instead of 462 now).

So, time to massage the data that goes into SOLAR:

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
m2 = read.csv('~/data/heritability_change/rsfmri_3min_assoc_n462.csv')
clin = read.csv('~/data/heritability_change/clinical_03132019.csv')
df = mergeOnClosestDate(m2, clin, unique(m2$Medical.Record...MRN),
                         x.date='record.date.collected...Scan',
                         x.id='Medical.Record...MRN')
b = read.csv('/Volumes/Shaw/Gustavo/desktop_backup/data/heritability_change/fmri_corr_tables/pearsonZ_3min_n462_power.csv')
var_names = colnames(b)[2:ncol(b)]
df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id', all.x=F)

# make sure we still have two scans for everyone
rm_subjs = names(which(table(df2$Medical.Record...MRN)<2))
rm_me = df2$Medical.Record...MRN %in% rm_subjs
df2 = df2[!rm_me, ]

library(MASS)
mres = df2
mres$SX_HI = as.numeric(as.character(mres$SX_hi))
mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))
for (t in var_names) {
    fm_str = sprintf('%s ~', t)
    fm_str = paste(fm_str,
                   'enormGoodTRs_fmri01 + I(enormGoodTRs_fmri01^2)')
    res.lm <- lm(as.formula(fm_str), data = mres, na.action=na.exclude)
    step <- stepAIC(res.lm, direction = "both", trace = F)
    mres[, t] = residuals(step)
}
res = c()
for (s in unique(mres$Medical.Record...MRN)) {
    idx = which(mres$Medical.Record...MRN == s)
    row = c(s, unique(mres[idx, 'Sex']))
    for (t in var_names) {
        if (sum(is.na(mres[idx, t])) > 0) {
            # if any of the variables is NA, make the slope NA
            row = c(row, NA)
        } else {
            fm_str = sprintf('%s ~ age_at_scan', t)
            fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
            row = c(row, coefficients(fit)[2])
        }
    }
    for (t in c('SX_inatt', 'SX_HI')) {
        fm_str = sprintf('%s ~ age_at_scan', t)
        fit = lm(as.formula(fm_str), data=mres[idx, ], na.action=na.exclude)
        row = c(row, coefficients(fit)[2])
    }
    # grabbing inatt and HI at baseline
    base_DOA = which.min(mres[idx, 'age_at_scan'])
    row = c(row, mres[idx[base_DOA], 'SX_inatt'])
    row = c(row, mres[idx[base_DOA], 'SX_HI'])
    # DX1 is DSMV definition, DX2 will make SX >=4 as ADHD
    if (mres[idx[base_DOA], 'age_at_scan'] < 16) {
        if ((row[length(row)] >= 6) || (row[length(row)-1] >= 6)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    } else {
        if ((row[length(row)] >= 5) || (row[length(row)-1] >= 5)) {
            DX = 'ADHD'
        } else {
            DX = 'NV'
        }
    }
    if ((row[length(row)] >= 4) || (row[length(row)-1] >= 4)) {
        DX2 = 'ADHD'
    } else {
        DX2 = 'NV'
    }
    row = c(row, DX)
    row = c(row, DX2)
    res = rbind(res, row)
    print(dim(res)[1])
}
colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI',
                                              'inatt_baseline',
                                              'HI_baseline', 'DX', 'DX2'))
write.csv(res, file='~/data/heritability_change/rsfmri_3min_assoc_n227_netSlopesZ.csv',
          row.names=F, na='', quote=F)

res_clean = res
# and remove outliers
for (t in var_names) {
    mydata = as.numeric(res[, t])
    # identifying outliers
    ul = mean(mydata) + 3 * sd(mydata)
    ll = mean(mydata) - 3 * sd(mydata)
    bad_subjs = c(which(mydata < ll), which(mydata > ul))

    # remove within-tract outliers
    res_clean[bad_subjs, t] = NA
}
write.csv(res_clean, file='~/data/heritability_change/rsfmri_3min_assoc_n227_slopesCleanZ.csv',
          row.names=F, na='', quote=F)

# and make sure every family has at least two people
good_nuclear = names(table(m2$Nuclear.ID...FamilyIDs))[table(m2$Nuclear.ID...FamilyIDs) >= 4]
good_extended = names(table(m2$Extended.ID...FamilyIDs))[table(m2$Extended.ID...FamilyIDs) >= 4]
keep_me = c()
for (f in good_nuclear) {
    keep_me = c(keep_me, m2[which(m2$Nuclear.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
for (f in good_extended) {
    keep_me = c(keep_me, m2[which(m2$Extended.ID...FamilyIDs == f),
                            'Medical.Record...MRN'])
}
keep_me = unique(keep_me)

fam_subjs = c()
for (s in keep_me) {
    fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
}
res2 = res[fam_subjs, ]
res2_clean = res_clean[fam_subjs, ]

write.csv(res2, file='~/data/heritability_change/rsfmri_3min_n111_netSlopesZ.csv',
          row.names=F, na='', quote=F)
write.csv(res2_clean, file='~/data/heritability_change/rsfmri_3min_n111_netSlopesCleanZ.csv',
          row.names=F, na='', quote=F)
write.table(var_names, file='~/data/heritability_change/power.txt',
            row.names=F, col.names=F, quote=F)
```

And I ran the code above for the regular R and Z version. I could potentially
explore only absolute correlation... maybe later.

We then run SOLAR, which should be fast enough locally, but I still use the
voxel-type of calls:

```bash
phen_file=rsfmri_3min_n111_netSlopes
tmp_dir=~/data/heritability_change
solar_dir=~/data/heritability_change
mkdir ${tmp_dir}/${phen_file}
mkdir /tmp/${phen_file}
for vox in `cat ~/data/heritability_change/power.txt`; do
    mkdir /tmp/${phen_file}/${vox};
    cp ${solar_dir}/pedigree.csv ${solar_dir}/procs.tcl ${solar_dir}/${phen_file}.csv /tmp/${phen_file}/${vox}/;
    cd /tmp/${phen_file}/${vox}/;
    ~/Downloads/solar842/solar run_phen_var $phen_file $vox;
    mv /tmp/${phen_file}/${vox}/i_${vox} ${tmp_dir}/${phen_file}/;
done
```

And as usual, I ran it for Clean and Z versions as well.

I'll need to re-run this because I'm having issues with relpairs in this version
of SOLAR. So, results might change a bit if I end up having to remove some
subjects. In any case, let's see what the compiled results look like here...

Only one connection (memory retrieval to ventral attention) was significant.
Whomp-whomp... it was consistent across all 4 modes, but the only thing
significant too. Gonna need to get more subjects, or try other things. At least,
it looks like the results across all 4 are somewhat consistent. Let's choose one
and play with it a bit.

I checked the alignment, and it looks fine. But the results using absolute
values only, negative as nan, or negative as zero didn't help. They were
actually worse.

There are a couple more things I want to try: not remove movement right away,
and MELODIC. Also, segregation, but I'd expect some hint of results without
using segregation to begin with... I also don't think movement is going to
affect it that much, because we're using stepAIC, so it wouldn't remove
moevement if it didn't need to. So, let's try MELODIC first.


# TODO

* try not removing movement, and then just checking later if the values in any
  significant variables are related to movement?
* try segregation
* ICA networks from MELODIC?
