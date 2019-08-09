# 2019-08-01 16:26:50

The 36P-despike pipeline seems to keep all 126 TRs, and it does a decent job
removing movement correlation overall. According to its paper
(https://jamanetwork.com/journals/jamapsychiatry/fullarticle/2734860):

```
"Task-free functional images were processed using a top-performing pipeline for removal of motion-related artifact.25 Preprocessing steps included (1) correction for distortions induced by magnetic field inhomogeneities using FSL’s FUGUE utility, (2) removal of the 4 initial volumes of each acquisition, (3) realignment of all volumes to a selected reference volume using MCFLIRT16 (4) removal of and interpolation over intensity outliers in each voxel’s time series using AFNI’s 3DDESPIKE utility, (5) demeaning and removal of any linear or quadratic trends, and (6) co-registration of functional data to the high-resolution structural image using boundary-based registration.26 The artifactual variance in the data was modelled using a total of 36 parameters, including the 6 framewise estimates of motion, the mean signal extracted from eroded white matter and cerebrospinal fluid compartments, the mean signal extracted from the entire brain, the derivatives of each of these 9 parameters, and quadratic terms of each of the 9 parameters and their derivatives. Both the BOLD-weighted time series and the artifactual model time series were temporally filtered using a first-order Butterworth filter with a passband between 0.01 and 0.08 Hz. Subjects included in this analysis had low motion as measured by mean frame wise displacement, specifically mean relative frame wise displacement less then 2.5 mm."
```

So, it's almost the same as 36P, but it does step 4. I actually think 2.5mm will
be too lenient in our dataset. Let's see:

![](images/2019-08-01-18-16-20.png)

That's the plot of all the mean FD across the best 2 scans of every kid. I think
we can even get away with 1mm. Let's run both then. For the MELODIC analysis,
all we need is to know the IDs of the two scans we'll be including. 

```bash
Rscript ~/research_code/fmri/make_aroma_condensed_data_FD.R
```

# 2019-08-02 10:36:29

I just noticed that the recommendation of one of the papers is to remove anyone
with mFD above .2. So, clearly there is no consensus in the field. Let's just
pick one that makes sense with our data (1), and then we can test later if
results hold with more stringent threshold (if necessary).


```bash
# desktop
cd ~/data/heritability_change/
mkdir xcp-36p_despike
cd xcp-36p_despike
mkdir masks;
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
fname=rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv;
awk '{FS=","; print $1}' ../$fname | tail -n +2 > ids_1.txt;
for maskid in `cat ids_1.txt`; do
    m=`printf %04d $maskid`;
    # this file is the target of the sub-???.nii.gz symlink
    3dAutomask -prefix masks/${m}_automask.nii \
        $mydir/sub-${m}/norm/sub-${m}_std.nii.gz;
done
cd masks
3dmask_tool -input ????_automask.nii -prefix ../group_epi_mask_inter.nii -frac 1
3dmask_tool -input ????_automask.nii -prefix ../group_epi_mask_fancy.nii \
    -dilate_input 5 -5 -frac 0.7 -fill_holes
```

Let's then run melodic, which we can do locally:

```bash
cd ~/data/heritability_change/xcp-36p_despike
for maskid in `cat ids_1.txt`; do
    m=`printf %04d $maskid`;
    echo $mydir/sub-${m}/norm/sub-${m}_std.nii.gz >> fd1_epi.txt;
done

melodic -i fd1_epi.txt -o groupmelodic_fancy.ica -v --nobet -m group_epi_mask_fancy.nii --tr=2.5 --report --Oall -a concat;

melodic -i fd1_epi.txt -o groupmelodic_inter.ica -v --nobet -m group_epi_mask_inter.nii --tr=2.5 --report --Oall -a concat;
```

Now we performt the dual regression to get each subject's values for the ICs:

```bash
pipe='fancy';
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_${pipe}.ica
mkdir dual
while read maskid; do
    m=`printf %04d $maskid`;
    echo ${pipe} $m;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d melodic_IC \
        -o dual/dr_stage1_${m}.txt --demean -m ../group_epi_mask_${pipe}.nii;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d dual/dr_stage1_${m}.txt \
        -o dual/dr_stage2_${m} --demean -m ../group_epi_mask_${pipe}.nii --des_norm \
        --out_z=dual/dr_stage2_${m}_Z;
done < ../ids_1.txt
```

Now, it's time to figure out which ICs we are going to use.

# 2019-08-05 10:52:29


We are already in MNI space, like the Yeo networks (I checke dit visually). But
the inter mask looked a bit funky, so I'll stay with the fancy mask for now.

```bash
cd ~/data/heritability_change/xcp-36p_despike
for i in {1..7}; do
    3dcalc -prefix Yeo_liberal_net${i}.nii \
        -a ~/data/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz -expr "amongst(a,${i})";
done
3dTcat -prefix Yeo_liberal_combined.nii Yeo_liberal_net1.nii \
    Yeo_liberal_net2.nii Yeo_liberal_net3.nii \
    Yeo_liberal_net4.nii Yeo_liberal_net5.nii \
    Yeo_liberal_net6.nii Yeo_liberal_net7.nii
3dresample -master groupmelodic_inter.ica/melodic_IC.nii.gz \
    -prefix Yeo_nets.nii -inset Yeo_liberal_combined.nii \
    -rmode NN -overwrite
```

So, let's figure out what are the best matching networks for each mask:

```bash
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_fancy.ica/
3dMatch -inset melodic_IC.nii.gz -refset ../Yeo_nets.nii \
    -mask ../group_epi_mask_fancy.nii -prefix matches -overwrite
cat matches_REF_coeff.vals
```

Keep in mind that the code is:

```
0: visual
1: somatomotor
2: DAN
3: VAN
4: limbic
5: cognitive (frontoparietal)
6: DMN
```

```
0               3               0.386           0.176
1               22              0.335           0.155
2               18              0.294           0.123
3               1               0.441           0.111
4               39              0.394           0.078
5               27              0.281           0.153
6               10              0.313           0.229
```

The maps are not great, but at least it's data-driven. Another idea is to use
the actual network maps for the dual regression:

```bash
cd ~/data/heritability_change/xcp-36p_despike/yeo_masks
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
mkdir dual
while read maskid; do
    m=`printf %04d $maskid`;
    echo yeo_masks $m;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d ../Yeo_nets.nii \
        -o dual/dr_stage1_${m}.txt --demean -m ../group_epi_mask_fancy.nii;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d dual/dr_stage1_${m}.txt \
        -o dual/dr_stage2_${m} --demean -m ../group_epi_mask_fancy.nii --des_norm \
        --out_z=dual/dr_stage2_${m}_Z;
done < ../ids_1.txt
```

Time to dump to R:

```bash
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_fancy.ica/
mkdir dumps
for m in `cat ../ids_1.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in 3 22 18 1 39 27 10; do
        3dmaskdump -mask ../group_epi_mask_fancy.nii \
            -o dumps/${maskid}_IC${i}_Z.txt dual/dr_stage2_${maskid}_Z.nii.gz[${i}];
    done;
done
```

```bash
cd ~/data/heritability_change/xcp-36p_despike/yeo_masks/
mkdir dumps
for m in `cat ../ids_1.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in {0..6}; do
        3dmaskdump -mask ../group_epi_mask_fancy.nii \
            -o dumps/${maskid}_net${i}_Z.txt dual/dr_stage2_${maskid}_Z.nii.gz[${i}];
    done;
done
```

Then, we collect our results in R:

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=231015
for (m in c(3, 22, 18, 1, 39, 27, 10)) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/groupmelodic_fancy.ica/dumps/%04d_IC%d_Z.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/xcp-36p_despike/melodic_fancy_IC%d.RData.gz', m)
     save(brain_data, file=fname, compress=T)
}
```

Then, repeat the same for the yeo_masks dual regression.

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_1.txt')[, 1]
nvox=231015
for (m in 0:6) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks/dumps/%04d_net%d_Z.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks_fancy_net%d.RData.gz', m)
     save(brain_data, file=fname, compress=T)
}
```

The downside of not using MELODIC, but the Yeo masks instead is that hopefully
ICA would also get rid of movement. So, let's see if doing ICA in the FD.25
group gives better results:

```bash
mydir=/Volumes/Labs/rsfmri_36p/xcpengine_output_fc-36p_despike/
cd ~/data/heritability_change/xcp-36p_despike
fname=rsfmri_fc-36p_despike_condensed_posOnly_FD0.25_scans292_08022019.csv;
awk '{FS=","; print $1}' ../$fname | tail -n +2 > ids_p25.txt;
for maskid in `cat ids_p25.txt`; do
    m=`printf %04d $maskid`;
    echo $mydir/sub-${m}/norm/sub-${m}_std.nii.gz >> fdp25_epi.txt;
done

melodic -i fdp25_epi.txt -o groupmelodic_fancyp25.ica -v --nobet \
    -m group_epi_mask_fancy.nii --tr=2.5 --report --Oall -a concat;

cd groupmelodic_fancyp25.ica
3dMatch -inset melodic_IC.nii.gz -refset ../Yeo_nets.nii \
    -mask ../group_epi_mask_fancy.nii -prefix matches -overwrite
cat matches_REF_coeff.vals
```

Keep in mind that the code is:

```
0: visual
1: somatomotor
2: DAN
3: VAN
4: limbic
5: cognitive (frontoparietal)
6: DMN
```

```
0               2               0.375           0.176
1               27              0.361           0.155
2               10              0.428           0.123
3               4               0.447           0.111
4               31              0.415           0.078
5               29              0.305           0.153
6               7               0.328           0.229
```

Yeah, it does look better in the p25 sample. So, maybe we should run that as
well.

```bash
cd ~/data/heritability_change/xcp-36p_despike/groupmelodic_fancyp25.ica
mkdir dual
while read maskid; do
    m=`printf %04d $maskid`;
    echo ${pipe} $m;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d melodic_IC \
        -o dual/dr_stage1_${m}.txt --demean -m ../group_epi_mask_${pipe}.nii;
    $FSLDIR/bin/fsl_glm -i $mydir/sub-${m}/norm/sub-${m}_std.nii.gz -d dual/dr_stage1_${m}.txt \
        -o dual/dr_stage2_${m} --demean -m ../group_epi_mask_${pipe}.nii --des_norm
        --out_z=dual/dr_stage2_${m}_Z;
done < ../ids_p25.txt

mkdir dumps
for m in `cat ../ids_p25.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in 2 27 10 4 31 29 7; do
        3dmaskdump -mask ../group_epi_mask_fancy.nii \
            -o dumps/${maskid}_IC${i}_Z.txt dual/dr_stage2_${maskid}_Z.nii.gz[${i}];
    done;
done
```

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_p25.txt')[, 1]
nvox=231015
for (m in c(2, 27, 10, 4, 31, 29, 7)) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/groupmelodic_fancyp25.ica/dumps/%04d_IC%d_Z.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/xcp-36p_despike/melodic_fancyp25_IC%d.RData.gz', m)
     save(brain_data, file=fname, compress=T)
}
```

Now it's just a matter of assigning MRNs, calculate slopes, and
prepare it for SOLAR voxelwise.

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv')
mydir='~/data/heritability_change/xcp-36p_despike/'
for (ic in c(3, 22, 18, 1, 39, 27, 10)) {
    fname = sprintf('%s/melodic_fancy_IC%d.RData.gz', mydir, ic)
    load(fname)
    b = brain_data
    var_names = colnames(b)[2:ncol(b)]
    df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id', all.x=F)

    # make sure we still have two scans for everyone
    rm_subjs = names(which(table(df2$Medical.Record...MRN)<2))
    rm_me = df2$Medical.Record...MRN %in% rm_subjs
    df2 = df2[!rm_me, ]

    mres = df2
    mres$SX_HI = as.numeric(as.character(mres$SX_hi))
    mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))

    res = c()
    for (s in unique(mres$Medical.Record...MRN)) {
        idx = which(mres$Medical.Record...MRN == s)
        row = c(s, unique(mres[idx, 'Sex']))
        y = mres[idx[2], var_names] - mres[idx[1], var_names]
        x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
        slopes = y / x
        row = c(row, slopes)
        for (t in c('SX_inatt', 'SX_HI', 'qc')) {
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
        print(nrow(res))
    }
    colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI', 'qc',
                                                'inatt_baseline',
                                                'HI_baseline', 'DX', 'DX2'))
    # we only open this in R, so it's OK to be RData to load faster
    fname = sprintf('%s/melodic_fancy_slopes_IC%d.rds', mydir, ic)
    saveRDS(res, file=fname)

    # and remove outliers
    res_clean = res
    for (t in var_names) {
        mydata = as.numeric(res_clean[, t])
        # identifying outliers
        ul = mean(mydata) + 3 * sd(mydata)
        ll = mean(mydata) - 3 * sd(mydata)
        bad_subjs = c(which(mydata < ll), which(mydata > ul))

        # remove within-variable outliers
        res_clean[bad_subjs, t] = NA
    }
    fname = sprintf('%s/melodic_fancy_slopesClean_IC%d.rds', mydir, ic)
    saveRDS(res_clean, file=fname)

    # and make sure every family has at least two people
    good_nuclear = names(table(df2$Nuclear.ID...FamilyIDs))[table(df2$Nuclear.ID...FamilyIDs) >= 4]
    good_extended = names(table(df2$Extended.ID...FamilyIDs))[table(df2$Extended.ID...FamilyIDs) >= 4]
    keep_me = c()
    for (f in good_nuclear) {
        keep_me = c(keep_me, df2[which(df2$Nuclear.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    for (f in good_extended) {
        keep_me = c(keep_me, df2[which(df2$Extended.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    keep_me = unique(keep_me)

    fam_subjs = c()
    for (s in keep_me) {
        fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
    }
    res2 = res[fam_subjs, ]
    res2_clean = res_clean[fam_subjs, ]

    fname = sprintf('%s/melodic_fancy_slopesFam_IC%d.csv', mydir, ic)
    write.csv(res2, file=fname, row.names=F, na='', quote=F)
    fname = sprintf('%s/melodic_fancy_slopesCleanFam_IC%d.csv', mydir, ic)
    write.csv(res2_clean, file=fname, row.names=F, na='', quote=F)
}
```

And for the Yeo masks:

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv')
mydir='~/data/heritability_change/xcp-36p_despike/'
for (ic in 0:6) {
    fname = sprintf('%s/yeo_masks_fancy_net%d.RData.gz', mydir, ic)
    load(fname)
    b = brain_data
    var_names = colnames(b)[2:ncol(b)]
    df2 = merge(df, b, by.x='Mask.ID', by.y='mask.id', all.x=F)

    # make sure we still have two scans for everyone
    rm_subjs = names(which(table(df2$Medical.Record...MRN)<2))
    rm_me = df2$Medical.Record...MRN %in% rm_subjs
    df2 = df2[!rm_me, ]

    mres = df2
    mres$SX_HI = as.numeric(as.character(mres$SX_hi))
    mres$SX_inatt = as.numeric(as.character(mres$SX_inatt))

    res = c()
    for (s in unique(mres$Medical.Record...MRN)) {
        idx = which(mres$Medical.Record...MRN == s)
        row = c(s, unique(mres[idx, 'Sex']))
        y = mres[idx[2], var_names] - mres[idx[1], var_names]
        x = mres[idx[2], 'age_at_scan'] - mres[idx[1], 'age_at_scan']
        slopes = y / x
        row = c(row, slopes)
        for (t in c('SX_inatt', 'SX_HI', 'qc')) {
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
        print(nrow(res))
    }
    colnames(res) = c('ID', 'sex', var_names, c('SX_inatt', 'SX_HI', 'qc',
                                                'inatt_baseline',
                                                'HI_baseline', 'DX', 'DX2'))
    # we only open this in R, so it's OK to be RData to load faster
    fname = sprintf('%s/yeo_masks_fancy_slopes_net%d.rds', mydir, ic)
    saveRDS(res, file=fname)

    # and remove outliers
    res_clean = res
    for (t in var_names) {
        mydata = as.numeric(res_clean[, t])
        # identifying outliers
        ul = mean(mydata) + 3 * sd(mydata)
        ll = mean(mydata) - 3 * sd(mydata)
        bad_subjs = c(which(mydata < ll), which(mydata > ul))

        # remove within-variable outliers
        res_clean[bad_subjs, t] = NA
    }
    fname = sprintf('%s/yeo_masks_fancy_slopesClean_net%d.rds', mydir, ic)
    saveRDS(res_clean, file=fname)

    # and make sure every family has at least two people
    good_nuclear = names(table(df2$Nuclear.ID...FamilyIDs))[table(df2$Nuclear.ID...FamilyIDs) >= 4]
    good_extended = names(table(df2$Extended.ID...FamilyIDs))[table(df2$Extended.ID...FamilyIDs) >= 4]
    keep_me = c()
    for (f in good_nuclear) {
        keep_me = c(keep_me, df2[which(df2$Nuclear.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    for (f in good_extended) {
        keep_me = c(keep_me, df2[which(df2$Extended.ID...FamilyIDs == f),
                                'Medical.Record...MRN'])
    }
    keep_me = unique(keep_me)

    fam_subjs = c()
    for (s in keep_me) {
        fam_subjs = c(fam_subjs, which(res[, 'ID'] == s))
    }
    res2 = res[fam_subjs, ]
    res2_clean = res_clean[fam_subjs, ]

    fname = sprintf('%s/yeo_masks_fancy_slopesFam_net%d.csv', mydir, ic)
    write.csv(res2, file=fname, row.names=F, na='', quote=F)
    fname = sprintf('%s/yeo_masks_fancy_slopesCleanFam_net%d.csv', mydir, ic)
    write.csv(res2_clean, file=fname, row.names=F, na='', quote=F)
}
```

**Note that I won' be removing movement here! Mostly because since we're using
ICA, the movement components should have been isolated already. But we can
always check any results later against correlation to movement.
**

# 2019-08-06 09:46:43

As the CSV files are getting ready, let's dig out the commands to run the
voxelwise SOLAR:

```bash
cd ~/data/heritability_change/xcp-36p_despike
# create master list of voxels
nvox=231015;
for i in `seq 1 $nvox`; do
    echo $i >> voxel_list.txt;
done
# our previous experiments had chunks of 5K voxels
split -da 2 -l 5000 voxel_list.txt vlist --additional-suffix=".txt";
```

Now we set up the swarm:

```bash
cd ~/data/heritability_change/xcp-36p_despike;
phen_file=melodic_fancyp25_slopesCleanFam_IC2;
jname=fancyp25_2c;
swarm_file=swarm.${jname};

rm -f $swarm_file;
for vlist in `ls $PWD/vlist*txt`; do  # getting full path to files
    echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
done;
swarm --gres=lscratch:100 -f $swarm_file --module solar -t 32 -g 20 \
        --logdir=trash_${jname} --job-name ${jname} --time=2:00:00 --merge-output \
        --partition quick,norm
```

For fancyp25 I ran 1000 voxels in a 32core node in 13 min. So, if I'm doing
bundles of 5000, I think I'll be safe to do allocate 2h for each. And even with
32 cores it's only taking 2.5Gb of RAM.

The swarms actually completed in about 1h, so 2h is quite conservative indeed. I
wonder if we can go to 16 cores (to get machines more quickly)? And cap it to 4h
just to stay conservative and still fit into the quick partititon, even though I
think we can probably do it in 3h... keep in mind this is the p25 set, so maybe
the bigger set will take longer.

It's taking about 2.5h in a 16 core node, which is pretty good and should be ok
for the FD1 sample as well?

# 2019-08-07 08:40:50

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in 2 27 10 4 31 29 7; do
    phen_file=melodic_fancyp25_slopesCleanFam_IC${i};
    jname=fancyp25_${i}c;
    swarm_file=swarm.${jname};

    rm -f $swarm_file;
    for vlist in `ls $PWD/vlist*txt`; do  # getting full path to files
        echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
    done;
    swarm --gres=lscratch:10 -f $swarm_file --module solar -t 16 -g 5 \
            --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
            --partition quick,norm
done
```

To compile, we do something like this:

```bash
#desktop
cd ~/data/heritability_change/xcp-36p_despike/
cut -d " " -f 1,2,3 \
    groupmelodic_fancy.ica/dumps/1351_IC3_Z.txt > group_mask_fancy_ijk.txt
```

```bash
module load afni

cd /lscratch/${SLURM_JOBID}
for i in 2 27 10 4 31 29 7; do
    phen=melodic_fancyp25_slopesFam_IC${i};
    mkdir $phen;
    cd $phen;
    cp ~/data/tmp/$phen/*.out ~/data/tmp/${phen}/*gz .;
    for f in `/bin/ls *gz`; do tar -zxf $f; done
    cd ..
    python ~/research_code/fmri/compile_solar_voxel_results.py \
        /lscratch/${SLURM_JOBID}/ $phen \
        ~/data/heritability_change/xcp-36p_despike/group_epi_mask_fancy.nii;
done
```

For FD1, we do:

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in 18 22 3 1 39 27 10; do
    phen_file=melodic_fancy_slopesCleanFam_IC${i};
    jname=fancy_${i}c;
    swarm_file=swarm.${jname};

    rm -f $swarm_file;
    for vlist in `ls $PWD/vlist*txt`; do  # getting full path to files
        echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
    done;
    swarm --gres=lscratch:10 -f $swarm_file --module solar -t 16 -g 5 \
            --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
            --partition quick,norm
done
```

And Yeo masks:

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in {0..6}; do
    phen_file=yeo_masks_fancy_slopesCleanFam_net${i};
    jname=yeo_${i}c;
    swarm_file=swarm.${jname};

    rm -f $swarm_file;
    for vlist in `ls $PWD/vlist*txt`; do  # getting full path to files
        echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
    done;
    swarm --gres=lscratch:10 -f $swarm_file --module solar -t 16 -g 5 \
            --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
            --partition quick,norm
done
```

As usual, do it for clean and nonClean, just for kicks...

# 2019-08-08 10:03:53

Good to know: some of the FD1 sets are not finishing... might need to increase
the time? Or run in 32 cores? Even the ones that completed took 3h:45min, so
that's cutting close. We either do more time and forgo the quick partition, or
increase it to 32 cores...

Let's check some of the results we already have:

```bash
#desktop
cd ~/data/heritability_change/xcp-36p_despike
rm fancyp25_clean_clusters.txt;
for f in `/bin/ls polygen_results_melodic_fancyp25_slopesCleanFam_IC*.nii`; do
    3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 \
        -NN1 30 $f >> fancyp25_clean_clusters.txt;
done
rm fancyp25_clusters.txt;
for f in `/bin/ls polygen_results_melodic_fancyp25_slopesFam_IC*.nii`; do
    3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 \
        -NN1 30 $f >> fancyp25_clusters.txt;
done
```

The unClean results look stronger... but now it's the chicken/egg situation. Run
permutations first or look at results? Since we have the cluster empty, let's
run some permutations bundled so we don't crush it. But first, let's generate
some permutations, for the more interesting networks first:

```r
# start it from lscratch
# 2 27 10 4 31 29 7
m = 7
nperms = 100
library(data.table)
dread = fread(sprintf('~/data/heritability_change/melodic_fancyp25_slopesFam_IC%d.csv', m), header = T, sep = ',')
d = as.data.frame(dread)  # just so we can index them a bit easier
vcols = c(which(grepl("v",colnames(d))), which(grepl("sex",colnames(d))))
d2 = d
for (p in seq(1, nperms, 2)) {
    d2[, vcols] = d[sample(nrow(d)), vcols]
    fname = sprintf('~/data/heritability_change/melodic_fancyp25_slopesFam_IC%d_p%03d.csv', m, p)
    print(fname)
    fwrite(d2, file=fname, row.names=F, quote=F)
}
```

<!-- while [ $cur_vox -lt $nvox ]; do
    let last_vox=${cur_vox}+${bundle}-1;
    # gets the min
    last_vox=$(($last_vox<$nvox?$last_vox:$nvox))
    echo "bash ~/research_code/run_solar_voxel_range.sh melodic_fancy_slopesClean_n111_IC${ic}_p0000_sexAndBrainShuffled ${cur_vox} ${last_vox}" >> $fname;
    let cur_vox=${last_vox}+1;
done;
# just copy it and rename IC and perm for the other ones
for n in `seq 1 $nperms`; do
    perm=`printf %04d $n`;
    cp $fname `printf net%d_p%04d $ic $n`.swarm;
    sed -i -- "s/p0000/p${perm}/g" `printf net%d_p%04d $ic $n`.swarm;
done

cd ~/data/heritability_change/fmri_same_space/
ic=5;
jstart=124;  # permutation to start with
jdeps=4;  # number of dependent jobs in each swarm
nswarms=16;  # number of independent swarms
for n in `seq 1 $nswarms`; do
    jname=`printf net%d_p%04d $ic $jstart`;
    cur_id=$(swarm --gres=lscratch:1 -f ${jname}.swarm --module solar -g 1 -t 1 \
        --logdir=${jname} --job-name ${jname} -p 2 --time=36:00:00 --merge-output);
    echo "Active swarm: ${jname} (${cur_id})"
    for d in `seq 1 $jdeps`; do
        let jstart=${jstart}+1;
        jname=`printf net%d_p%04d $ic $jstart`;
        job_id=$(swarm --gres=lscratch:1 -f ${jname}.swarm --module solar -g 1 -t 1 \
            --logdir=${jname} --job-name ${jname} -p 2 \
            --time=36:00:00 --merge-output --dependency afterany:$cur_id);
        echo "Dependent swarm: ${jname} (${job_id}, on ${cur_id})"
        cur_id=$job_id;
    done;
    let jstart=${jstart}+1;
done-->

While that's being generated, let's see if the current clusters look nice. Going
backwards to start with DMN:

```bash
#desktop
3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -orient LPI \
    -savemask mymask.nii -NN1 60 \
    polygen_results_melodic_fancyp25_slopesFam_IC7.nii
```

IC7 (DMN) doesn't look that great, mostly sensory cortex? It does have 2
clusters above 60 (whatever that means, as we don't have significance thresholds
yet), and the second one is in IFG.

![](images/2019-08-08-11-49-17.png)
![](images/2019-08-08-11-50-27.png)

The only cluster for cognitive control was in the inferior temporal gyrus, and
it's only 50 voxels:

![](images/2019-08-08-11-52-50.png)

For affective we had a 61 cluster (postcentral), and a 58 cluster (postcentral too):

![](images/2019-08-08-12-03-36.png)
![](images/2019-08-08-12-04-19.png)

For VAN (4), we got 3 clusters above 60: SFG, lingual, and the last one looks to
be in a ventricle... never good.

![](images/2019-08-08-12-13-05.png)
![](images/2019-08-08-12-13-46.png)
![](images/2019-08-08-12-14-50.png)

Finally, for DAN (10) only saw a single cluster above 50, and it's falling in
the white matter too... not good:

![](images/2019-08-08-12-17-17.png)

I'm not too excited about these results, so I stopped the p25 permutations.
Let's see if the FD1 or ye-Masks results are better.

# 2019-08-09 09:54:47

Let's then investigate the Yeo mask results and the MELODIC FD1 to see if
anything looks interesting. 

```bash
#desktop
cd ~/data/heritability_change/xcp-36p_despike
rm ym_clean_clusters.txt;
for f in `/bin/ls polygen_results_yeo_masks_fancy_slopesCleanFam_net*.nii`; do
    3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 \
        -NN1 100 $f >> ym_clean_clusters.txt;
done
rm ym_clusters.txt;
for f in `/bin/ls polygen_results_yeo_masks_fancy_slopesFam_net*.nii`; do
    3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 \
        -NN1 100 $f >> ym_clusters.txt;
done
```

I got some huge clusters, especially compared to what I was seeing for p25.
Maybe it'll be different for FD1... 