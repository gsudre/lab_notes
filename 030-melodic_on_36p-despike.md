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

Now, it's time to figure out which ICs we are going to use:

We are already in MNI space, like the Yeo networks. So, let's go with that:






```bash
#bw
cd ~/data/heritability_change/fmri_same_space/epi/;
for i in {1..7}; do
    3dcalc -prefix Yeo_liberal_inMNI_net${i}.nii \
        -a /data/NCR_SBRB/software/Yeo_JNeurophysiol11_MNI152/Yeo2011_7Networks_MNI152_FreeSurferConformed1mm_LiberalMask.nii.gz -expr "amongst(a,${i})";
done
3dTcat -prefix Yeo_liberal_inMNI_combined.nii Yeo_liberal_inMNI_net1.nii \
    Yeo_liberal_inMNI_net2.nii Yeo_liberal_inMNI_net3.nii \
    Yeo_liberal_inMNI_net4.nii Yeo_liberal_inMNI_net5.nii \
    Yeo_liberal_inMNI_net6.nii Yeo_liberal_inMNI_net7.nii
3dresample -master groupmelodic_inter.ica/melodic_IC.nii.gz \
    -prefix Yeo_nets.nii -inset Yeo_liberal_inMNI_combined.nii \
    -rmode NN -overwrite
```

So, let's figure out what are the best matching networks for each mask:

```bash
#bw
cd ~/data/heritability_change/fmri_same_space/epi/groupmelodic_inter.ica/
3dMatch -inset melodic_IC.nii.gz -refset ../Yeo_nets.nii \
    -mask ../../group_epi_mask_inter.nii -prefix matches -overwrite
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

inter:

```
0               41              0.173           0.251
1               48              0.240           0.144
2               35              0.325           0.083
3               23              0.214           0.170
4               214             0.261           0.025
5               228             0.174           0.191
6               0               0.179           0.296
```

I didn't like the DMN component in inter.

fancy:

```
0               60              0.364           0.138
1               15              0.344           0.122
2               2               0.381           0.097
3               8               0.365           0.086
4               62              0.413           0.074
5               5               0.350           0.122
6               0               0.322           0.181
```

The one from fancy looks much better:

![](images/2019-05-16-11-31-29.png)

And this is cognitive:

![](images/2019-05-16-11-43-11.png)

Let's go with fancy for now.

union:

```
0               107             0.340           0.054
1               138             0.295           0.048
2               70              0.310           0.038
3               34              0.310           0.033
4               67              0.405           0.029
5               98              0.256           0.048
6               48              0.267           0.073
```

Union mask has too much crap...

```bash
pipe=fancy;
cd ~/data/heritability_change/fmri_same_space/epi/groupmelodic_${pipe}.ica/dual
mkdir dumps
for m in `cat ../../../../3min_mni.txt`; do
    maskid=`printf %04d $m`;
    echo $maskid;
    rm dumps/${maskid}_*.txt
    for i in 60 15 2 8 62 5 0; do  # fancy
        3dmaskdump -mask ../../../group_epi_mask_${pipe}.nii \
            -o dumps/${maskid}_IC${i}_Z.txt dr_stage2_${maskid}_Z.nii.gz[${i}];
    done;
done
```

Then, we collect our results in R:

```r
maskids = read.table('~/data/heritability_change/3min_mni.txt')[, 1]
nvox=154058
for (m in c(60, 15, 2, 8, 62, 5, 0)) {
    print(m)
    brain_data = matrix(nrow=length(maskids), ncol=nvox)
    for (s in 1:nrow(brain_data)) {
        fname = sprintf('~/data/heritability_change/fmri_same_space/epi/groupmelodic_fancy.ica/dual/dumps/%04d_IC%d_Z.txt', maskids[s], m)
        a = read.table(fname)
        brain_data[s, ] = a[,4]
     }
     brain_data = cbind(maskids, brain_data)
     cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
     colnames(brain_data) = cnames
     fname = sprintf('~/data/heritability_change/fmri_same_space/melodic_fancy_IC%d.RData', m)
     save(brain_data, file=fname)
}
```

Now that the data is into CSV, we need to assign MRNs, calculate slopes, and
prepare it for SOLAR voxelwise.

**Note that I won' be removing movement here! Mostly because since we're using
ICA, the movement components should have been isolated already. But we can
always check any results later against correlation to movement.
**
