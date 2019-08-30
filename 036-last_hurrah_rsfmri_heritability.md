# 2019-08-30 08:20:52

After the disappointing results I got yesterday with the p<.01 threshold, I
wondered whether I could push them over the edge by considering everyone in the
sample (instead of only kids with sibs) when calculating heritability. Or,
another option would be to use a more stringent threshold (more than FD=1) to
get cleaner data for SOLAR. Let's try those options.

But I'll only run them for the results that were close:

yeo
'6': 60 (p=.081, n=98)
'2': 57 (p=.1, n=100)
'4': 59 (p=.05, n=100)

My last resort will actually be to try a 400x7 connectivity matrix (or less than
400... we'll see).

```r
mydir='~/data/heritability_change/xcp-36p_despike/'
suf = ''
for (ic in c(2, 4, 6)) {
    for (suf in c('', '_Z')) {
        print(ic)
        print(suf)
        fname = sprintf('%s/yeo_masks_gray_slopes_net%d%s.rds', mydir, ic, suf)
        res = readRDS(fname)
        fname = sprintf('%s/yeo_masks_gray_slopes_net%d%s.csv', mydir, ic, suf)
        write.csv(res, file=fname, row.names=F, na='', quote=F)
    }
}
```

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in 2 4 6; do
    for suf in '' '_Z'; do
        phen_file=yeo_masks_gray_slopes_net${i}${suf};
        jname=ymAll_${i}${suf};
        swarm_file=swarm.${jname};

        rm -f $swarm_file;
        for vlist in `ls $PWD/vlistg*txt`; do  # getting full path to files
            echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
        done;
        swarm --gres=lscratch:10 -f $swarm_file --module solar -t 32 -g 10 \
                --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 \
                --merge-output --partition quick,norm
    done;
done
```

Then, while that's running, let's set up a few situations with cleaner data.
Because this is not MELODIC, we don't need to recalculate masks.

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

Because we have already dumped everyone with FD < 1, we just need to collect the
best ones in R:

```r
maskids = read.table('~/data/heritability_change/xcp-36p_despike/ids_p25.txt')[, 1]
nvox=155301
for (m in 0:6) {
    for (suf in c('', '_Z')) {
        print(m)
        print(suf)
        brain_data = matrix(nrow=length(maskids), ncol=nvox)
        for (s in 1:nrow(brain_data)) {
            fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks_gray/dumps/%04d_net%d%s.txt', maskids[s], m, suf)
            a = read.table(fname)
            brain_data[s, ] = a[,4]
        }
        brain_data = cbind(maskids, brain_data)
        cnames = c('mask.id', sapply(1:nvox, function(d) sprintf('v%06d', d)))
        colnames(brain_data) = cnames
        fname = sprintf('~/data/heritability_change/xcp-36p_despike/yeo_masks_grayp25_net%d%s.rds', m, suf)
        saveRDS(brain_data, file=fname)
    }
}
```

Finally, make the slopes:

```r
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
df = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_condensed_posOnly_FD1.00_scans520_08022019.csv')
mydir='~/data/heritability_change/xcp-36p_despike/'
ic = 5
suf = ''

fname = sprintf('%s/yeo_masks_grayp25_net%d%s.rds', mydir, ic, suf)
b = readRDS(fname)
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
fname = sprintf('%s/yeo_masks_grayp25_slopes_net%d%s.rds', mydir, ic, suf)
saveRDS(res, file=fname)
# in case we want to run everyone, not just family through solar
fname = sprintf('%s/yeo_masks_grayp25_slopes_net%d%s.csv', mydir, ic, suf)
write.csv(res, file=fname, row.names=F, na='', quote=F)

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

fname = sprintf('%s/yeo_masks_grayp25_slopesFam_net%d%s.csv', mydir, ic, suf)
write.csv(res2, file=fname, row.names=F, na='', quote=F)
```

And as usual, we set up the swarms:

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in {0..6}; do
    for suf in '' '_Z'; do
        for f in '' 'Fam'; do
            phen_file=yeo_masks_grayp25_slopes${f}_net${i}${suf};
            jname=ym${f}p25_${i}${suf};
            swarm_file=swarm.${jname};

            rm -f $swarm_file;
            for vlist in `ls $PWD/vlistg*txt`; do  # getting full path to files
                echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
            done;
            echo "ERROR" > swarm_wait;
            while grep -q ERROR swarm_wait; do
                echo "Trying $jname"
                swarm --gres=lscratch:10 -f $swarm_file --module solar -t 32 -g 10 \
                        --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
                        --partition quick,norm 2> swarm_wait;
                if grep -q ERROR swarm_wait; then
                    echo -e "\tError, sleeping..."
                    sleep 30m;
                fi;
            done;
        done;
    done;
done
```

And while we wait for SOLAR to run, let's go ahead and generate the perms, for
the non-Fam cases, and then for p25:

```r
suf = '0'
m = 'yeo_masks'
start=1
nperms = 25
step=1

library(data.table)
set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
dread = fread(sprintf('~/data/heritability_change/%s_grayp25_slopes_net%s.csv', m, suf),
              header = T, sep = ',')
d = as.data.frame(dread)  # just so we can index them a bit easier
vcols = c(which(grepl("v",colnames(d))), which(grepl("sex",colnames(d))),
          which(grepl("qc",colnames(d))))
d2 = d
for (p in seq(start, nperms, step)) {
    d2[, vcols] = d[sample(nrow(d)), vcols]
    fname = sprintf('~/data/heritability_change/%s_grayp25_slopes_net%s_p%03d.csv', m, suf, p)
    print(fname)
    fwrite(d2, file=fname, row.names=F, quote=F)
}
```

OK, let's see if the clusters are any bigger using this bigger sample:

```bash
module load afni

cd /lscratch/${SLURM_JOBID}
for i in 2 4 6; do
    for suf in '' '_Z'; do
        phen=yeo_masks_gray_slopes_net${i}${suf};
        mkdir $phen;
        cd $phen;
        cp ~/data/tmp/${phen}/*gz .;
        for f in `/bin/ls *gz`; do tar -zxf $f; done
        cd ..
        python ~/research_code/fmri/compile_solar_voxel_results.py \
            /lscratch/${SLURM_JOBID}/ $phen \
            ~/data/heritability_change/xcp-36p_despike/gray_matter_mask.nii;
        rm -rf $phen;
    done;
done
cp polygen*yeo_masks_gray_slopes_net*nii ~/data/heritability_change/xcp-36p_despike/

cd ~/data/heritability_change/xcp-36p_despike/
for i in 2 4 6; do
    for suf in '' '_Z'; do
        phen=yeo_masks_gray_slopes_net${i}${suf};
        3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.99 -orient LPI \
            -savemask ${phen}_grayAll_NN1_clusters.nii -NN1 40 \
            polygen_results_${phen}.nii >> NN1_yeo_masks_grayAll_results.txt;
    done
done
```

There wasn't again a clear difference between Z and non-Z, but the size of the
clusters wasn't much different, as expected based on how SOLAR works:

yeo
'6': 56 (before 60)
'2': 61 (before 57)
'4': 50 (before 59)

And I don't think the perms will make much difference here, so I won't even run
them. Let's put our eggs in the p25 basket for now.


<!--

# 2019-08-16 11:30:03

While the perms for melodic are running, let's take a look at the yeo mask
results. Are they worth checking on through the small permutation approach? Are
Z results better than regular?

```bash
module load afni

cd /lscratch/${SLURM_JOBID}
for i in {0..6}; do
    for suf in '' '_Z'; do
        phen=yeo_masks_gray_slopesFam_net${i}${suf};
        mkdir $phen;
        cd $phen;
        cp ~/data/tmp/${phen}/*gz .;
        for f in `/bin/ls *gz`; do tar -zxf $f; done
        cd ..
        python ~/research_code/fmri/compile_solar_voxel_results.py \
            /lscratch/${SLURM_JOBID}/ $phen \
            ~/data/heritability_change/xcp-36p_despike/gray_matter_mask.nii;
        rm -rf $phen;
    done;
done
cp polygen*yeo_masks_gray_slopesFam_net*nii ~/data/heritability_change/xcp-36p_despike/

cd ~/data/heritability_change/xcp-36p_despike/
for i in {0..6}; do
    for suf in '' '_Z'; do
        phen=yeo_masks_gray_slopesFam_net${i}${suf};
        3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -orient LPI \
            -savemask ${phen}_NN1_clusters.nii -NN1 50 \
            polygen_results_${phen}.nii >> NN1_yeo_masks_gray_results.txt;
    done
done
```

The difference wasn't huge, for for most of the ICs that matter the best result
was with the nonZ results. Let's make a few figures to see if they are worth
exploring, as we cannot submit much new stuff in the queue at this moment
anyways. I'm centering my figures for now on max intensity, because the clusters
are looking better this way, compared to using COM.

cognitive:
![](images/2019-08-16-12-44-03.png)
That's a very nice precuneus hit, which is normally associated with DMN. That
seems to be the second hit for _Z. There, the first hit is quite frontal.

DMN: 
got a hit on frontal pole... could be interesting if it survives
permutations
![](images/2019-08-16-12-39-34.png)
Note that the best cluster for DMN_Z is in mFG/iFG, so that could potentially be
more interesting...
![](images/2019-08-16-12-42-00.png)

But let's not get too excited about these results. It'll be cool if they
survive, but let's run some perms before we do much else:


```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in '6' '6_Z' '5' '5_Z' '2' '3' '4'; do
    for p in {1..25}; do
        perm=`printf %03d $p`;
        phen_file=yeo_masks_gray_slopesFam_net${i}_p${perm};
        swarm_file=swarm.ymg${i}_p${perm};

        for vlist in `ls $PWD/vlistg*txt`; do  # getting full path to files
            echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
        done;
    done;
done

for i in '6' '6_Z' '5' '5_Z' '2' '3' '4'; do
    for p in {1..25}; do
        perm=`printf %03d $p`;
        jname=ymg${i}_p${perm};
        swarm_file=swarm.${jname};
        echo "ERROR" > swarm_wait;
        while grep -q ERROR swarm_wait; do
            echo "Trying $jname"
            swarm --gres=lscratch:10 -f $swarm_file --module solar -t 32 -g 10 \
                    --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
                    --partition quick,norm 2> swarm_wait;
            if grep -q ERROR swarm_wait; then
                echo -e "\tError, sleeping..."
                sleep 30m;
            fi;
        done;
    done;
done
```

# 2019-08-19 10:01:45

Let's compile all permutations we've been waiting on:

```bash
module load afni

cd /lscratch/${SLURM_JOBID}
for i in '6' '6_Z' '5' '5_Z' '2' '3' '4'; do
    for p in {1..25}; do
        perm=`printf %03d $p`;
        phen=yeo_masks_gray_slopesFam_net${i}_p${perm};
        mkdir $phen;
        cd $phen;
        cp ~/data/tmp/${phen}/*gz .;
        for f in `/bin/ls *gz`; do tar -zxf $f; done
        cd ..
        python ~/research_code/fmri/compile_solar_voxel_results.py \
            /lscratch/${SLURM_JOBID}/ $phen \
            ~/data/heritability_change/xcp-36p_despike/gray_matter_mask.nii;
        rm -rf $phen;
    done;
done
```

```bash
module load afni

cd /lscratch/${SLURM_JOBID}
for i in '8' '9_Z' '27' '6_Z' '10_Z'; do
    for p in {1..25}; do
        perm=`printf %03d $p`;
        phen=melodic_gray_slopesFam_IC${i}_p${perm};
        mkdir $phen;
        cd $phen;
        cp ~/data/tmp/${phen}/*gz .;
        for f in `/bin/ls *gz`; do tar -zxf $f; done
        cd ..
        python ~/research_code/fmri/compile_solar_voxel_results.py \
            /lscratch/${SLURM_JOBID}/ $phen \
            ~/data/heritability_change/xcp-36p_despike/gray_matter_mask.nii;
        rm -rf $phen;
    done;
done
```

Because we have many different combinations, I don't want to use the loop code
to figure out optimal cluster size. Also, with only 25 permutations it won't
tell me much. I much rather just check the chances of our buggest cluster and
see what I get:

```bash
cd ~/data/heritability_change/xcp-36p_despike/perms
froot=polygen_results_yeo_masks_gray_slopesFam_net6
csize=145;
res=`3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.95 -NN1 $csize \
    -quiet ${froot}_p*.nii | grep CLUSTERS | wc -l`
nperms=`ls -1 ${froot}_p*.nii | wc -l`;
p=$(bc <<<"scale=3;($nperms - $res)/$nperms")
echo negatives=${res}, perms=${nperms}, pval=$p
```

melodic:
'8': 78 (p=.24)
'9_Z': 80 (p=.24)
'27': 77 (p=.24)
'6_Z': 72 (p=.24)
'10_Z': 53 ...

yeo
'6': 145 (p = .434)
'6_Z': 149 (p=.56)
'5': 233 (p=.08)
'5_Z': 166 (p=.28)
'2': 193 (p=.28)
'3': 114 ...
'4': 188 ...

Again, these are using only 25 perms. But they don't look too promising. Yeo 5
might work out, but it will be though. And that's just one too...

The IC results look weird. Probably an error somewhere, but not promising enough
to look into it. 

# 2019-08-21 16:53:22

Let's just check on p<.01:

```bash
cd ~/data/heritability_change/xcp-36p_despike/
for i in {0..6}; do
    for suf in '' '_Z'; do
        phen=yeo_masks_gray_slopesFam_net${i}${suf};
        3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.99 -orient LPI \
            -savemask ${phen}_NN1_clusters_p01.nii -NN1 20 \
            polygen_results_${phen}.nii >> NN1_yeo_masks_gray_results_p01.txt;
    done
done
```

```bash
cd ~/data/heritability_change/xcp-36p_despike/perms
froot=polygen_results_yeo_masks_gray_slopesFam_net4
csize=59;
res=`3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.99 -NN1 $csize \
    -quiet ${froot}_p*.nii | grep CLUSTERS | wc -l`
nperms=`ls -1 ${froot}_p*.nii | wc -l`;
p=$(bc <<<"scale=3;($nperms - $res)/$nperms")
echo negatives=${res}, perms=${nperms}, pval=$p
```

yeo
'6': 60 (p=0)
'6_Z': 46 (p=.16)
'5': 39 (p=.32)
'5_Z': 32 (p=.24)
'2': 57 (p=.12)
'3': 23 ...
'4': 59 (p=0)

Well, there might be some interesting stuff here. Let's leave some more perms
running then:

```r
m = 'yeo_masks'
suf = '6'
start = 27 # 26
nperms = 100
step = 2

library(data.table)
set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
dread = fread(sprintf('~/data/heritability_change/%s_gray_slopesFam_net%s.csv', m, suf),
              header = T, sep = ',')
d = as.data.frame(dread)  # just so we can index them a bit easier
vcols = c(which(grepl("v",colnames(d))), which(grepl("sex",colnames(d))),
          which(grepl("qc",colnames(d))))
d2 = d
for (p in seq(start, nperms, step)) {
    d2[, vcols] = d[sample(nrow(d)), vcols]
    fname = sprintf('~/data/heritability_change/%s_gray_slopesFam_net%s_p%03d.csv', m, suf, p)
    print(fname)
    fwrite(d2, file=fname, row.names=F, quote=F)
}
```

```bash
cd ~/data/heritability_change/xcp-36p_despike;
for i in '6' '4' '2' '3' '5'; do
    for p in {26..100}; do
        perm=`printf %03d $p`;
        phen_file=yeo_masks_gray_slopesFam_net${i}_p${perm};
        swarm_file=swarm.ymg${i}_p${perm};

        for vlist in `ls $PWD/vlistg*txt`; do  # getting full path to files
            echo "bash ~/research_code/run_solar_voxel_parallel.sh $phen_file $vlist" >> $swarm_file;
        done;
    done;
done

# just because I couldn't wait before starting this...
sleep 5h;
for i in '6' '4' '2' '3' '5'; do
    for p in {26..100}; do
        perm=`printf %03d $p`;
        jname=ymg${i}_p${perm};
        swarm_file=swarm.${jname};
        echo "ERROR" > swarm_wait;
        while grep -q ERROR swarm_wait; do
            echo "Trying $jname"
            swarm --gres=lscratch:10 -f $swarm_file --module solar -t 32 -g 10 \
                    --logdir=trash_${jname} --job-name ${jname} --time=4:00:00 --merge-output \
                    --partition quick,norm 2> swarm_wait;
            if grep -q ERROR swarm_wait; then
                echo -e "\tError, sleeping..."
                sleep 30m;
            fi;
        done;
    done;
done
```

# 2019-08-29 09:32:50

Let's go ahead and try compiling these new permutations:

```bash
module load afni

cd /lscratch/${SLURM_JOBID}
for i in '6' '4' '2' '3' '5'; do
    for p in {26..100}; do
        perm=`printf %03d $p`;
        phen=yeo_masks_gray_slopesFam_net${i}_p${perm};
        mkdir $phen;
        cd $phen;
        cp ~/data/tmp/${phen}/*gz .;
        for f in `/bin/ls *gz`; do tar -zxf $f; done
        cd ..
        python ~/research_code/fmri/compile_solar_voxel_results.py \
            /lscratch/${SLURM_JOBID}/ $phen \
            ~/data/heritability_change/xcp-36p_despike/gray_matter_mask.nii;
        rm -rf $phen;
    done;
done
```

```bash
cd ~/data/heritability_change/xcp-36p_despike/perms
froot=polygen_results_yeo_masks_gray_slopesFam_net4
csize=59;
res=`3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.99 -NN1 $csize \
    -quiet ${froot}_p*.nii | grep CLUSTERS | wc -l`
nperms=`ls -1 ${froot}_p*.nii | wc -l`;
p=$(bc <<<"scale=3;($nperms - $res)/$nperms")
echo negatives=${res}, perms=${nperms}, pval=$p
```

yeo
'6': 60 (p=.081, n=98)
'5': 39 (p=.32, n=25)
'2': 57 (p=.1, n=100)
'3': 23 (p=.83, n=66)
'4': 59 (p=.05, n=100)

We're hovering, but definitely not good enough. Let's try the other options... -->