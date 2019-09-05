# 2019-09-05 11:30:09

I was curious about whether using a stricter voxelwise p-value would help the
DTI results a bit. Let's see what happens:

```bash
# bw
cd ~/data/heritability_change
3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.99 -NN3 2 polygen_results_dti_rd_residNoSex_OLS_naSlopes133Clean.nii
```

(NN1; NN2; NN3)
FA: 3, 4, 4
AD: 4, 5, 5
RD: 4, 5, 5

OK, these clusters are quite small. But how big do they get when doing
permutations?

```bash
cd ~/data/heritability_change/perms
froot=polygen_results_dti_ad_residNoSex_OLS_naSlopes133Clean;
csize=4;
res=`3dclust -1Dformat -nosum -1dindex 0 -1tindex 1 -1thresh 0.99 -NN1 $csize \
    -quiet ${froot}_p*.nii | grep CLUSTERS | wc -l`
nperms=`ls -1 ${froot}_p*.nii | wc -l`;
p=$(bc <<<"scale=3;($nperms - $res)/$nperms")
echo negatives=${res}, perms=${nperms}, pval=$p
```

The AD cluster was .056 using the 500 perms in NN2, and .068 in NN3. The RD cluster
was .164 and .218, respectively. So, nothing very exciting. If we try NN1, we
get .1 for AD and for .192 RD. So, nothing worth coming down from p = .05...