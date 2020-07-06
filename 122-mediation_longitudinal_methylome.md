# 2020-06-25 16:14:02

Just some notes for running the script. I ran 3 Xs in 4 minutes using a single
thread. If I split the files into groups of 15, that should give me 20 minutes
per Rscript call, which means 480 Xs every 20 minutes if using 32 cores. Even
the HI dataset, with 3.3K Xs, should take tops 3h. 

```bash
cd ~/data/longitudinal_methylome
split -l 15 IN_xs inatt_
ls -1 inatt_* > IN_sets.txt;
cat IN_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        INATT_ROC_methyl_sx_dti_82.csv {} IN_ms IN_ys res_{}.csv;
```

And for HI is very similar:

```bash
cd ~/data/longitudinal_methylome
split -l 15 HI_xs hi_
ls -1 hi_* > HI_sets.txt;
cat HI_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        HI_ROC_methyl_sx_dti_82.csv {} HI_ms HI_ys res_{}.csv;
```

Then, to compile: 

```bash
cd ~/data/longitudinal_methylome/results
head -n 1 res_inatt_bj.csv > inatt_compiled.csv;
for f in `ls -1 res_inatt_??.csv`; do tail -n +2 $f >> inatt_compiled.csv ; done
```

And I added the FDR, within and across Ms, in R:

```r
a = read.csv('~/data/longitudinal_methylome/results/inatt_compiled.csv')
a$tot_p_FDR = p.adjust(a$tot_p, method='fdr')
a$acme_p_FDR = p.adjust(a$acme_p, method='fdr')
a$ade_p_FDR = p.adjust(a$ade_p, method='fdr')
a$prop_p_FDR = p.adjust(a$prop_p, method='fdr')
idx = a$M=='AD_left_unc'
a$tot_p_FDR_withinM = NA
a$acme_p_FDR_withinM = NA
a$ade_p_FDR_withinM = NA
a$prop_p_FDR_withinM = NA
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
idx = a$M=='RD_left_unc'
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
write.csv(a, file='~/data/longitudinal_methylome/results/inatt_compiled_FDR.csv', row.names=F)
```

# 2020-06-29 20:27:08

I had ot make a few changes to the script to make it all run nicely. So, now
here's how to call it. I'll change it to only 8 cores so I can leave all
permutations running over night (HI each direction, plus without baseline)

```bash
cd ~/data/longitudinal_methylome
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} IN_ms_jhu IN_ys F \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_{}.csv;
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} IN_ms_jhu IN_ys T \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_{}.csv;
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv IN_ys IN_ms_jhu {} F \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_xSX_{}.csv;
cat IN_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv IN_ys IN_ms_jhu {} T \
        INATT_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_xSX_{}.csv;
```

And its similar for HI, just that it takes 3 times longer...

```bash
cd ~/data/longitudinal_methylome
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_{}.csv;
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys T \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_{}.csv;
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv HI_ys HI_ms_jhu {} F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_xSX_{}.csv;
cat HI_sets.txt | parallel -j 8 --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        final_dti_file_for_mediation_number_tracts.csv HI_ys HI_ms_jhu {} T \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_withBaseline_xSX_{}.csv;
```

And to compile it's similar:

```r
a = read.csv('~/data/longitudinal_methylome/results/hi_JHU_compiled.csv')
a$tot_p_FDR = p.adjust(a$tot_p, method='fdr')
a$acme_p_FDR = p.adjust(a$acme_p, method='fdr')
a$ade_p_FDR = p.adjust(a$ade_p, method='fdr')
a$prop_p_FDR = p.adjust(a$prop_p, method='fdr')
idx = a$M=='rd_11_rate'
a$tot_p_FDR_withinM = NA
a$acme_p_FDR_withinM = NA
a$ade_p_FDR_withinM = NA
a$prop_p_FDR_withinM = NA
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
idx = a$M=='rd_13_rate'
a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
write.csv(a, file='~/data/longitudinal_methylome/results/hi_JHU_compiled_FDR.csv', row.names=F)
```

# 2020-06-30 21:33:59

Philip asked me to re-run the mediations for probe (X) to HI (.Y) through the
DTI. Up to 10000, maybe even 100k bootstraps? Based on the results, I think he
meant the option without baseline correction. I'll run that first, only 10K
bootstraps, and if it all looks good I can change the code to run 10K, because
it will likely need to be swarmed.

```bash
head -n 109 HI_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_10K_{}.csv;
# and also
tail -n 110 HI_sets.txt | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 \
    Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv {} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_10K_{}.csv;
```

# 2020-07-01 08:19:36

So, this is not quite working out, because using 2 interactive machines of 32
cores each I'm still on 64 result files out of 219, after 10h of computation.
Let's think of how to best warm it.

It looks like the first result files were generated at 3:30am. Assuming we
started at 9:30, it's taking 6h to run a single set, which has 15 cgs. But
that's a single core.

So, I could potentially fire up all 219 at the same time, and take only 6h total
if I can allocate all swarms at the same time. And hopefully loading 219
instances of R won't be too hard on the cluster.

```bash
mydir=~/data/longitudinal_methylome
sfile=swarm.hi_10k
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        final_dti_file_for_mediation_number_tracts.csv ${s} HI_ms_jhu HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_JHU_10K2_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 --job-name med_hi --time 8:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

Philip asked me to run the cognitve data on the HI probes, which take a while.
So, let's swarm it. First, using 1K. I also change the cg_ csv he sent to cog_
and modified the baseline variables to have the same name as the ROC ones and
end with baseline.

# 2020-07-02 07:04:36

Philip is not responding about the cognitive probes, so I'll go ahead and
construct them myself. This way we can see if there is at least something here.
And I can re-run with the same data he ran with later.

```r
a = readRDS('~/data/longitudinal_methylome/UPDATED_PROCESSED_WITH_DIFF.rds')
tmp = read.table('~/data/longitudinal_methylome/HI_xs')[,1]
probes = gsub(x=tmp, pattern='_ROC', replacement='')
cog = read.csv('~/data/longitudinal_methylome/cog_mediation_for_gs.csv')
rocs = c()
for (s in unique(a$PersonID)[1:5]) {
    print(s)
    idx = which(a$PersonID == s)
    subj_ROC = c()
    for (p in probes[1:10]) {
        fit = lm(as.formula(sprintf('%s ~ ageACQ', p)),
                 data=a[idx, c('ageACQ', p)])
        subj_ROC = c(subj_ROC, summary(fit)$coefficients['ageACQ', 'Estimate'])
    }
    rocs = rbind(rocs, subj_ROC)
}
```

Nevermind, just got the file.

```r
library(data.table)
dread = fread('~/data/longitudinal_methylome/ROC_methyl.csv', header = T, sep = ',')
d = as.data.frame(dread)
colnames(d)[1] = 'PersonID'
tmp = read.table('~/data/longitudinal_methylome/HI_xs')[,1]
probes = gsub(x=tmp, pattern='_ROC', replacement='')
a = d[, c('PersonID', probes)]
roc = sapply(probes, function(x) sprintf('%s_ROC', x))
colnames(a)[2:ncol(a)] = roc
# grab baseline from the other file
b = readRDS('~/data/longitudinal_methylome/UPDATED_PROCESSED_WITH_DIFF.rds')
b2 = b[, c('PersonID', 'ageACQ', probes)]
b3 = b2[b2$PersonID %in% a$PersonID, ]
is_base = c()
for (s in unique(b3$PersonID)) {
    idx = which(b3$PersonID == s)
    if (diff(b3[idx, 'ageACQ']) > 0) {
        is_base = c(is_base, idx[1])
    } else {
        is_base = c(is_base, idx[2])
    }
}
b4 = b3[is_base, ]
bases = sapply(probes, function(x) sprintf('%s_baseline', x))
colnames(b4)[3:ncol(b4)] = bases
# grabing other variables from yet another file...
f = read.csv('~/data/longitudinal_methylome/ROC_data_pheno_file_yun_ching_used.csv')
m = merge(a, b4, by='PersonID')
m = merge(m, f, by='PersonID')
write.csv(m, file='~/data/longitudinal_methylome/ROC_data_160.csv',
          row.names=F, quote=F)
```

Now let's go ahead and swarm it. I tested with one_probe and it took 2 min. Each X
file has 15 of those, so I'll budget 1h for each subjob just to be safe.

```bash
mydir=~/data/longitudinal_methylome
sfile=swarm.cog_1k
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        cog_mediation_for_gs.csv ${s} cog_ms HI_ys F \
        ROC_data_160.csv res_cog_1K_${s}.csv;" >> $sfile;
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        cog_mediation_for_gs.csv ${s} cog_ms HI_ys T \
        ROC_data_160.csv res_cog_1K_withBase_${s}.csv;" >> $sfile;
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        cog_mediation_for_gs.csv HI_ys cog_ms ${s} F \
        ROC_data_160.csv res_cog_1K_xSX_${s}.csv;" >> $sfile;
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim_xSX.R \
        cog_mediation_for_gs.csv HI_ys cog_ms ${s} T \
        ROC_data_160.csv res_cog_1K_xSX_withBase_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name cog_hi --time 1:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

1h was cutting quite close... might need to go a bit higher if running this
again.

To compile, we do:

```bash
cd ~/data/longitudinal_methylome/results/
head -n 1 res_cog_1K_hi_bj.csv > cog_1K_compiled.csv
for f in `ls -1 res_cog_1K_hi_??.csv`; do tail -n +2 $f >> cog_1K_compiled.csv; done
```

```r
fname = '~/data/longitudinal_methylome/results/cog_1K_compiled.csv'
ms = c('win_Raw_score_SS_WISC_ROC', 'win_Raw_score_SSF_WISC_ROC',
       'win_Raw_score_SSB_WISC_ROC')
a = read.csv(fname)
a$tot_p_FDR = p.adjust(a$tot_p, method='fdr')
a$acme_p_FDR = p.adjust(a$acme_p, method='fdr')
a$ade_p_FDR = p.adjust(a$ade_p, method='fdr')
a$prop_p_FDR = p.adjust(a$prop_p, method='fdr')
a$tot_p_FDR_withinM = NA
a$acme_p_FDR_withinM = NA
a$ade_p_FDR_withinM = NA
a$prop_p_FDR_withinM = NA
for (m in ms) {
    idx = a$M==m
    a[idx,]$tot_p_FDR_withinM = p.adjust(a[idx,]$tot_p, method='fdr')
    a[idx,]$acme_p_FDR_withinM = p.adjust(a[idx,]$acme_p, method='fdr')
    a[idx,]$ade_p_FDR_withinM = p.adjust(a[idx,]$ade_p, method='fdr')
    a[idx,]$prop_p_FDR_withinM = p.adjust(a[idx,]$prop_p, method='fdr')
}
out_fname = gsub(fname, pattern='.csv', replacement='_FDR.csv')
write.csv(a, file=out_fname, row.names=F)
```
