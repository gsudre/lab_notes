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

# 2020-07-06 19:57:07

Running some methyl regressions to be the icing on the cake in the paper.

```r
cog = read.csv('~/data/longitudinal_methylome/cog_mediation_for_gs.csv')
a = readRDS('~/data/longitudinal_methylome/UPDATED_PROCESSED_WITH_DIFF.rds')
idx = a$PersonID %in% cog$PersonID
a2 = a[idx, ]
saveRDS(a2, file='~/data/longitudinal_methylome/UPDATED320.rds')
cg_vars = colnames(a2)[grepl(colnames(a2), pattern='^cg')]
a3 = a2[, c('PersonID', 'ageACQ.1', cg_vars)]
tmp = c()
i = 1
while (i < nrow(a3)) {
    if (a3[i, 'ageACQ.1'] < a3[i+1, 'ageACQ.1']) {
        tmp = c(tmp, 1, 2)
    } else {
        tmp = c(tmp, 2, 1)
    }
    i = i + 2
}
a3 = cbind(a3, tmp)
a4 = reshape(a3, idvar = "PersonID", v.names=cg_vars,
           timevar = "tmp",
           direction = "wide")
```

This is not working... it's taking forever to make it wide.

Let's try something else. Can we reshape per cg? Then I can use all cores to do
what's needed.

```r
a = readRDS('~/data/longitudinal_methylome/UPDATED320.rds')
tmp = a$ageACQ
for (s in unique(a$PersonID)) {
    subj_rows = which(a$PersonID==s)
    if (a[subj_rows[1], 'ageACQ'] < a[subj_rows[2], 'ageACQ']) {
        tmp[subj_rows] = c(1, 2)
    } else {
        tmp[subj_rows] = c(2, 1)
    }
}
a = cbind(a, tmp)
f = read.csv('~/data/longitudinal_methylome/ROC_data_pheno_file_yun_ching_used.csv')

fit_wrapper = function(x) {
    a_slim = a[, c('PersonID', 'tmp', x)]
    a_wide = reshape(a_slim, idvar = "PersonID",
                     timevar = "tmp",
                     direction = "wide")
    m = merge(a_wide, f, by='PersonID')
    # f_str = '%s.2 ~ %s.1 + age.diff + ageACQ.1 + sample_type + PC1 + PC2 + PC3 + PC4 + PC5 + SV.one.m2 + CD8T.diff + CD4T.diff + NK.diff + Bcell.diff + Mono.diff + Gran.diff + sex'
    f_str = '%s.2 ~ %s.1 + age.diff + ageACQ.1 + PC1 + PC2 + PC3 + PC4 + PC5 + SV.one.m2 + CD8T.diff + CD4T.diff + NK.diff + Bcell.diff + Mono.diff + Gran.diff + sex'
    # f_str = '%s.2 ~ %s.1 + sample_type + age.diff + ageACQ.1 + PC1 + PC2 + PC3 + PC4 + PC5 + SV.one.m2 + CD8T.diff + CD4T.diff + NK.diff + Bcell.diff + Mono.diff + Gran.diff + sex'
    fit = lm(as.formula(sprintf(f_str, x, x)), data=m)
    res = summary(fit)$coefficients
    myrow = which(grepl(rownames(res), pattern=sprintf('^%s', x)))
    return(res[myrow, ])
}

cg_vars = colnames(a)[grepl(colnames(a), pattern='^cg')]
Ms = cg_vars[1:100]
m1_res = lapply(Ms, fit_wrapper)
all_res = do.call(rbind, m1_res)

# taking too long to export the big a variable to the cluster. Will try to run it
# in individual CPUs to see how long it'll take.
library(parallel)
cl <- makeCluster(2, type='FORK')
m1_res2 = parLapply(cl, Ms, fit_wrapper)
all_res2 = do.call(rbind, m1_res2)
stopCluster(cl)
```

This is still taking a long time. Let's swarm it using the script
methyl_lm_wrapper.R.

It took me 32s to run 96 probes. Rounding up to one minute, I can do 23K probes
in 4h, so I'll split 20K probes per file.

```bash
mydir=~/data/longitudinal_methylome
cd ${mydir}
for v in '' '_blood' '_saliva' '_inter'; do
    sfile=swarm.lm${v}
    rm -rf $sfile
    for s in `ls methyls_??`; do
        echo "cd ${mydir}; Rscript ~/research_code/methyl_lm_wrapper${v}.r $s 32 ~/tmp/adjusted.csv" >> $sfile;
    done;
    swarm -g 100 -t 32 --job-name lm${v} --time 4:00:00 -f ${sfile} \
        -m R --partition quick,norm --logdir trash --gres=lscratch:40
done
```

Actually it didn't take very long to run it in interactive... only an hour.
Let's keep running it that way (with the parallel package).

Now, just some clean up and calculating FDR:

```bash
grep -v sample interactions.csv > inter.csv
sed -i -e 's/\.1\,/\,/g' inter.csv
```

```r
a = read.csv('~/data/longitudinal_methylome/blood.csv')
a$FDR_q = p.adjust(a[, 5], method='fdr')
write.csv(a, file='~/data/longitudinal_methylome/blood_FDR.csv', row.names=F, quote=F)
```

# 2020-08-03 12:51:22

Philip sent a new file to run. Let's go through the same steps as above. I also
changed the script to run a version without x_base, just in case.

```bash
mydir=~/data/longitudinal_methylome
sfile=swarm.hiNB2_1k
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        HI_ROC_methyl_sx_dti_82.csv res_2_1K_noBaseX_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name hiNB2 --time 4:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
sfile=swarm.inattNB2_1k
cd ${mydir}
rm -rf $sfile
for s in `cat IN_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        INATT_ROC_methyl_sx_dti_82.csv res_2_1K_noBaseX_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name inattNB2 --time 4:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

Let's work on compiling the results:

```bash
cd ~/data/longitudinal_methylome/results/
head -n 1 res_2_1K_hi_bj.csv > res_2_1K_hi_compiled.csv
for f in `ls -1 res_2_1K_hi_??.csv`; do tail -n +2 $f >> res_2_1K_hi_compiled.csv; done
```

```r
fname = '~/data/longitudinal_methylome/results/res_2_1K_hi_compiled.csv'
ms = c('AD_left_ifo_rate', 'RD_left_ifo_rate', 'AD_left_ilf_rate',
       'AD_left_slf_rate', 'RD_left_slf_rate', 'RD_right_ilf_rate',
       'AD_right_slf_rate')
ms = c('AD_left_unc_rate', 'AD_right_unc_rate', 'RD_right_unc_rate')
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

# 2020-08-04 08:10:00

The cg file I was using didn't have everyone, so I'll have to re-create the file
like above, but with both sets of cgs (above only had HI):

```r
library(data.table)
dread = fread('~/data/longitudinal_methylome/ROC_methyl.csv', header = T, sep = ',')
d = as.data.frame(dread)
colnames(d)[1] = 'PersonID'
tmp1 = read.table('~/data/longitudinal_methylome/HI_probes_3294.txt')[,1]
tmp2 = read.table('~/data/longitudinal_methylome/IN_probes_1117.txt')[,1]
probes = unique(c(tmp1, tmp2))
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
write.csv(m, file='~/data/longitudinal_methylome/ROC_data_inattAndHI_160.csv',
          row.names=F, quote=F)
```

Now we need to re-split the new probe files that contian the ch. probes. I also
changed the mediation code to accomodate those.

```bash
cd ~/data/longitudinal_methylome
sed -e "s/$/_ROC/g" HI_probes_3294.txt > HI_probes_3294_ROC.txt
sed -e "s/$/_ROC/g" IN_probes_1117.txt > IN_probes_1117_ROC.txt
split -l 15 IN_probes_1117_ROC.txt inatt_
ls -1 inatt_* > IN_sets.txt;
split -l 15 HI_probes_3294_ROC.txt hi_
ls -1 hi_* > HI_sets.txt;
```

And finally we swarm everything:

```bash
mydir=~/data/longitudinal_methylome
jname=inatt90_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat IN_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_inattAndHI_160.csv res90_2_1K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=hi90_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_inattAndHI_160.csv res90_2_1K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=inattNB90_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat IN_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_inattAndHI_160.csv resNB90_2_1K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=hiNB90_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_inattAndHI_160.csv resNB90_2_1K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
    -m R --partition quick,norm --logdir trash
```

For 10K I'll have to increase the wall time:

```bash
mydir=~/data/longitudinal_methylome
jname=inatt90_10k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat IN_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_inattAndHI_160.csv res90_2_10K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 16:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=hi90_10k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_inattAndHI_160.csv res90_2_10K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 16:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=inattNB90_10k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat IN_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_inattAndHI_160.csv resNB90_2_10K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 16:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=hiNB90_10k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat HI_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_inattAndHI_160.csv resNB90_2_10K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 16:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

Philip asked me to change the code to spit out the regression values as well.
Since that doesn't depend on the perms, I'll re-run everything with a couple
boostraps just so it goes fast and we can have a file that is just useful for
the regressions values. In the future, when we run the code we'll have
everything in the same file though.

```bash
cd ~/data/longitudinal_methylome
Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv HI_probes_3294_ROC.txt HI_ms2 HI_ys F \
        ROC_data_inattAndHI_160.csv resNB90_2_lmValuesOnly.csv
Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv IN_probes_1117_ROC.txt IN_ms2 IN_ys F \
        ROC_data_inattAndHI_160.csv resNB90_2_inatt_lmValuesOnly.csv
```

Then change the code to run:

```bash
cd ~/data/longitudinal_methylome
Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv HI_probes_3294_ROC.txt HI_ms2 HI_ys F \
        ROC_data_inattAndHI_160.csv res90_2_hi_lmValuesOnly.csv
`Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv IN_probes_1117_ROC.txt IN_ms2 IN_ys F \
        ROC_data_inattAndHI_160.csv res90_2_inatt_lmValuesOnly.csv
```

# 2020-08-06 16:15:27

Philip asked me to run all probes this time. 1K first, no baseline.

```r
library(data.table)
dread = fread('~/data/longitudinal_methylome/ROC_methyl.csv', header = T, sep = ',')
d = as.data.frame(dread)
colnames(d)[1] = 'PersonID'
probes = colnames(d)[2:ncol(d)]
a = d
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
saveRDS(m, file='~/data/longitudinal_methylome/ROC_data_ALL_160.rds')
```

```bash
cd ~/data/longitudinal_methylome
split -l 15 all_probes.txt probes_
ls -1 probes_* > ALL_sets.txt;
```

```bash
mydir=~/data/longitudinal_methylome
jname=inatt90allProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inatt90allProbes_2_1K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=hi90allProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90allProbes_2_1K_${s}.csv;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
    -m R --partition norm --logdir trash
```

Biowulf is estimating 4.5days because of the autobundle... let me see if I can
compact this a bit more. So, 1K for HI, which has 7Ms, took a little over 2h for
15 probes. Say it took 2.5h, so 30 probes in 5h, 45 comfortably in 8h.

```bash
cd ~/data/longitudinal_methylome
rm probes_*
split -l 45 all_probes.txt probes_
ls -1 probes_* > ALL_sets.txt;
```

We will still have about 19K lines in the swarm file, but it might help...
```bash
mydir=~/data/longitudinal_methylome
jname=hi90allProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90allProbes_2_1K_${s}.csv;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 8:00:00 -f ${f} \
            -m R --partition norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Keeping track of last submitted because I mistakenly ran the loop in an
interactive node... last submitted: ak

The inattention runs, with 15 probes but only 2 Ms, took less than 1h. So, I
could resplit the probes to keep it at 8h, or just run as many but set it to 4h
instead... 

# 2020-08-07 14:08:44

Philip asked me to run a random subsample of about 100K methyls, including the
good ones. Let's pre-select them. So, I created a file in R called
probes100KandGood.txt, which lists 100K probes selected at random plus 97 of the
108 good probes that hadn't been randomly selected. Let's run those probes now.

```bash
cd ~/data/longitudinal_methylome
rm probes_*
split -l 45 probes100KandGood.txt probes_
ls -1 probes_* > ALL_sets.txt;
```

We now have about 2225 lines in the swarm file, but it might help...

```bash
mydir=~/data/longitudinal_methylome
jname=hi90SomeProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90SomeProbes_2_1K_${s}.csv;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 16:00:00 -f ${f} \
            -m R --partition norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

Since now we only have 3 split swarm files, I can run inatt as well. And since
it only has 2 Ms, it can be only 4h.

```bash
mydir=~/data/longitudinal_methylome
jname=inatt90SomeProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inatt90SomeProbes_2_1K_${s}.csv;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 8:00:00 -f ${f} \
            -m R --partition norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

# 2020-08-09 12:18:26

I made a few changes to the code to make it even faster, because now we're
trimming the big file right away. Let's first run associations that include the
x_base. I'll run them in the laptop, and start with the good probes, then 100K,
and finally 800K.

```bash
Rscript ~/research_code/mediation_code_for_methylation_slim.R \
    dti_2_for_sam_slim.csv probes.txt IN_ms2 IN_ys F \
    ROC_data_ALL_160.rds res_inatt90SigProbes_linearOnly.csv 2

Rscript ~/research_code/mediation_code_for_methylation_slim.R \
    dti_2_for_sam_slim.csv probes.txt HI_ms2 HI_ys F \
    ROC_data_ALL_160.rds res_hi90SigProbes_linearOnly.csv 2
```

For the 100K randomly selected probes, I'll have to use the chunks like before
so that it's fast inside the script:

```bash
cd ~/data/longitudinal_methylome
rm probes_*
split -l 45 probes100KandGood.txt probes_
ls -1 probes_* > SOME_sets.txt;
```

We now have about 2225 lines in the swarm file, but it might help...

```bash
mydir=~/data/longitudinal_methylome
jname=hi90SomeProbes
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90SomeProbes_linearOnly_${s}.csv 2;" >> $sfile;
done
cat $sfile | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 {};
```

```bash
mydir=~/data/longitudinal_methylome
jname=inatt90SomeProbes
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inatt90SomeProbes_linearOnly_${s}.csv 2;" >> $sfile;
done
cat $sfile | parallel -j $SLURM_CPUS_PER_TASK --max-args=1 {};
```

And now we should be fine to re-run the mediations within the amount of time we
expect them to take:

```bash
mydir=~/data/longitudinal_methylome
jname=inatt90SomeProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inatt90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 8:00:00 -f ${f} \
            -m R --partition norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
mydir=~/data/longitudinal_methylome
jname=hi90SomeProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 16:00:00 -f ${f} \
            -m R --partition norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

And compiling the lmOnly results:

```bash
cd ~/data/longitudinal_methylome/results/
head -n 1 res_inatt90SomeProbes_linearOnly_probes_aa.csv > res_inatt90SomeProbes_linearOnly_compiled.csv
for f in `ls -1 res_inatt90SomeProbes_linearOnly_probes_*.csv`; do echo $f; tail -n +2 $f >> res_inatt90SomeProbes_linearOnly_compiled.csv; done
```

I ran 3 probes for 1K perms in inattention, and it took 7 min. Multiply it by 15
for the 45 probes we have, so we have 105 min for inattention. We can add some delay
of handling the bigger file and any other shenanigans, and 3 to 4h should be
enough. Actually, some have already finished before I ran this test, and they
took about 2.5h. For hi, 3 probes ran in 14min, so 45 would be 210min, which is
3.5h. Let's go for 6h to be safe.

Some of the runs failed... let's see which ones.

```bash
mydir=~/data/longitudinal_methylome
jname=inatt_redo
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    if [ ! -e res_inatt90SomeProbes_1K_${s}.csv ]; then
        echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inatt90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
    fi;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${sfile} \
            -m R --partition quick,norm --logdir trash
```

And I created a noXBase version of the code so I don't have to keep waiting for
everything to get queued:

```bash
mydir=~/data/longitudinal_methylome
jname=inattNXB90SomeProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim_noXBase.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inattNXB90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${f} \
            -m R --partition quick,norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

```bash
mydir=~/data/longitudinal_methylome
jname=hiNXB90SomeProbes_1k
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim_noXBase.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hiNXB90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
done
split -l 1000 $sfile ${jname}_split;
for f in `/bin/ls ${jname}_split??`; do
    echo "ERROR" > swarm_wait_${USER}
    while grep -q ERROR swarm_wait_${USER}; do
        echo "Trying $f"
        swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 6:00:00 -f ${f} \
            -m R --partition norm --logdir trash 2> swarm_wait_${USER};
        if grep -q ERROR swarm_wait_${USER}; then
            echo -e "\tError, sleeping..."
            sleep 10m;
        fi;
    done;
done
```

And Philip asked me to run the linear regression for all probe sets. I want to
keep it to a single swarm file, to let's keep only 1000 probe sets.

```bash
cd ~/data/longitudinal_methylome
rm all_probes_*
split -l 900 all_probes.txt all_probes_
ls -1 all_probes_* > ALL_sets.txt;

mydir=~/data/longitudinal_methylome
jname=inatt90ALLlmOnly
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inatt90All_lmOnly_${s}.csv 1;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${f} \
            -m R --partition quick,norm --logdir trash

jname=hi90ALLlmOnly
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat ALL_sets.txt`; do
    echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90All_lmOnly_${s}.csv 1;" >> $sfile;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 4:00:00 -f ${f} \
            -m R --partition quick,norm --logdir trash
```

That took only 8min for the entire file in inatt, so 1h should be plenty, and 2
for hi.

# 2020-08-10 06:40:12

Figuring out what else left to run:

```bash
mydir=~/data/longitudinal_methylome
jname=hi_redo
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    if [ ! -e res_hi90SomeProbes_1K_${s}.csv ]; then
        echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hi90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
    fi;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 8:00:00 -f ${sfile} \
            -m R --partition norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=inattNXB_redo
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    if [ ! -e res_inattNXB90SomeProbes_1K_${s}.csv ]; then
        echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim_noXBase.R \
        dti_2_for_sam_slim.csv ${s} IN_ms2 IN_ys F \
        ROC_data_ALL_160.rds res_inattNXB90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
    fi;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 6:00:00 -f ${sfile} \
            -m R --partition norm --logdir trash
```

```bash
mydir=~/data/longitudinal_methylome
jname=hiNXB_redo
sfile=swarm.${jname}
cd ${mydir}
rm -rf $sfile
for s in `cat SOME_sets.txt`; do
    if [ ! -e res_hiNXB90SomeProbes_1K_${s}.csv ]; then
        echo "cd ${mydir}; Rscript ~/research_code/mediation_code_for_methylation_slim_noXBase.R \
        dti_2_for_sam_slim.csv ${s} HI_ms2 HI_ys F \
        ROC_data_ALL_160.rds res_hiNXB90SomeProbes_1K_${s}.csv 1000;" >> $sfile;
    fi;
done
swarm -g 12 -t 1 -p 2 --job-name ${jname} --time 10:00:00 -f ${sfile} \
            -m R --partition norm --logdir trash
```

And we compile the ones that finished:

```bash
cd ~/data/longitudinal_methylome/
head -n 1 res_inatt90SomeProbes_1K_probes_aa.csv > res_inatt90SomeProbes_1K_compiled.csv
for f in `ls -1 res_inatt90SomeProbes_1K_probes_*.csv`; do tail -n +2 $f >> res_inatt90SomeProbes_1K_compiled.csv; done
```

```r
fname = '~/data/longitudinal_methylome/results/res_inatt90SomeProbes_1K_compiled.csv'
# hi
ms = c('AD_left_ifo_rate', 'RD_left_ifo_rate', 'AD_left_ilf_rate',
       'AD_left_slf_rate', 'RD_left_slf_rate', 'RD_right_ilf_rate',
       'AD_right_slf_rate')
# inatt
ms = c('AD_left_unc_rate', 'AD_right_unc_rate', 'RD_right_unc_rate')
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
