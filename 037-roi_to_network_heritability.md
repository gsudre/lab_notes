# 2019-08-30 16:11:55

As I mentioned in note 36, let me see how goo the results are if I use a ROI to
netowkr matrix. I'll try 100 and 400 ROIs, just for starters.

```R
source('~/research_code/fmri/make_schaefer_roi2nets_data_FD.R')
```

# 2019-09-05 13:06:34

Now that we have those CSV files, let's run SOLAR to see if anything is
heritable. First:

 ```bash
# desktop
# remove the double quotes from the files otherwise SOLAR can't run
cd ~/data/heritability_change/
for f in `/bin/ls rsfmri_fc-36p_despike_schaefer100roi2nets_*_09052019.csv`; do
    sed -i -e "s/\"//g" $f;
done
```

```bash
#sin
cd ~/data/heritability_change/
sfile=swarm.roi2nets100;
rm $sfile
for f in `/bin/ls rsfmri_fc-36p_despike_schaefer100roi2nets_*_09052019.csv`; do
    phen=`echo $f | sed "s/\.csv//"`;
    echo "bash ~/research_code/run_solar_parallel.sh $phen " \
        "~/data/heritability_change/roi2nets100_conns.txt" >> $sfile;
done
bash $sfile
```

And collect everything:

```bash
cd ~/data/tmp;
for f in `/bin/ls ~/data/heritability_change/rsfmri_fc-36p_despike_schaefer100roi2nets_*_09052019.csv`; do
    pheno=`echo $f | sed "s/\.csv//" | cut -d"/" -f 6`;
    echo "Working on $pheno";
    cd $pheno;
    tar -zxf *tgz;
    echo "  Compiling...";
    python ~/research_code/compile_solar_multivar_results.py $pheno;
    echo "  Cleaning up...";
    rm *.out;
    cd ..;
done
```

I can make several pictures, but for now let's just spit out how many nominally
significant connections each file has:

```r
mydir = '~/data/heritability_change/'
nnets = 7
nrois = 100
fnames = list.files(mydir, pattern='polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_.*s.*\\.csv')
roi_names = sapply(1:nrois, function(x) sprintf('roi%03d', x))
map_names = c()
sig_conns = c()
for (fname in fnames) {
    # read in the results
    res = read.csv(sprintf('%s/%s', mydir, fname))
    # figuring out possible connections
    nets = sapply(as.character(res$phen), function(x) strsplit(x, 'TO')[[1]][1])
    nets = unique(nets)        
    vals = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
    stats = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
    for (r in 1:nrow(res)) {
        ij = strsplit(as.character(res$phen[r]), 'TO')[[1]]
        rname = sprintf('roi%s', ij[2])
        cname = ij[1]
        vals[rname, cname] = res[r, 'h2r']
        stats[rname, cname] = res[r, 'h_pval']
    }
    p2 = p.adjust(stats, method='fdr')
    phen = strsplit(strtrim(fname, nchar(fname)-4), '/')[[1]]
    # sig_conns = c(sig_conns, sum(p2 < .1, na.rm=T))
    sig_conns = c(sig_conns, sum(stats < .05, na.rm=T))
    map_names = c(map_names, phen)
}
s = sort(sig_conns, index.return=T, decreasing=T)
for (i in 1:10) {
    cat(sprintf('%s: %d\n', map_names[s$ix[i]], s$x[i]))
}
```

```
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD1.00_slopes_n260_09052019: 188
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.75_slopes_n238_09052019: 153
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.50_slopes_n210_09052019: 136
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD2.50_slopes_n296_09052019: 134
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.25_slopes_n146_09052019: 130
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_FD0.25_residSlopes_n146_09052019: 64
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_FD0.25_slopes_n146_09052019: 61
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.25_residSlopes_n146_09052019: 44
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.50_residSlopes_n210_09052019: 32
polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.75_residSlopes_n238_09052019: 30
```

OK, so this is what we see for nominal p-values, out of possible 700. 

## Everything

That's a
little under 10%, and similar to what we're seeing with the voxel results, the
.25 results seem a bit more robust. Let's plot them and see if there is any sort
of interesting pattern going on.

```r
mydir = '~/data/heritability_change/'
nnets = 7
nrois = 100
fname = 'polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_FD0.25_residSlopes_n146_09052019.csv'
res = read.csv(sprintf('%s/%s', mydir, fname))
roi_names = sapply(1:nrois, function(x) sprintf('roi%03d', x))
nets = sapply(as.character(res$phen), function(x) strsplit(x, 'TO')[[1]][1])
nets = unique(nets)        
vals = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
stats = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
for (r in 1:nrow(res)) {
    ij = strsplit(as.character(res$phen[r]), 'TO')[[1]]
    rname = sprintf('roi%s', ij[2])
    cname = ij[1]
    vals[rname, cname] = res[r, 'h2r']
    stats[rname, cname] = res[r, 'h_pval']
}
p2 = p.adjust(stats, method='fdr')
phen = strsplit(strtrim(fname, nchar(fname)-4), '/')[[1]]
# sig_conns = c(sig_conns, sum(p2 < .1, na.rm=T))
sig_conns = c(sig_conns, sum(stats < .05, na.rm=T))
map_names = c(map_names, phen)
library(corrplot)
corrplot(vals, method='color', p.mat = stats, sig.level = .05,
        insig = "blank", is.corr=F, tl.cex=.8)
title(phen, cex.main=.8)
```
![](images/2019-09-05-15-12-29.png)

Well, kinda hard to see... what if we do some filtering?

```r
keep_rows = which(rowSums(stats < .05) > 0)
vals2 = vals[keep_rows, ]
stats2 = stats[keep_rows, ]
corrplot(vals2, method='color', p.mat = stats2, sig.level = .05,
        insig = "blank", is.corr=F, tl.cex=.8)
```

![](images/2019-09-05-15-16-25.png)

Maybe... we'll see. If anything, it'd be good to remove SomMot and Vis, even
though they might actually be helping with the overall distribution. We could
also keep it to only connections that all subjects have, and try Meff?

```r
corrplot(vals2[, 1:5], method='color',
         p.mat = stats2[, 1:5], sig.level = .05,
        insig = "blank", is.corr=F, tl.cex=.8)
```

![](images/2019-09-05-15-22-13.png)

OK, so let's figure out the rois that are not represented in everyone, and then
calculate Meff:

```r
df = readRDS('~/data/heritability_change/rsfmri_fc-36p_despike_schaefer100roi2nets_FD0.25_scans292_09052019.rds')
pipe_dir = '/Volumes/Shaw/rsfmri_36P/xcpengine_output_fc-36p_despike/'
nrois = 100
cat(sprintf('Reading connectivity data from %s\n', pipe_dir))
for (s in df$Mask.ID) {
    subj = sprintf('sub-%04d', s)
    fname = sprintf('%s/%s/fcon/schaefer%d/%s_schaefer%d.net',
                                pipe_dir, subj, nrois, subj, nrois)
    data = read.table(fname, skip=2)
    b = matrix(nrow=nrois, ncol=nrois)
    for (r in 1:nrow(data)) {
        b[data[r,1], data[r,2]] = data[r,3]
        b[data[r,2], data[r,1]] = data[r,3]
    }
    print(which(rowSums(is.na(b)) > 1))
}
```

Apparently, no connections were NA... that sucks. I guess the 100 Yeo ROIs are
big enough to always get a bit of the brain. Well, let's look into Meff then.

```r
fname = '~/data/heritability_change/rsfmri_fc-36p_despike_schaefer100roi2nets_FD0.25_residSlopes_n146_09052019.csv'
data = read.csv(fname)
cnames = colnames(data)[grepl(colnames(data), pattern='TO')]
cc = cor(data[, cnames])
svd = eigen(cc)
absev = abs(svd$values)
meff = (sum(sqrt(absev))^2)/sum(absev)
cat(sprintf('Galwey Meff = %.2f\n', meff))
```

So, if we include all 700 variables, we're looking at 65.22 (p = .05/65.22 =
0.000766581). If we ignore the visual and somatomotor networks (500 variables
only), we get 61.63 (p = 0.0008112449). Only one connection survives that
threshold, and it's in somatomotor anyways...

## Positive correlation only

Those results are **MUCH** stronger. 

```r
mydir = '~/data/heritability_change/'
nnets = 7
nrois = 100
fname = 'polygen_results_rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD1.00_slopes_n260_09052019.csv'
res = read.csv(sprintf('%s/%s', mydir, fname))
roi_names = sapply(1:nrois, function(x) sprintf('roi%03d', x))
nets = sapply(as.character(res$phen), function(x) strsplit(x, 'TO')[[1]][1])
nets = unique(nets)        
vals = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
stats = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
for (r in 1:nrow(res)) {
    ij = strsplit(as.character(res$phen[r]), 'TO')[[1]]
    rname = sprintf('roi%s', ij[2])
    cname = ij[1]
    vals[rname, cname] = res[r, 'h2r']
    stats[rname, cname] = res[r, 'h_pval']
}
p2 = p.adjust(stats, method='fdr')
phen = strsplit(strtrim(fname, nchar(fname)-4), '/')[[1]]
# sig_conns = c(sig_conns, sum(p2 < .1, na.rm=T))
sig_conns = c(sig_conns, sum(stats < .05, na.rm=T))
map_names = c(map_names, phen)
library(corrplot)
corrplot(vals, method='color', p.mat = stats, sig.level = .05,
        insig = "blank", is.corr=F, tl.cex=.8)
title(phen, cex.main=.8)
```

![](images/2019-09-05-16-28-49.png)

Lots of NAs as expected because of the negative connections, but I think I have
something in FDR???

```
> sum(p2 < .05, na.rm=T)
[1] 160
```

Interesting... let's see what they are...

# 2019-09-06 10:26:38

The 400 ROI analysis ran overnight. Let's see if the results are better or
worse...

```bash
cd ~/data/tmp;
for f in `/bin/ls ~/data/heritability_change/rsfmri_fc-36p_despike_schaefer400roi2nets_*_09052019.csv`; do
    pheno=`echo $f | sed "s/\.csv//" | cut -d"/" -f 6`;
    echo "Working on $pheno";
    cd $pheno;
    tar -zxf *tgz;
    echo "  Compiling...";
    python ~/research_code/compile_solar_multivar_results.py $pheno;
    echo "  Cleaning up...";
    rm *.out;
    cd ..;
done
```

As usual, let's ran the different files by how many significant connections each
one has:

```r
mydir = '~/data/heritability_change/'
nnets = 7
nrois = 400
fnames = list.files(mydir, pattern='polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_.*s.*\\.csv')
roi_names = sapply(1:nrois, function(x) sprintf('roi%03d', x))
map_names = c()
sig_conns = c()
for (fname in fnames) {
    # read in the results
    res = read.csv(sprintf('%s/%s', mydir, fname))
    # figuring out possible connections
    nets = sapply(as.character(res$phen), function(x) strsplit(x, 'TO')[[1]][1])
    nets = unique(nets)        
    vals = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
    stats = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
    for (r in 1:nrow(res)) {
        ij = strsplit(as.character(res$phen[r]), 'TO')[[1]]
        rname = sprintf('roi%s', ij[2])
        cname = ij[1]
        vals[rname, cname] = res[r, 'h2r']
        stats[rname, cname] = res[r, 'h_pval']
    }
    p2 = p.adjust(stats, method='fdr')
    phen = strsplit(strtrim(fname, nchar(fname)-4), '/')[[1]]
    sig_conns = c(sig_conns, sum(p2 < .05, na.rm=T))
    # sig_conns = c(sig_conns, sum(stats < .05, na.rm=T))
    map_names = c(map_names, phen)
}
s = sort(sig_conns, index.return=T, decreasing=T)
for (i in 1:10) {
    cat(sprintf('%s: %d\n', map_names[s$ix[i]], s$x[i]))
}
```

```
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD1.00_slopes_n260_09052019: 756
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.75_slopes_n238_09052019: 588
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.50_slopes_n210_09052019: 565
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.25_slopes_n146_09052019: 544
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD2.50_slopes_n296_09052019: 510
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.25_slopes_n146_09052019: 236
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.25_residSlopes_n146_09052019: 224
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.25_residSlopes_n146_09052019: 209
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.50_residSlopes_n210_09052019: 138
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.50_residSlopes_n210_09052019: 137
```

And that's out of 2800 connections. So, our top has 27% of significant
connections. If I look at FDR correct values instead, then I get:

```
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD1.00_slopes_n260_09052019: 615
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.75_slopes_n238_09052019: 407
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.50_slopes_n210_09052019: 386
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD0.25_slopes_n146_09052019: 333
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_posOnly_FD2.50_slopes_n296_09052019: 184
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.25_residSlopes_n146_09052019: 0
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.25_slopes_n146_09052019: 0
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.50_residSlopes_n210_09052019: 0
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.50_slopes_n210_09052019: 0
polygen_results_rsfmri_fc-36p_despike_schaefer400roi2nets_FD0.75_residSlopes_n238_09052019: 0
```

Oh wait... are we removing movement here? Yep... we were when significant! OK
then. Now we need to decide which FD threshold we'll keep, and whether we'll use
the 100 or the 400 set. We could make pictures of the heritability, then
association, then intersection, and see if there is a trend, and what looks
best.

The last stage, of course, will be to make plots of the brain with network and
ROI colors, so that we can best interpret the results.

## Association analysis

```r
library(nlme)
dd = read.csv('~/data/heritability_change/rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD1.00_slopes_n260_09052019.csv')
dd$Medical.Record...MRN = as.numeric(as.character(dd$ID))
# to get famID
tmp = read.csv('~/data/heritability_change/resting_demo_07032019.csv')
tmp$famID = sapply(1:nrow(tmp), function(x)
                                if (is.na(tmp$Extended.ID...FamilyIDs[x])) {
                                    tmp$Nuclear.ID...FamilyIDs[x]
                                }
                                else {
                                    tmp$Extended.ID...FamilyIDs[x]
                                }
                  )
tmp2 = tmp[, c('Medical.Record...MRN', 'famID')]
tmp3 = tmp2[!duplicated(tmp2[, 'Medical.Record...MRN']), ]
data = merge(dd, tmp3, by='Medical.Record...MRN', all.x=T, all.y=F)
targets = colnames(data)[grepl(colnames(data), pattern='TO')]
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline' )
for (t in predictors) {
    data[, t] = as.numeric(as.character(data[, t]))
}
out_fname = '~/data/heritability_change/assoc_LME_100posOnlyFD1.csv'

predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline', 'DX',
               'DX2')
hold=NULL
for (i in targets) {
    cat(sprintf('%s\n', i))
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex + qc', i, j)
        model1<-try(lme(as.formula(fm_str), data, ~1|famID, na.action=na.omit))
        if (length(model1) > 1) {
            temp<-summary(model1)$tTable
            a<-as.data.frame(temp)
            a$formula<-fm_str
            a$target = i
            a$predictor = j
            a$term = rownames(temp)
            hold=rbind(hold,a)
        } else {
            hold=rbind(hold, NA)
        }
    }
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX=='ADHD', ]
out_fname = gsub('.csv', x=out_fname, '_dx1.csv')
predictors = c('SX_inatt', 'SX_HI', 'inatt_baseline', 'HI_baseline' )
hold=NULL
for (i in targets) {
    cat(sprintf('%s\n', i))
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex + qc', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|famID, na.action=na.omit))
        if (length(model1) > 1) {
            temp<-summary(model1)$tTable
            a<-as.data.frame(temp)
            a$formula<-fm_str
            a$target = i
            a$predictor = j
            a$term = rownames(temp)
            hold=rbind(hold,a)
        } else {
            hold=rbind(hold, NA)
        }
    }
}
write.csv(hold, out_fname, row.names=F)

data2 = data[data$DX2=='ADHD', ]
out_fname = gsub('dx1', x=out_fname, 'dx2')
hold=NULL
for (i in targets) {
    cat(sprintf('%s\n', i))
    for (j in predictors) {
        fm_str = sprintf('%s ~ %s + sex + qc', i, j)
        model1<-try(lme(as.formula(fm_str), data2, ~1|famID, na.action=na.omit))
        if (length(model1) > 1) {
            temp<-summary(model1)$tTable
            a<-as.data.frame(temp)
            a$formula<-fm_str
            a$target = i
            a$predictor = j
            a$term = rownames(temp)
            hold=rbind(hold,a)
        } else {
            hold=rbind(hold, NA)
        }
    }
}
write.csv(hold, out_fname, row.names=F)
```

I'm running the code above in separate R instances to go faster. But then, the
idea is to compile the data in the same ways we were compiling the heritability
results. Except that now we have multiple different predictors to account for:

```r
mydir = '~/data/heritability_change/'
nnets = 7
nrois = 100
fnames = list.files(mydir, pattern='assoc_LME_100posOnlyFD.*\\.csv')
roi_names = sapply(1:nrois, function(x) sprintf('roi%03d', x))
for (pred in c('SX_inatt', 'SX_HI')) {
    cat(sprintf('predictor: %s\n', pred))
    map_names = c()
    sig_conns = c()
    sig_conns_fdr = c()
    for (fname in fnames) {
        # read in the results
        res = read.csv(sprintf('%s/%s', mydir, fname))
        # figuring out possible connections
        nets = sapply(as.character(res$target), function(x) strsplit(x, 'TO')[[1]][1])
        nets = unique(nets)
        nets = nets[!is.na(nets)]
        vals = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
        stats = matrix(nrow=nrois, ncol=nnets, dimnames=list(roi_names, nets))
        mres = res[which(res$predictor==pred & res$term==pred), ]
        for (r in 1:nrow(mres)) {
            ij = strsplit(as.character(mres$target[r]), 'TO')[[1]]
            rname = sprintf('roi%s', ij[2])
            cname = ij[1]
            vals[rname, cname] = mres[r, 't.value']
            stats[rname, cname] = mres[r, 'p.value']
        }
        p2 = p.adjust(stats, method='fdr')
        phen = strsplit(strtrim(fname, nchar(fname)-4), '/')[[1]]
        sig_conns_fdr = c(sig_conns_fdr, sum(p2 < .05, na.rm=T))
        sig_conns = c(sig_conns, sum(stats < .05, na.rm=T))
        map_names = c(map_names, phen)
    }
    cat('Nominal\n')
    s = sort(sig_conns, index.return=T, decreasing=T)
    for (i in 1:10) {
        cat(sprintf('%s: %d\n', map_names[s$ix[i]], s$x[i]))
    }
    cat('FDR\n')
    s = sort(sig_conns_fdr, index.return=T, decreasing=T)
    for (i in 1:10) {
        cat(sprintf('%s: %d\n', map_names[s$ix[i]], s$x[i]))
    }
}
```

Association after FDR is getting a bit difficult, but maybe we can restrict it
to only the heritable connections? This way we'll have less comparison to deal
with.

In any case, it still worked for p25, which might be what we end up with anyways
because of the stricter threshold on movement:

```
predictor: SX_inatt
Nominal
assoc_LME_100posOnlyFDp25: 124
assoc_LME_100posOnlyFD1: 101
assoc_LME_100posOnlyFDp5: 77
assoc_LME_100posOnlyFDp75: 60
assoc_LME_100posOnlyFD1_dx2: 41
assoc_LME_100posOnlyFD1_dx1: 38
assoc_LME_100posOnlyFDp75_dx2: 18
assoc_LME_100posOnlyFDp5_dx2: 16
assoc_LME_100posOnlyFDp25_dx1: 14
assoc_LME_100posOnlyFDp75_dx1: 14
FDR
assoc_LME_100posOnlyFDp25: 7
assoc_LME_100posOnlyFDp25_dx1: 1
assoc_LME_100posOnlyFD1_dx1: 0
assoc_LME_100posOnlyFD1_dx2: 0
assoc_LME_100posOnlyFD1: 0
assoc_LME_100posOnlyFDp25_dx2: 0
assoc_LME_100posOnlyFDp5_dx1: 0
assoc_LME_100posOnlyFDp5_dx2: 0
assoc_LME_100posOnlyFDp5: 0
assoc_LME_100posOnlyFDp75_dx1: 0
predictor: SX_HI
Nominal
assoc_LME_100posOnlyFDp75: 91
assoc_LME_100posOnlyFD1: 64
assoc_LME_100posOnlyFD1_dx2: 55
assoc_LME_100posOnlyFDp5: 39
assoc_LME_100posOnlyFDp25: 32
assoc_LME_100posOnlyFDp75_dx2: 32
assoc_LME_100posOnlyFD1_dx1: 27
assoc_LME_100posOnlyFDp5_dx2: 26
assoc_LME_100posOnlyFDp5_dx1: 22
assoc_LME_100posOnlyFDp25_dx2: 17
FDR
assoc_LME_100posOnlyFD1_dx1: 0
assoc_LME_100posOnlyFD1_dx2: 0
assoc_LME_100posOnlyFD1: 0
assoc_LME_100posOnlyFDp25_dx1: 0
assoc_LME_100posOnlyFDp25_dx2: 0
assoc_LME_100posOnlyFDp25: 0
assoc_LME_100posOnlyFDp5_dx1: 0
assoc_LME_100posOnlyFDp5_dx2: 0
assoc_LME_100posOnlyFDp5: 0
assoc_LME_100posOnlyFDp75_dx1: 0
```

With 400 connections, we get:

```
predictor: SX_inatt
Nominal
assoc_LME_400posOnlyFD1: 350
assoc_LME_400posOnlyFDp25: 342
assoc_LME_400posOnlyFDp75: 271
assoc_LME_400posOnlyFD1_dx2: 179
assoc_LME_400posOnlyFD1_dx1: 116
assoc_LME_400posOnlyFDp5: 101
assoc_LME_400posOnlyFDp75_dx2: 82
assoc_LME_400posOnlyFDp25_dx1: 63
assoc_LME_400posOnlyFDp75_dx1: 47
assoc_LME_400posOnlyFDp5_dx2: 41
FDR
assoc_LME_400posOnlyFDp25: 4
assoc_LME_400posOnlyFDp25_dx1: 3
assoc_LME_400posOnlyFDp25_dx2: 2
assoc_LME_400posOnlyFD1_dx1: 1
assoc_LME_400posOnlyFDp75_dx1: 1
assoc_LME_400posOnlyFD1_dx2: 0
assoc_LME_400posOnlyFD1: 0
assoc_LME_400posOnlyFDp5_dx1: 0
assoc_LME_400posOnlyFDp5_dx2: 0
assoc_LME_400posOnlyFDp5: 0
predictor: SX_HI
Nominal
assoc_LME_400posOnlyFDp75: 265
assoc_LME_400posOnlyFD1: 259
assoc_LME_400posOnlyFD1_dx2: 208
assoc_LME_400posOnlyFD1_dx1: 132
assoc_LME_400posOnlyFDp25: 112
assoc_LME_400posOnlyFDp75_dx2: 110
assoc_LME_400posOnlyFDp75_dx1: 68
assoc_LME_400posOnlyFDp5: 64
assoc_LME_400posOnlyFDp25_dx1: 62
assoc_LME_400posOnlyFDp5_dx2: 55
FDR
assoc_LME_400posOnlyFDp75_dx1: 3
assoc_LME_400posOnlyFDp75: 1
assoc_LME_400posOnlyFD1_dx1: 0
assoc_LME_400posOnlyFD1_dx2: 0
assoc_LME_400posOnlyFD1: 0
assoc_LME_400posOnlyFDp25_dx1: 0
assoc_LME_400posOnlyFDp25_dx2: 0
assoc_LME_400posOnlyFDp25: 0
assoc_LME_400posOnlyFDp5_dx1: 0
assoc_LME_400posOnlyFDp5_dx2: 0
```

OK, so let's do a more principled approach. First, remove all the NAs. When
dealing with positive only, remove all the NAs. Then, keep only the rows that
are significant using FDR. Then, look at association only for those connections.

```r
fd = 1
nrois = 100

mydir = '~/data/heritability_change/'
fmask = sprintf('polygen_results_.*er%droi2nets_posOnly_FD%.2f_slop.*s.*\\.csv', nrois, fd)
fname = list.files(mydir, pattern=fmask)[1]  # has SOLAR results

res = read.csv(sprintf('%s/%s', mydir, fname))
# figuring out possible connections
bad_conns = is.na(res$h_pval)
res = res[!bad_conns, ]
res$p2 = p.adjust(res$h_pval, method='fdr')
cat(sprintf('Nominal p<.05 heritable: %d out of %d\n', sum(res$h_pval<.05), nrow(res)))
cat(sprintf('FDR q<.05 heritable: %d\n', sum(res$p2<.05)))

nets = sapply(as.character(res$phen), function(x) strsplit(x, 'TO')[[1]][1])
nets = unique(nets)
roi_names = sapply(as.character(res$phen), function(x) strsplit(x, 'TO')[[1]][2])
roi_names = sapply(unique(roi_names), function(x) sprintf('roi%s', x))

vals = matrix(nrow=length(roi_names), ncol=length(nets), dimnames=list(roi_names, nets))
stats = matrix(nrow=length(roi_names), ncol=length(nets), dimnames=list(roi_names, nets))
for (r in 1:nrow(res)) {
    ij = strsplit(as.character(res$phen[r]), 'TO')[[1]]
    rname = sprintf('roi%s', ij[2])
    cname = ij[1]
    vals[rname, cname] = res[r, 'h2r']
    stats[rname, cname] = res[r, 'p2']
}
phen = strsplit(strtrim(fname, nchar(fname)-4), '/')[[1]]

s = sort(sig_conns, index.return=T, decreasing=T)
for (i in 1:10) {
    cat(sprintf('%s: %d\n', map_names[s$ix[i]], s$x[i]))
}
```

I had to stop this analysis because I was getting several results with
heritability==1, which is quite weird (e.g. ContTO012 and ContTO020 in
rsfmri_fc-36p_despike_schaefer100roi2nets_posOnly_FD0.25_slopes_n146_09052019.csv).
Need to check what's going on there. Yes, it looks like lots of NAs!

So, for posOnly we get lots of NAs whenever one of the two scans had NA for the
connection. Also, that's why it's still possible to have negative data, because
it's the slope of two positive connections over time. In any case, if we have
many NAs for any connection, then the heritability estimate might be unstable.
In that situation we're getting too high heritability values. So, we'll need to
figure out the proportion of NAs we can still rely on for the heritability
estimate. So, far I'm seeing more than 63% causing instability, but that's in
the 146 dataset. 

But it doesn't look like SOLAR is pulling the correct number of individuals
either... OK, I was missing 3 of them from the pedigree file. Actually, let's
make sure everyone in the FD2.5 file is there, just in case.

These were missing:

7766622
7212288
7762197
7746465
4540621
7600811

At least for DTI they're all there.

OK, so I fixed the pedigree file, and I noticed that I need to make sure all my
csv fileas have NA in them, not NaN! Also, I should probably change the
compilation function as well to grab the number of individuals we're using, just
in case... the H2 is still high, but the p-value isn't as good anymore.




