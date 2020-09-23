# 2020-09-21 17:21:09

Let's check the Freesurfer reliability results. But first I need to process all
scans for the same IDs.

```bash
cd /mnt/NCR/sudregp/MR_data_by_maskid/
grep rage_ 2022/2*/RE*
scp -rq 2022/20160722-03043/mr_0003 helix:/scratch/sudregp/mprage/56579_v1r1
scp -rq 2022/20160722-03044/mr_0005 helix:/scratch/sudregp/mprage/56579_v1r2
grep rage_ 2097/2*/RE*
scp -rq 2097/20170512-25270/mr_0006 helix:/scratch/sudregp/mprage/56579_v2r2
scp -rq 2097/20170512-25268/mr_0011 helix:/scratch/sudregp/mprage/56579_v2r1
```

Just doing a lot of these manually, because it doesn't make sense to write a
script for it.

There are a total of 62 scans. Let's then use our normal freesurfer pipeline to
see what we get:

```bash
# biowulf
cd /scratch/sudregp/mprage
/bin/ls -1 > ~/tmp/rel.txt
cd ~/freesurfer_logs/
# do some cleanup if needed
while read s; do f=`/bin/ls /scratch/sudregp/mprage/${s}/*0001.dcm | tr -d '\n'`; echo "source /usr/local/apps/freesurfer/5.3.0/SetUpFreeSurfer.sh; recon-all -i $f -subjid ${s} -all -openmp 8 | tee -a ${s}_freesurfer.log" >> reliability.swarm; done < ~/tmp/rel.txt

swarm -g 8 -t 8 --job-name freesurfer --time 48:00:00 -f reliability.swarm -m freesurfer/5.3.0 --logdir trash
```

Then, to check who finished:

```bash
cd ~/freesurfer_logs/
rm finished.txt
while read s; do 
    if grep -q "finished without error" ~/data/MEG_structural/freesurfer/${s}/scripts/recon-all.log; then
        echo $s >> finished.txt;
    else
        echo $s;
    fi;
done < ~/tmp/rel.txt
```

I then created Shaw/reliability_scans/freesurfer_rois.csv. Now it's just a
matter of making a few plots.

```r
library(ggplot2)
df = read.csv('/Volumes/Shaw/reliability_scans/freesurfer_rois.csv')
df$subject=factor(df$subject)
ggplot(df, aes(x=1:nrow(df), y=TotalGrayVol, shape=visit, color=subject, size=2)) + geom_point()
```

![](images/2020-09-22-07-02-34.png)

![](images/2020-09-22-07-03-12.png)

![](images/2020-09-22-07-04-28.png)

![](images/2020-09-22-07-05-57.png)

I then created two files so we could run ICC. One by doubling scans that didn't
have duplicates, and one without them. The new files are sorted in ausch a way
that if I grab every other entry I'll be seeing v1 and v2 as diffeent raters.

```r
library(psych)
df = read.csv('/Volumes/Shaw/reliability_scans/freesurfer_rois_slim.csv')
x='TotalGrayVol'
tmp = cbind(df[seq(1,nrow(df),2), x], df[seq(2,nrow(df),2), x])
res = ICC(tmp)
```

I also created a one-per-subjet file.

And we can use CCC and Bland-Altman plots:

```r
x = "rh_parsopercularis_thickness"
scanner1 = df[seq(1,nrow(df),2), x]
scanner2 = df[seq(2,nrow(df),2), x]

plot_CC_BA = function(scanner1, scanner2, x) {
    tmp.ccc <- CCC(scanner1, scanner2, ci = "z-transform", conf.level = 0.95)

    lab <- paste(sprintf("CCC %s: ", x), round(tmp.ccc$rho.c[,1], digits = 2), " (95% CI ", 
    round(tmp.ccc$rho.c[,2], digits = 2), " - ",
    round(tmp.ccc$rho.c[,3], digits = 2), ")", sep = "")
    z <- lm(scanner2 ~ scanner1)

    par(mfrow=c(1, 2), pty = "s", oma = c(0, 0, 2, 0))
    plot(scanner1, scanner2, xlab = "Old scanner", 
    ylab = "New scanner", pch = 16)
    abline(a = 0, b = 1, lty = 2)
    abline(z, lty = 1)
    legend(x = "topleft", legend = c("Line of perfect concordance", 
    "Reduced major axis"), lty = c(2,1), lwd = c(1,1), bty = "n")
    mtext(lab, outer = T, cex = 1)

    tmp.ccc <- CCC(scanner1, scanner2, ci = "z-transform", conf.level = 0.95)
    tmp.mean <- mean(tmp.ccc$blalt$delta)
    tmp.sd <- sqrt(var(tmp.ccc$blalt$delta))

    plot(tmp.ccc$blalt$mean, tmp.ccc$blalt$delta, pch = 16, 
    xlab = "Average", ylab = "Difference") 
    abline(h = tmp.mean, lty = 1, col = "gray")
    abline(h = tmp.mean - (2 * tmp.sd), lty = 2, col = "gray")
    abline(h = tmp.mean + (2 * tmp.sd), lty = 2, col = "gray")
    legend(x = "topleft", legend = c("Mean difference", 
    "Mean difference +/ 2SD"), lty = c(1,2), bty = "n")
    legend(x = 0, y = 125, legend = c("Difference"), pch = 16, 
        bty = "n")
}
```

# 2020-09-23 07:17:57

The plots were useful, but didn't show any new info that wasn't there before.
What's the ICC within scanner?

```r
library(psych)
x = "lh_precentral_volume"
df = read.csv('/Volumes/Shaw/reliability_scans/freesurfer_rois_slim.csv')
tmp = cbind(df[seq(1,nrow(df),2), x], df[seq(2,nrow(df),2), x])
print('Across visits')
ICC(tmp)

df = read.csv('/Volumes/Shaw/reliability_scans/freesurfer_rois.csv')
for (v in c('v1', 'v2')) {
    print(v)
    df2 = df[df$visit==v,]
    subj_count = table(df2$subject)
    use_subjs = names(subj_count)[subj_count==2]
    v1 = c()
    v2 = c()
    for (s in use_subjs) {
        v1 = c(v1, df2[df2$subject==s & df2$run=='r1', x])
        v2 = c(v2, df2[df2$subject==s & df2$run=='r2', x])
    }
    print(ICC(cbind(v1, v2)))
}
```


