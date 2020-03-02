# 2020-03-02 09:01:41

Philip suggested I should try ML on the imputed data. If we're doing that, we
don't actually need to stick with the stacked classifier, and we can allow for
cross-domain interactions. Does that elevate our results?

```bash
g=2
cd ~/data/baseline_prediction/prs_start
my_script=~/research_code/baseline_prediction/nonstacked_${g}group_dataImpute.R;
out_file=swarm.${g}group_di
for clf in hdda rda stepLDA glmStepAIC dwdLinear bayesglm earth LogitBoost \
    kernelpls cforest; do
    for sx in inatt hi; do
        for cd in 1 2 3; do
            for cm in "T F" "T T" "F F"; do
                echo "Rscript $my_script $sx $clf $cd $cm ~/tmp/resids_${g}group_di.csv;" >> $out_file;
            done;
        done;
    done
done

swarm -g 10 -t 1 --job-name group${g} --time 20:00 -f $out_file \
    -m R --partition quick --logdir trash
```

And we can also compare it to the actual stacked classifier from before:

```bash
g=2
cd ~/data/baseline_prediction/prs_start
my_script=~/research_code/baseline_prediction/stacked_${g}group_dataImpute.R;
out_file=swarm.${g}group_diStacked
for clf in hdda rda stepLDA glmStepAIC dwdLinear bayesglm earth LogitBoost \
    kernelpls cforest; do
    for ens in rpart glm glmStepAIC rpart2 C5.0Tree; do
        for sx in inatt hi; do
            for cd in 1 2 3; do
                for cm in "T F" "T T" "F F"; do
                    echo "Rscript $my_script $sx $clf $ens $cd $cm ~/tmp/resids_${g}group_diStack.csv;" >> $out_file;
                done;
            done;
        done
    done;
done

swarm -g 10 -t 1 --job-name group${g} --time 20:00 -f $out_file \
    -m R --partition quick --logdir trash
```

Since I only ran it for the 2 group so far, let's see what's fairing best:

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/resids_2group_di.csv', header=F)
for (clf in unique(res$V2)) {
    for (cd in unique(res$V3)) {
        for (uc in unique(res$V4)) {
            for (um in unique(res$V5)) {
                idx = (res$V2 == clf & res$V3 == cd & res$V4 == uc &
                        res$V5 == um)
                pos = which(idx)
                if (length(pos) > 0) {
                    my_str = paste(c(clf, cd, uc, um), collapse='_')
                    params = c(params, my_str)
                    scores = c(scores, mean(res[pos, 'V8']))
                }
            }
        }
    }
}
print(params[which.max(scores)])
```

Our best result for nonstacked is:

```
hi	bayesglm	3	1	TRUE	2	0.956395	0.853741
inatt	bayesglm	3	1	TRUE	2	0.942115	0.616
```

Note that it's 3 year diff!

For the stacked and imputed, we get:

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/resids_2group_diStack.csv', header=F)
for (clf in unique(res$V2)) {
    for (ens in unique(res$V3)) {
        for (cd in unique(res$V4)) {
            for (uc in unique(res$V5)) {
                for (um in unique(res$V6)) {
                    idx = (res$V2 == clf & res$V4 == cd & res$V5 == uc &
                            res$V6 == um & res$V3 == ens)
                    pos = which(idx)
                    if (length(pos) > 0) {
                        my_str = paste(c(clf, ens, cd, uc, um), collapse='_')
                        params = c(params, my_str)
                        scores = c(scores, mean(res[pos, 'V9']))
                    }
                }
            }
        }
    }
}
print(params[which.max(scores)])
```

Our best one now is "stepLDA_glmStepAIC_1_TRUE_FALSE":

```
inatt	stepLDA	glmStepAIC	1	TRUE	FALSE	2	0.832644	0.708929
hi	stepLDA	glmStepAIC	1	TRUE	FALSE	2	0.807823	0.766804
```

These are a little less overfitty, and go back to clinDiff1.


