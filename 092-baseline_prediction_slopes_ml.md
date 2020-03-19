# 2020-03-17 20:53:35

Using the new files Philip has generated. He described the predictors:

```
It has 8 anatomic predictors (all residualized); lots of PRS (keep them all I
think); 5 cognitive (residualized and around 10% imputed); and 3 demographics
(SES_group3, population_self2 and sex_numeric).
```

```r
library(caret)
data = read.csv('~/Downloads/gf_impute_based_anatomy_272.csv')
# so they don't get rescaled
data$sex_numeric = as.factor(data$sex_numeric)
data$SES_group3 = as.factor(data$SES_group3)
var_names = colnames(data)[c(10:17, 18:29, 30:34, 4:6)]
phen = "slope_inatt_res_trim.x"
reg_model = 'elasticnet'

set.seed(42)
fitControl <- trainControl(method = "repeatedcv",
                           number = 10,
                           repeats = 10)

# if imputation is fair game before we split, so is scaling
scale_me = c()
for (v in var_names) {
    if (!is.factor(data[, v])) {
        scale_me = c(scale_me, v)
    } else {
        data[, v] = as.numeric(data[, v])
    }
}
data[, scale_me] = scale(data[, scale_me])

set.seed(42)
fit <- train(x=data[, var_names],
             y=data[, phen],
             method = reg_model,
             trControl = fitControl,
             tuneLength = 10)

line=sprintf("%s,%s,%f,%f", phen, reg_model,
             mean(fit$results$RMSE), sd(fit$results$RMSE))
print(line)
```

This seems to be working. Now, there are several things we can script out, so
let's create the script before we farm it out.

So, we'd do something like:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/nonstacked_slope_dataImpute.R;
out_file=swarm.slope_impInter
rm $out_file
for clf in `cat all_reg.txt`; do
    for sx in inatt hi; do
        for fname in anatomy_272 dti_165; do 
            for fold in "10 10" "5 5" "3 10"; do
                echo "Rscript $my_script ${my_dir}/gf_impute_based_${fname}.csv $sx $clf $fold ~/tmp/residsFixed_slope_impInter.csv;" >> $out_file;
            done
        done;
    done;
done

swarm -g 10 -t 1 --job-name interFixedSlope --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

Now we check our best results:

```r
res = read.csv('~/data/baseline_prediction/prs_start/residsFixed_slope_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'fname', 'nfolds', 'nreps', 'meanRMSE', 'sdRMSE')
res[which.min(res$meanRMSE),]
```

And we can also look at R2:

```bash
my_dir=~/data/baseline_prediction/prs_start
cd $my_dir
my_script=~/research_code/baseline_prediction/nonstacked_slope_dataImpute_R2.R;
out_file=swarm.slope_impInterR2
rm $out_file
for clf in `cat all_reg.txt`; do
    for sx in inatt hi; do
        for fname in anatomy_272 dti_165; do 
            for fold in "10 10" "5 5" "3 10"; do
                echo "Rscript $my_script ${my_dir}/gf_impute_based_${fname}.csv $sx $clf $fold ~/tmp/residsR2_slope_impInter.csv;" >> $out_file;
            done
        done;
    done;
done

swarm -g 10 -t 1 --job-name interR2Slope --time 4:00:00 -f $out_file \
    -m R --partition quick --logdir trash
```

# 2020-03-19 07:22:57

```r
res = read.csv('~/data/baseline_prediction/prs_start/residsFixed_slope_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'fname', 'nfolds', 'nreps', 'meanRMSE', 'sdRMSE')
res[which.min(res$meanRMSE),]
```

Our best hi result is using the dti data: 

```
> res[which.min(res$meanRMSE),]
    sx          model
537 hi blassoAveraged
                                                                           fname
537 /home/sudregp/data/baseline_prediction/prs_start/gf_impute_based_dti_165.csv
    nfolds nreps meanRMSE sdRMSE
537     10    10 0.459838     NA
```

For inatt it's a different on, but still dti:

```
> res2[which.min(res2$meanRMSE),]
       sx      model
614 inatt blackboost
                                                                           fname
614 /home/sudregp/data/baseline_prediction/prs_start/gf_impute_based_dti_165.csv
    nfolds nreps meanRMSE sdRMSE
614     10    10 0.562351      0
```

Can we check combined, both for anat and DTI datasets?

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsFixed_slope_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'fname', 'nfolds', 'nreps', 'meanRMSE', 'sdRMSE')
for (reg in unique(res$model)) {
    for (nf in unique(res$nfolds)) {
        for (nr in unique(res$nreps)) {
            for (fn in unique(res$fname)) {
                idx = (res$model == reg &
                        res$fname == fn & res$nfolds == nf &
                        res$nreps == nr)
                pos = which(idx)
                if (length(pos) == 2) {
                    my_str = paste(c(reg, fn, nf, nr), collapse='_')
                    params = c(params, my_str)
                    scores = c(scores, mean(res[pos, 'meanRMSE']))
                }
            }
        }
    }
}
a = sort(scores, decreasing=F, index.return=T)
print(params[a$ix[1]])
```

So, for DTI blackboost did best:

```
> res[res$model=='blackboost' & res$nfolds==10 & res$nreps==10, c(1, 2, 4:6)]
       sx      model nfolds nreps meanRMSE
612 inatt blackboost     10    10 0.627078
613    hi blackboost     10    10 0.527887
614 inatt blackboost     10    10 0.562351
615    hi blackboost     10    10 0.460291
```

The top two are anat, but our best result is just dti. What's the best anat
result?

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsFixed_slope_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'fname', 'nfolds', 'nreps', 'meanRMSE', 'sdRMSE')
for (reg in unique(res$model)) {
    for (nf in unique(res$nfolds)) {
        for (nr in unique(res$nreps)) {
            idx = (res$model == reg & grepl(res$fname, pattern='anatomy') &
                   res$nfolds == nf & res$nreps == nr)
            pos = which(idx)
            if (length(pos) == 2) {
                my_str = paste(c(reg, nf, nr), collapse='_')
                params = c(params, my_str)
                scores = c(scores, mean(res[pos, 'meanRMSE']))
            }
        }
    }
}
a = sort(scores, decreasing=F, index.return=T)
print(params[a$ix[1]])
```

Then we're looking at conditional forest at 10x10:

```
> res[res$model=='cforest' & res$nfolds==10 & res$nreps==10, c(1:2, 4:7)]
       sx   model nfolds nreps meanRMSE   sdRMSE
414 inatt cforest     10    10 0.575755 0.006415
416    hi cforest     10    10 0.468978 0.003739
541 inatt cforest     10    10 0.619933 0.000505
545    hi cforest     10    10 0.525780 0.000603
```

Result is the bottom 2 again. How do they change if we lok at R2?

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsR2_slope_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'fname', 'nfolds', 'nreps', 'meanRsquared', 'sdRsquared')
for (reg in unique(res$model)) {
    for (nf in unique(res$nfolds)) {
        for (nr in unique(res$nreps)) {
            for (fn in unique(res$fname)) {
                idx = (res$model == reg &
                        res$fname == fn & res$nfolds == nf &
                        res$nreps == nr)
                pos = which(idx)
                if (length(pos) == 2) {
                    my_str = paste(c(reg, fn, nf, nr), collapse='_')
                    params = c(params, my_str)
                    scores = c(scores, mean(res[pos, 'meanRsquared']))
                }
            }
        }
    }
}
a = sort(scores, decreasing=T, index.return=T)
print(params[a$ix[1]])
```

So, for DTI kernelpls did best:

```
> res[res$model=='kernelpls' & res$nfolds==10 & res$nreps==10,][c(2,4),c(1,2,4:7)]
       sx     model nfolds nreps meanRsquared sdRsquared
312 inatt kernelpls     10    10     0.081295   0.012537
329    hi kernelpls     10    10     0.094364   0.011303
```

What's the best anat result?

```r
params = c()
scores = c()
res = read.csv('~/data/baseline_prediction/prs_start/residsR2_slope_impInter.csv', header=F)
colnames(res) = c('sx', 'model', 'fname', 'nfolds', 'nreps', 'meanRsquared', 'sdRsquared')
for (reg in unique(res$model)) {
    for (nf in unique(res$nfolds)) {
        for (nr in unique(res$nreps)) {
            idx = (res$model == reg & grepl(res$fname, pattern='anatomy') &
                   res$nfolds == nf & res$nreps == nr)
            pos = which(idx)
            if (length(pos) == 2) {
                my_str = paste(c(reg, nf, nr), collapse='_')
                params = c(params, my_str)
                scores = c(scores, mean(res[pos, 'meanRsquared']))
            }
        }
    }
}
a = sort(scores, decreasing=T, index.return=T)
print(params[a$ix[1]])
```

Then we're looking at bagEarth at 10x10:

```
> res[res$model=='bagEarth' & res$nfolds==10 & res$nreps==10, ][c(3,4),c(1,2,4:7)]
       sx    model nfolds nreps meanRsquared sdRsquared
571    hi bagEarth     10    10     0.093050   0.004943
573 inatt bagEarth     10    10     0.060669   0.011501
```

But how do all of these compare with a simple prediction of the mean? Can we
also check the feature weights?




