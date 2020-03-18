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