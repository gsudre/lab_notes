# 2019-01-23 15:04:55

The more I played with the logistic regression results, the less I liked them.
The comparison to NV doesn't make that much sense to me.

So, let's take a couple hours to go back to the prediction framework, this time
using a few simpler algorithms.

```r
if (Sys.info()['sysname'] == 'Darwin') {
  max_mem = '16G'
  base_name = '~/data/'
} else {
  max_mem = paste(Sys.getenv('SLURM_MEM_PER_NODE'),'m',sep='')
  base_name = '/data/NCR_SBRB/'
}
clin = read.csv(sprintf('%s/baseline_prediction/long_clin_11302018.csv', base_name))
load(sprintf('%s/baseline_prediction/struct_thickness_11142018_260timeDiff12mo.RData.gz', base_name))
df = merge(clin, data, by='MRN')
df$OLS_inatt_categ = NULL
df[df$OLS_inatt_slope <= -.33, 'OLS_inatt_categ'] = 'marked'
df[df$OLS_inatt_slope > -.33 & df$OLS_inatt_slope <= 0, 'OLS_inatt_categ'] = 'mild'
df[df$OLS_inatt_slope > 0, 'OLS_inatt_categ'] = 'deter'
df = df[df$DX != 'NV', ]
df$OLS_inatt_categ = as.factor(df$OLS_inatt_categ)
x = colnames(df)[grepl(pattern = '^v_rh', colnames(df))]

library(caret)
set.seed(42)
control <- trainControl(method="repeatedcv", number=10, repeats=3)
# train the models
set.seed(42)
modelLinearSvm <- train(df[, x], df$OLS_inatt_categ, method="svmLinear", trControl=control)
set.seed(42)
modelNB <- train(df[, x], df$OLS_inatt_categ, method="nb", trControl=control)


```

Philip doens't want to go this route yet.