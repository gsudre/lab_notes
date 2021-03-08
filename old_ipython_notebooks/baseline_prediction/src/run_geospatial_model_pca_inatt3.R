args = commandArgs(trailingOnly=TRUE)
if (length(args) != 5) {
  stop("Arguments: seed[-1] cpuDiff tuneLength mymod root_fname", call.=FALSE)
} else {
  if (args[1] == '-1') {
    myseed = as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31)
  } else {
    myseed = as.integer(args[1])
  }
  cpuDiff = as.integer(args[2])
  tuneLength = as.integer(args[3])
  mymod = args[4] # AdaBoost.M1, AdaBag, ada
  root_fname = args[5]
}

###########

fname = sprintf('%s_%04d.log', root_fname, myseed)
sink(fname, append=FALSE, split=TRUE)
source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')

geo_data = read.csv('~/data/baseline_prediction/stripped/geospatial.csv')
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf = gf[gf$BASELINE=='BASELINE', ]
merged = merge(gf, geo_data, by='MRN')
# some variables are being read as numeric...
merged$Home_Price = as.numeric(merged$Home_Price)
merged$Fam_Income = as.numeric(merged$Fam_Income)
merged$Crime_Rate = as.numeric(merged$Crime_Rate)

phen_vars = c('SES', 'poverty_education', 'social_environment', 'physical_health_environment', 
              'physical_envronment')
keep_me = c()
for (v in phen_vars) {
  keep_me = c(keep_me, which(colnames(merged) == v))
}
X = merged[, keep_me]
y = merged$inatt3_named
y = factor(y, levels=c('low', 'medium', 'high'))

source('~/ncr_notebooks/baseline_prediction/src/do_multi_classification.R')

sink()
