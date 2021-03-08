args = commandArgs(trailingOnly=TRUE)
if (length(args) != 6) {
  stop("Arguments: seed[-1] cpuDiff tuneLength mymod root_fname sx_target", call.=FALSE)
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
  y_target = args[6]
}

###########

fname = sprintf('%s_%04d.log', root_fname, myseed)
sink(fname, append=FALSE, split=TRUE)

source('~/ncr_notebooks/baseline_prediction/src/aux_functions.R')

tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
rm_me = (tract_data$fa_avg < .4 | tract_data$ad_avg < 1.18 | tract_data$rd_avg > .65 | tract_data$rd_avg < .5 |
           tract_data$norm.trans > .45 | tract_data$norm.rot > .008 | tract_data$goodSlices < 45 |
           tract_data$goodSlices > 70)
tract_data = tract_data[!rm_me, ]
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)

phen_vars = c(which(grepl("^FA_", colnames(tract_data))),
              which(grepl("^AD_", colnames(tract_data))),
              which(grepl("^RD_", colnames(tract_data))),
              which(grepl("^MO_", colnames(tract_data)))
)
Xraw = tract_data[, phen_vars]
my_ids = unique(gf$MRN)

# MAKE SURE nrow(tract_data) == nrow(Xraw)!!!!

get_delta = function (d) {
  # if we have too many NAs, return NA
  if (sum(is.na(d)) >= (length(d)-1)) {
    return(NA)
  }
  else {
    lm(d ~ tract_data[idx, ]$age_at_scan)$coefficients[2]
  }
}

X = c()
y = c()
target_col = which(colnames(gf)==y_target)
for (s in my_ids) {
  idx = tract_data$MRN==s
  # proceed if we have more than one observation in the data
  if (sum(idx) >= 2) {
    slopes = sapply(Xraw[idx, ], get_delta)
    names(slopes) = colnames(Xraw)
    X = rbind(X, slopes)
    idxy = gf$MRN==s
    y = c(y, unique(as.character(gf[idxy, target_col])))
  }
}
y = factor(y)

source('~/ncr_notebooks/baseline_prediction/src/do_multi_classification.R')

sink()
