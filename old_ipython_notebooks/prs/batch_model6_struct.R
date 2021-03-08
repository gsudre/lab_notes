# loading clinical data
gf = read.csv('~/data/prs/clinical_06192017.csv')
gf = gf[gf$ADHD_current_yes_no!='exclude',]
gf$ADHD_current_yes_no = factor(gf$ADHD_current_yes_no)

# loading PRS data
pgc = read.csv('~/data/prs/PRS2017_original_clump_default.csv')
df = merge(gf, pgc)
# remove duplicated MRNs
df = df[!duplicated(df$MRN),]

# loading neuropsych
neuropsych = read.csv('~/data/prs/neuropsych_07282017.csv')
df = merge(df, neuropsych, by='MRN')

# loading brain structural
struct = read.csv('~/data/prs/struct_08042017.csv')
rois = merge(df, struct, by='MRN')
# filtering on QC
keep_me = which(rois$avg_freesurfer_score <= 2 & rois$MPRAGE_QC <= 2)
# in the end, all data needs to be in a matrix called mydata!
mydata = rois[keep_me, ]

# choosing mediators
args = commandArgs(trailingOnly=TRUE)
m1_name = args[1]
M2s = c(33:44, 47:59)

out_fname = sprintf('~/data/prs/results/struct/model6_p3_%s_neuropsych_DX_QCse2Both.csv', m1_name)
X = mydata$PROFILES.0.3.profile
Y = mydata$ADHD_current_yes_no
nboot = 1000
ncpus = 4


run_model6 = function(X, M1, M2, Y, nboot=1000) {
  library(lavaan)
  idx = is.na(X) | is.na(Y) | is.na(M1) | is.na(M2)
  Y = Y[!idx]
  if (!is.factor(Y)) {
    Y = scale(Y)
  }
  run_data = data.frame(X = scale(X[!idx]),
                        Y = Y,
                        M1 = scale(M1[!idx]),
                        M2 = scale(M2[!idx]))
  
  hayes6 <- ' Y ~ b1*M1
  Y ~ b2*M2
  Y ~ cdash*X
  # direct effect of X on Y
  directXonY := cdash
  
  M1 ~ a1*X
  M2 ~ a2*X
  M2 ~ d1*M1
  
  # Specific indirect effect of X on Y via M1
  indXonYvM1 := a1*b1
  # Specific indirect effect of X on Y via M2
  indXonYvM2 := a2*b2
  # Specific indirect effect of X on Y via M1 and M2
  indXonYvM1M2 := a1*d1*b2
  # Total indirect effect of X on Y via M1, M2
  totalIndXonYvM1M2 := a1*b1 + a2*b2 + a1*d1*b2
  # Total effect of X on Y
  totalXonY := a1*b1 + a2*b2 + a1*d1*b2 + cdash'
  
  # fit model
  sem6 <- sem(model = hayes6,
              data = run_data,
              se = "bootstrap",
              bootstrap = nboot, parallel='multicore', ncpu=4)
  
  res = parameterEstimates(sem6,
                           boot.ci.type = "bca.simple",
                           level = .95, ci = TRUE,
                           standardized = FALSE)
  res2 = res[11:16, 5:10]
  rownames(res2) = res[11:16,]$label
  # adding n
  res2 = cbind(rep(nrow(run_data), nrow(res2)), res2)
  colnames(res2)[1] = 'n'
  return(res2)
}

run_wrapper = function(m2, run_model, mydata, m1_name, nboot, X, Y, m1) {
  m2_name = colnames(mydata)[m2]
  cat('\t', sprintf('M2=%s', m2_name), '\n')
  tmp = run_model(X, mydata[, m1], mydata[, m2], Y, nboot=nboot)
  res = c(m1_name, m2_name)
  names(res) = c('M1', 'M2')
  for (i in 1:nrow(tmp)) {
    for (j in 1:ncol(tmp)) {
      res = c(res, tmp[i, j])
      names(res)[length(res)] = sprintf('%s_%s', rownames(tmp)[i], colnames(tmp)[j])
    }
  }
  return(res)
}

# no need to change anything below here. The functions remove NAs and zscore variables on their own
library(parallel)
cl <- makeCluster(ncpus)

m1 = which(colnames(mydata)==m1_name)
m1_res = parLapply(cl, M2s, run_wrapper, run_model6, mydata, m1_name, nboot, X, Y, m1)
# m1_res = lapply(M2s, run_wrapper, run_model6, mydata, m1_name, nboot, X, Y, m1)
all_res = do.call(rbind, m1_res)
  
stopCluster(cl)
write.csv(all_res, file=out_fname, row.names=F)
