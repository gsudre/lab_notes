run_model4 = function(X, M, Y, nboot=1000) {
  library(lavaan)
  idx = is.na(X) | is.na(Y) | is.na(M)
  run_data = data.frame(X = scale(X[!idx]),
                        Y = scale(Y[!idx]),
                        M = scale(M[!idx]))
  
  hayes4 <- ' # direct effect
  Y ~ c*X
  direct := c
  
  # regressions
  M ~ a*X
  Y ~ b*M
  
  # indirect effect (a*b)
  indirect := a*b
  
  # total effect
  total := c + (a*b)'
  
  # fit model
  sem4 <- sem(model = hayes4,
              data = run_data,
              se = "bootstrap",
              bootstrap = nboot)
  
  res = parameterEstimates(sem4,
                           boot.ci.type = "bca.simple",
                           level = .95, ci = TRUE,
                           standardized = FALSE)
  res2 = res[7:9, 5:10]
  rownames(res2) = res[7:9,]$label
  return(res2)
}

run_model6 = function(X, M1, M2, Y, nboot=1000) {
  library(lavaan)
  idx = is.na(X) | is.na(Y) | is.na(M1) | is.na(M2)
  run_data = data.frame(X = scale(X[!idx]),
                        Y = scale(Y[!idx]),
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
              bootstrap = nboot)
  
  res = parameterEstimates(sem6,
                           boot.ci.type = "bca.simple",
                           level = .95, ci = TRUE,
                           standardized = FALSE)
  res2 = res[11:16, 5:10]
  rownames(res2) = res[11:16,]$label
  return(res2)
}


# CHANGES START HERE

# loading clinical data
gf = read.csv('~/data/prs/clinical_06192017.csv')
gf = gf[gf$ADHD_current_yes_no!='exclude',]

# loading PRS data
pgc = read.csv('~/data/prs/PRS2017_original_clump_default.csv')
df = merge(gf, pgc)
# remove duplicated MRNs
df = df[!duplicated(df$MRN),]

# loading neuropsych
neuropsych = read.csv('~/data/prs/neuropsych_07072017.csv')
mydata = merge(df, neuropsych, by='MRN')

# choosing mediators
Ms = c(33:39, 42)

out_fname = '~/data/prs/model4_inatt.csv'
X = mydata$PROFILES.0.3.profile
Y = mydata$SX_inatt
nboot = 100

# no need to change anything below here. The functions remove NAs and zscore variables on their own
all_res = c()
for (m in Ms) {
  m_name = colnames(mydata)[m]
  print(sprintf('Running %s', m_name))
  tmp = run_model4(X, mydata[, m], Y, nboot=nboot)
  res = c()
  for (i in 1:nrow(tmp)) {
    for (j in 1:ncol(tmp)) {
      res = c(res, tmp[i, j])
      names(res)[length(res)] = sprintf('%s_%s', rownames(tmp)[i], colnames(tmp)[j])
    }
  }
  all_res = rbind(all_res, res)
  rownames(all_res)[nrow(all_res)] = m_name
}
write.csv(all_res, file=out_fname)
