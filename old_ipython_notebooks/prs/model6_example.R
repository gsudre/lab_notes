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

# CHANGE STARTING HERE

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
df = merge(df, neuropsych, by='MRN')

# loading brain structural
struct = read.csv('~/data/prs/struct_07112017.csv')
rois = merge(df, struct, by='MRN')
# filtering on QC
keep_me = which(rois$avg_freesurfer_score < 2.5 & rois$MPRAGE_QC <= 2)
# in the end, all data needs to be in a matrix called mydata!
mydata = rois[keep_me, ]

# choosing mediators
M1s = which(grepl("lh_[a-z]*_thickness", colnames(mydata)))[1:10]
M2s = c(33:39, 42)

out_fname = '~/data/prs/model6_inatt.csv'
X = mydata$PROFILES.0.3.profile
Y = mydata$SX_inatt
nboot = 100

# no need to change anything below here. The functions remove NAs and zscore variables on their own
all_res = c()
for (m1 in M1s) {
  m1_name = colnames(mydata)[m1]
  for (m2 in M2s) {
    m2_name = colnames(mydata)[m2]
    print(sprintf('Running M1=%s, M2=%s', m1_name, m2_name))
    tmp = run_model6(X, mydata[, m1], mydata[, m2], Y, nboot=nboot)
    res = c(m1_name, m2_name)
    names(res) = c('M1', 'M2')
    for (i in 1:nrow(tmp)) {
      for (j in 1:ncol(tmp)) {
        res = c(res, tmp[i, j])
        names(res)[length(res)] = sprintf('%s_%s', rownames(tmp)[i], colnames(tmp)[j])
      }
    }
    all_res = rbind(all_res, res)
  }
}
write.csv(all_res, file=out_fname, row.names=F)
