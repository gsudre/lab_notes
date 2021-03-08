# let's cram all regressions in there
library(nlme)
gf = read.csv('~/data/prs/clinical_06192017.csv')
gf = gf[gf$ADHD_current_yes_no!='exclude',]
pgc = read.csv('~/data/prs/PRS2017_original_clump_default.csv')

df = merge(gf, pgc)
# remove duplicated MRNs
df = df[!duplicated(df$MRN),]

pheno = read.csv('~/data/prs/neuropsych_07282017.csv')
mydata = merge(df, pheno, by='MRN')

indep_vars = which(grepl('^PROFILES', colnames(mydata)))
dep_var_names = colnames(pheno)[2:28]
all_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  y_res = c()
  for (y_str in dep_var_names) {
    fm_str = sprintf('%s ~ %s', y_str, x_str) # make sure the interesting variable comes last!
    fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
    last_row = nrow(summary(fit)$tTable)
    res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
    y_res = rbind(y_res, res)
  }
  all_res = cbind(all_res, y_res)
}

cnames = c()
for (i in indep_vars) {
  cnames = c(cnames, c(sprintf('Tstat_%s', colnames(mydata)[i]),
                       sprintf('pval_%s', colnames(mydata)[i])))
}
colnames(all_res) = cnames
rownames(all_res) = dep_var_names


# now we just run the exact same thing for the brain variables
pheno = read.csv('~/data/prs/dti_07062017.csv')
mydata = merge(df, pheno, by='MRN')
rm_me = (mydata$numVolumes < 60 | mydata$norm.rot > .003 | mydata$norm.trans > .3 |
           mydata$mean_fa < .3 | mydata$mean_ad < .97 | mydata$mean_rd < .5)
mydata = mydata[!rm_me, ]
dep_var_names = colnames(pheno)[22:68]
pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  y_res = c()
  for (y_str in dep_var_names) {
    fm_str = sprintf('%s ~ %s', y_str, x_str) # make sure the interesting variable comes last!
    fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
    last_row = nrow(summary(fit)$tTable)
    res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
    y_res = rbind(y_res, res)
  }
  pheno_res = cbind(pheno_res, y_res)
}
cnames = c()
for (i in indep_vars) {
  cnames = c(cnames, c(sprintf('Tstat_%s', colnames(mydata)[i]),
                       sprintf('pval_%s', colnames(mydata)[i])))
}
colnames(pheno_res) = cnames
rownames(pheno_res) = dep_var_names
all_res = rbind(all_res, pheno_res)

pheno = read.csv('~/data/prs/struct_08042017.csv')
pheno = pheno[pheno$scanner != '1p5t',]
mydata = merge(df, pheno, by='MRN')
keep_me = which(mydata$avg_freesurfer_score <= 2 & mydata$MPRAGE_QC <= 2)
# in the end, all data needs to be in a matrix called mydata!
mydata = mydata[keep_me, ]
dep_var_names = colnames(pheno)[11:280]
pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  y_res = c()
  for (y_str in dep_var_names) {
    fm_str = sprintf('%s ~ %s', y_str, x_str) # make sure the interesting variable comes last!
    fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
    if (length(fit) > 1) {
      last_row = nrow(summary(fit)$tTable)
      res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
    } else {
      res = c(NA, NA)
    }
    y_res = rbind(y_res, res)
  }
  pheno_res = cbind(pheno_res, y_res)
}
cnames = c()
for (i in indep_vars) {
  cnames = c(cnames, c(sprintf('Tstat_%s', colnames(mydata)[i]),
                       sprintf('pval_%s', colnames(mydata)[i])))
}
colnames(pheno_res) = cnames
rownames(pheno_res) = dep_var_names
all_res = rbind(all_res, pheno_res)


# finally, we try different clinical metrics
mydata = df
clin_res = c()
pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  fm_str = sprintf('%s ~ ADHD_current_yes_no', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  fm_str = sprintf('SX_inatt ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  fm_str = sprintf('SX_HI ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

mydata = df[df$ADHD_current_yes_no!='no',]
pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  fm_str = sprintf('SX_inatt ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  fm_str = sprintf('SX_HI ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  y = pmax(mydata[, 'SX_inatt'], mydata[, 'inatt_child'], na.rm=T)
  fm_str = sprintf('y ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

pheno_res = c()
for (i in 1:length(indep_vars)) {
  x_str = colnames(mydata)[indep_vars][i]
  y = pmax(mydata[, 'SX_HI'], mydata[, 'HI_child'], na.rm=T)
  fm_str = sprintf('y ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = c(pheno_res, res)
}
clin_res = rbind(clin_res, pheno_res)

rownames(clin_res) = c('DX', 'inatt', 'HI', 'inattADHDonly',
                       'HIADHDOnly', 'maxInatt', 'maxHI')
colnames(clin_res) = cnames
all_res = rbind(all_res, pheno_res)

write.csv(all_res, file='~/data/prs/results/all_model6_fromX_no1p5t.csv')


#####################
# Now we need to work on the M1 to M2 interactions
#####################
pheno = read.csv('~/data/prs/dti_07062017.csv')
M1s = colnames(pheno)[22:68]
mydata = merge(df, pheno, by='MRN')
rm_me = (mydata$numVolumes < 60 | mydata$norm.rot > .003 | mydata$norm.trans > .3 |
           mydata$mean_fa < .3 | mydata$mean_ad < .97 | mydata$mean_rd < .5)
mydata = mydata[!rm_me, ]

pheno = read.csv('~/data/prs/neuropsych_07282017.csv')
M2s = colnames(pheno)[2:28]
mydata = merge(mydata, pheno, by='MRN')

all_res = c()
for (x_str in M1s) {
  y_res = c()
  for (y_str in M2s) {
    fm_str = sprintf('%s ~ %s', y_str, x_str) # make sure the interesting variable comes last!
    fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
    if (length(fit) > 1) {
      last_row = nrow(summary(fit)$tTable)
      res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
    } else {
      res = c(NA, NA)
    }
    y_res = c(y_res, res)
  }
  all_res = rbind(all_res, y_res)
}
rownames(all_res) = M1s


pheno = read.csv('~/data/prs/struct_08042017.csv')
pheno = pheno[pheno$scanner != '1p5t',]
mydata = merge(df, pheno, by='MRN')
keep_me = which(mydata$avg_freesurfer_score <= 2 & mydata$MPRAGE_QC <= 2)
mydata = mydata[keep_me, ]
M1s = colnames(pheno)[11:280]
pheno = read.csv('~/data/prs/neuropsych_07282017.csv')
mydata = merge(mydata, pheno, by='MRN')
pheno_res = c()
for (x_str in M1s) {
  y_res = c()
  for (y_str in M2s) {
    fm_str = sprintf('%s ~ %s', y_str, x_str) # make sure the interesting variable comes last!
    fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
    if (length(fit) > 1) {
      last_row = nrow(summary(fit)$tTable)
      res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
    } else {
      res = c(NA, NA)
    }
    y_res = c(y_res, res)
  }
  pheno_res = rbind(pheno_res, y_res)
}
rownames(pheno_res) = M1s
all_res = rbind(all_res, pheno_res)
cnames = c()
for (i in M2s) {
  cnames = c(cnames, c(sprintf('Tstat_%s', i), sprintf('pval_%s', i)))
}
colnames(all_res) = cnames

write.csv(all_res, file='~/data/prs/results/all_model6_M1toM2_no1p5t.csv')

################
# Finally, let's run everything going to Y
################
pheno = read.csv('~/data/prs/dti_07062017.csv')
M1s = colnames(pheno)[22:68]
mydata = merge(df, pheno, by='MRN')
rm_me = (mydata$numVolumes < 60 | mydata$norm.rot > .003 | mydata$norm.trans > .3 |
           mydata$mean_fa < .3 | mydata$mean_ad < .97 | mydata$mean_rd < .5)
mydata = mydata[!rm_me, ]

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('%s ~ ADHD_current_yes_no', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  if (length(fit) > 1) {
    last_row = nrow(summary(fit)$tTable)
    res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  } else {
    res = c(NA, NA)
  }
  pheno_res = rbind(pheno_res, res)
}
all_res = pheno_res

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_inatt ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res = cbind(all_res, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_HI ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res = cbind(all_res, pheno_res)

mydata = mydata[df$ADHD_current_yes_no!='no',]
pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_inatt ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res = cbind(all_res, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_HI ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res = cbind(all_res, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  y = pmax(mydata[, 'SX_inatt'], mydata[, 'inatt_child'], na.rm=T)
  fm_str = sprintf('y ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res = cbind(all_res, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  y = pmax(mydata[, 'SX_HI'], mydata[, 'HI_child'], na.rm=T)
  fm_str = sprintf('y ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res = cbind(all_res, pheno_res)

rownames(all_res) = M1s
cols = c('DX', 'inatt', 'HI', 'inattADHDonly',
         'HIADHDOnly', 'maxInatt', 'maxHI')
cnames = c()
for (i in cols) {
  cnames = c(cnames, c(sprintf('Tstat_%s', i), sprintf('pval_%s', i)))
}
colnames(all_res) = cnames

pheno = read.csv('~/data/prs/struct_08042017.csv')
pheno = pheno[pheno$scanner != '1p5t',]
mydata = merge(df, pheno, by='MRN')
keep_me = which(mydata$avg_freesurfer_score <= 2 & mydata$MPRAGE_QC <= 2)
mydata = mydata[keep_me, ]
M1s = colnames(pheno)[11:280]

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('%s ~ ADHD_current_yes_no', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  if (length(fit) > 1) {
    last_row = nrow(summary(fit)$tTable)
    res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  } else {
    res = c(NA, NA)
  }
  pheno_res = rbind(pheno_res, res)
}
all_res2 = pheno_res

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_inatt ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res2 = cbind(all_res2, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_HI ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res2 = cbind(all_res2, pheno_res)

mydata = mydata[df$ADHD_current_yes_no!='no',]
pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_inatt ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res2 = cbind(all_res2, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  fm_str = sprintf('SX_HI ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res2 = cbind(all_res2, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  y = pmax(mydata[, 'SX_inatt'], mydata[, 'inatt_child'], na.rm=T)
  fm_str = sprintf('y ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res2 = cbind(all_res2, pheno_res)

pheno_res = c()
for (x_str in M1s) {
  y = pmax(mydata[, 'SX_HI'], mydata[, 'HI_child'], na.rm=T)
  fm_str = sprintf('y ~ %s', x_str) # make sure the interesting variable comes last!
  fit = try(lme(as.formula(fm_str), random=~1|NuclearFamID, data=mydata, na.action = na.omit))
  last_row = nrow(summary(fit)$tTable)
  res = c(summary(fit)$tTable[last_row,'t-value'], summary(fit)$tTable[last_row,'p-value'])
  pheno_res = rbind(pheno_res, res)
}
all_res2 = cbind(all_res2, pheno_res)

rownames(all_res2) = M1s
cols = c('DX', 'inatt', 'HI', 'inattADHDonly',
         'HIADHDOnly', 'maxInatt', 'maxHI')
cnames = c()
for (i in cols) {
  cnames = c(cnames, c(sprintf('Tstat_%s', i), sprintf('pval_%s', i)))
}
colnames(all_res2) = cnames

all_res = rbind(all_res, all_res2)

write.csv(all_res, file='~/data/prs/results/all_model6_toY_no1p5t.csv')

