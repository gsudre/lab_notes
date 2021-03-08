# loading clinical data
gf = read.csv('~/data/prs/clinical_08232017.csv')

idx = gf$Race=='[White]' & gf$Ethnicity=='Not Hispanic or Latino'
gf = gf[idx, ]

# loading PRS data
pgc = read.csv('~/data/prs/PRS2017_noInversion_eur_all.csv')
df = merge(gf, pgc)
# remove duplicated MRNs
df = df[!duplicated(df$MRN),]

# loading mediator data
pheno = read.csv('~/data/prs/struct_5vars_08232017.csv')

merged = merge(df, pheno, by='MRN')
# filtering on QC
rm_me = (merged$avg_freesurfer_score > 2 & merged$MPRAGE_QC > 2)
# rm_me = (merged$numVolumes < 60 | merged$norm.rot > .003 | merged$norm.trans > .3 |
#          merged$mean_fa < .3 | merged$mean_ad < .97 | merged$mean_rd < .5)
# in the end, all data needs to be in a matrix called mydata!
mydata = merged[!rm_me, ]

# choosing mediators
Ms = c(57)#:60)
# Ms = c(47:50, 56)
Xs = c('PROFILES.0.05.profile')#, 'PROFILES.0.1.profile', 'PROFILES.0.2.profile',
       # 'PROFILES.0.3.profile', 'PROFILES.0.4.profile', 'PROFILES.0.5.profile')
Ys = c('SX_inatt')#, 'SX_HI', 'SX_total')
mydata$SX_total = mydata$SX_inatt + mydata$SX_HI

nboot = 1000
ncpus = 1
mixed = T
fname_root = '~/data/prs/results/TEST_struct_%s_lme_sex_%s_eur_5K.csv'

# no need to change anything below here. The functions remove NAs and zscore variables on their own
run_model4 = function(X, M, Y, nboot=1000, short=T, data2) {
  library(mediation)
  idx = is.na(X) | is.na(Y) | is.na(M)
  Y = Y[!idx]
  imdiscrete = T
  if (!is.factor(Y)) {
    Y = scale(Y)
    imdiscrete = F
  }
  run_data = data.frame(X = scale(X[!idx]),
                        Y = Y,
                        M = scale(M[!idx]),
                        FAMID = data2[!idx,]$NuclearFamID,
                        age= data2[!idx,]$AGE,
                        sex = data2[!idx,]$Sex)
  
  if (!is.na(run_data[1,]$FAMID)) {
    library(lme4)
    fm = as.formula('M ~ X + sex + (1|FAMID)')
    fy = as.formula('Y ~ X + M + sex + (1|FAMID)')
    model.M <- lmer(fm, data=run_data)
    if (imdiscrete) {
      model.Y <- glmer(fy, data=run_data, family=binomial(link='logit'))
    } else {
      model.Y <- lmer(fy, data=run_data)
    }
    results <- mediate(model.M, model.Y, treat='X', mediator='M', boot=F, sims=nboot)
    print(summary(results))
  } else {
    fm = as.formula('M ~ X')
    fy = as.formula('Y ~ X + M')
    model.M <- lm(fm, data=run_data)
    if (imdiscrete) {
      model.Y <- glm(fy, data=run_data, family=binomial(link='logit'))
    } else {
      model.Y <- lm(fy, data=run_data)
    }
    results <- mediate(model.M, model.Y, treat='X', mediator='M', boot=T, sims=nboot, boot.ci.type='bca')
  }
  
  if (short) {
    res = c(results$mediator, results$nobs, results$tau.coef, results$tau.p, results$d.avg, results$d.avg.p,
            results$z.avg, results$z.avg.p, results$n.avg, results$n.avg.p)
    names(res) = c('M', 'nobs', 'tot', 'tot_p', 'acme', 'acme_p', 'ade', 'ade_p', 'prop', 'prop_p')
    return(res)
  } else {
    return(results)
  }
}

run_wrapper = function(m, run_model, mydata, nboot, X, Y) {
  cat('\t', sprintf('M=%s', colnames(mydata)[m]), '\n')
  tmp = run_model(X, mydata[, m], Y, nboot=nboot, data2=mydata)
  tmp$M = colnames(mydata)[m]
  return(tmp)
}

if (!mixed) {
  mydata$NuclearFamID = NA
}

for (x_str in Xs) {
  for (y_str in Ys) {
    X = mydata[, x_str]
    Y = mydata[, y_str]
    out_fname = sprintf(fname_root,
                        y_str, x_str)
    print(out_fname)
    
    if (ncpus > 1) {
      library(parallel)
      cl <- makeCluster(ncpus)
      m1_res = parLapply(cl, Ms, run_wrapper, run_model4, mydata, nboot, X, Y)
      stopCluster(cl)
    } else {
      m1_res = lapply(Ms, run_wrapper, run_model4, mydata, nboot, X, Y)
    }
    all_res = do.call(rbind, m1_res)
    
    write.csv(all_res, file=out_fname, row.names=F)
  }
}


