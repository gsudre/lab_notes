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

beery_data = read.csv('~/data/baseline_prediction/stripped/beeryVMI.csv')
gf_fname = '~/data/baseline_prediction/stripped/clinical.csv'
gf = read.csv(gf_fname)
gf = gf[gf$BASELINE=='BASELINE', ]
# we only need to keep MRN and DOA for now, to avoid duplicated
gf = gf[, c('MRN', 'DOA')]
my_ids = intersect(gf$MRN, beery_data$Medical.Record...MRN)
mbeery = mergeOnClosestDate(gf, beery_data, my_ids, y.date='record.date.collected', y.id='Medical.Record...MRN')
rm_me = abs(mbeery$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mbeery), nrow(mbeery)-sum(rm_me)))
mbeery = mbeery[!rm_me, ]
mbeery$dateClinical.minus.dateBeery.months = mbeery$dateX.minus.dateY.months
mbeery$dateX.minus.dateY.months = NULL

cpt_data = read.csv('~/data/baseline_prediction/stripped/cpt.csv')
my_ids = intersect(gf$MRN, cpt_data$MRN)
mcpt = mergeOnClosestDate(gf, cpt_data, my_ids)
rm_me = abs(mcpt$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mcpt), nrow(mcpt)-sum(rm_me)))
mcpt = mcpt[!rm_me, ]
mcpt$dateClinical.minus.dateCPT.months = mcpt$dateX.minus.dateY.months
mcpt$dateX.minus.dateY.months = NULL
colnames(mcpt)[colnames(mcpt)=="DOA.x"] <- "DOA.gf"
colnames(mcpt)[colnames(mcpt)=="DOA.y"] <- "DOA.cpt"

iq_data = read.csv('~/data/baseline_prediction/stripped/iq.csv')
my_ids = intersect(gf$MRN, iq_data$Medical.Record...MRN)
miq = mergeOnClosestDate(gf, iq_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(miq$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(miq), nrow(miq)-sum(rm_me)))
miq = miq[!rm_me, ]
miq$dateClinical.minus.dateIQ.months = miq$dateX.minus.dateY.months
miq$dateX.minus.dateY.months = NULL

wisc_data = read.csv('~/data/baseline_prediction/stripped/wisc.csv')
my_ids = intersect(gf$MRN, wisc_data$Medical.Record...MRN)
mwisc = mergeOnClosestDate(gf, wisc_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(mwisc$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mwisc), nrow(mwisc)-sum(rm_me)))
mwisc = mwisc[!rm_me, ]
mwisc$dateClinical.minus.dateWISC.months = mwisc$dateX.minus.dateY.months
mwisc$dateX.minus.dateY.months = NULL

wj_data = read.csv('~/data/baseline_prediction/stripped/wj.csv')
my_ids = intersect(gf$MRN, wj_data$Medical.Record...MRN)
mwj = mergeOnClosestDate(gf, wj_data, my_ids, y.id='Medical.Record...MRN', y.date='record.date.collected')
rm_me = abs(mwj$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d tests', nrow(mwj), nrow(mwj)-sum(rm_me)))
mwj = mwj[!rm_me, ]
mwj$dateClinical.minus.dateWJ.months = mwj$dateX.minus.dateY.months
mwj$dateX.minus.dateY.months = NULL

merged = merge(mwj, mbeery, by='MRN', all.x = T, all.y = T, suffixes=c('.wj', '.beery'))
merged = merge(merged, miq, by='MRN', all.x = T, all.y = T)
merged = merge(merged, mwisc, by='MRN', all.x = T, all.y = T, suffixes=c('.iq', '.wisc'))
merged = merge(merged, mcpt, by='MRN', all.x = T, all.y = T)

clin = read.csv(gf_fname)
clin = clin[clin$BASELINE=='BASELINE', ]
my_ids = gf$MRN
merged = merge(merged, clin, by='MRN')

phen_vars = c('FSIQ',
              # CPT
              'N_of_omissions', 'N_commissions', 'hit_RT', 'hit_RT_SE', 'variability_of_SE', 'N_perservations',
              'hit_RT_block_change', 'hit_RT_SE_block_change', 'hit_RT_ISI_change', 'hit_RT_SE_ISI_change',
              # WISC
              'Raw.score..DSF', 'Raw.score..DSB', 'Raw.score..SSF', 'Raw.score..SSB',
              # WJ
              'PS',
              # Beery
              'Standard.score')
keep_me = c()
for (v in phen_vars) {
  keep_me = c(keep_me, which(colnames(merged) == v))
}
X = merged[, keep_me]
y = merged$DX_BASELINE
y[y != 'NV'] = 'ADHD'
y = factor(y, levels = c('NV', 'ADHD'))

source('~/ncr_notebooks/baseline_prediction/src/do_binary_classification.R')

sink()
