args = commandArgs(trailingOnly=TRUE)
if (length(args) == 0) {
    myseed = 1234
    cpuDiff = 1
    tuneLength = 1
    mymod = 'AdaBag' # AdaBoost.M1, AdaBag, ada
    root_fname = '~/tmp/TMP'
} else if (length(args) != 5) {
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

merged = merge(mwj, mbeery, by='MRN', all.x = T, all.y = T)
merged = merge(merged, miq, by='MRN', all.x = T, all.y = T)
merged = merge(merged, mwisc, by='MRN', all.x = T, all.y = T)
merged = merge(merged, mcpt, by='MRN', all.x = T, all.y = T)

tract_data = read.csv('~/data/baseline_prediction/stripped/dti.csv')
rm_me = (tract_data$fa_avg < .4 | tract_data$ad_avg < 1.18 | tract_data$rd_avg > .65 | tract_data$rd_avg < .5 |
tract_data$norm.trans > .45 | tract_data$norm.rot > .008 | tract_data$goodSlices < 45 | tract_data$goodSlices > 70)
print(sprintf('Reducing from %d to %d scans', nrow(tract_data), nrow(tract_data)-sum(rm_me)))
tract_data = tract_data[!rm_me, ]
my_ids = intersect(gf$MRN, tract_data$MRN)
mdti = mergeOnClosestDate(gf, tract_data, my_ids)
rm_me = abs(mdti$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d scans', nrow(mdti), nrow(mdti)-sum(rm_me)))
mdti = mdti[!rm_me, ]

geo_data = read.csv('~/data/baseline_prediction/stripped/geospatial.csv')
mgeo = merge(gf, geo_data, by='MRN')
# some variables are being read as numeric...
mgeo$Home_Price = as.numeric(mgeo$Home_Price)
mgeo$Fam_Income = as.numeric(mgeo$Fam_Income)
mgeo$Crime_Rate = as.numeric(mgeo$Crime_Rate)

merged = merge(merged, mdti, by='MRN', all.x = T, all.y = T)
merged = merge(merged, mgeo, by='MRN', all.x = T, all.y = T)

struct_data = read.csv('~/data/baseline_prediction/stripped/structural.csv')
rm_me = (struct_data$mprage_score > 2)
struct_data = struct_data[!rm_me, ]
my_ids = intersect(gf$MRN, struct_data$MRN)
mstruct = mergeOnClosestDate(gf, struct_data, my_ids)
rm_me = abs(mstruct$dateX.minus.dateY.months) > 12
print(sprintf('Reducing from %d to %d scans', nrow(mstruct), nrow(mstruct)-sum(rm_me)))
mstruct = mstruct[!rm_me, ]
merged = merge(merged, mstruct, by='MRN', all.x = T, all.y = T)

# putting back clinical data
clin = read.csv(gf_fname)
my_ids = gf$MRN
merged = mergeOnClosestDate(merged, clin, my_ids)

phen_vars = c('FSIQ',
              # CPT
              'N_of_omissions', 'N_commissions', 'hit_RT', 'hit_RT_SE', 'variability_of_SE', 'N_perservations',
              'hit_RT_block_change', 'hit_RT_SE_block_change', 'hit_RT_ISI_change', 'hit_RT_SE_ISI_change',
              # WISC
              'Raw.score..DSF', 'Raw.score..DSB', 'Raw.score..SSF', 'Raw.score..SSB',
              # WJ
              'PS',
              # Beery
              'Standard.score',
              #GeoSpatial
              'SES', 'Home_Type', 'Home_Price', 'Fam_Income', 'Pop_BPL', 'Fam_BPL', 'Pub_School',
              'Crime_Rate', 'Green_Space', 'Park_Access', 'Air_Quality', 'Obesity_Rate',
              'Food_Index', 'Exercise_Access', 'Excessive_Drinking',
              # DTI
              colnames(merged)[grepl("^FA_", colnames(merged))],
              colnames(merged)[grepl("^AD_", colnames(merged))],
              colnames(merged)[grepl("^RD_", colnames(merged))],
              colnames(merged)[grepl("^MO_", colnames(merged))],
              colnames(merged)[grepl("^lh_", colnames(merged))],
              colnames(merged)[grepl("^rh_", colnames(merged))],
              colnames(merged)[grepl("^Left", colnames(merged))],
              colnames(merged)[grepl("^Right", colnames(merged))],
              colnames(merged)[grepl("^CC_", colnames(merged))]
)
keep_me = sapply(phen_vars, function(d) which(colnames(merged) == d))
X = merged[, keep_me]
y = as.character(merged$DX_BASELINE)
y[y != 'NV'] = 'ADHD'
y[y == 'ADHD' & merged$GAS > 68.91] = 'ADHD_high'
y[y == 'ADHD'] = 'ADHD_low'
y = factor(y, levels=c('NV', 'ADHD_low', 'ADHD_high'))

source('~/ncr_notebooks/baseline_prediction/src/do_multi_classification.R')

sink()
