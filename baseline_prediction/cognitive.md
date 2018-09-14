# 2018-09-14 14:15:53

The idea here is that not everybody has all cognitive data, so let's split it into the different tests. We can run one with everything as well, but for now let's be conservative:

```r
df = read.csv('~/data/baseline_prediction/cog_baseline_merge_clean.csv')
x = c("N_of_omissions", "PC_omissions", "N_commissions", "PC_commissions",
      "hit_RT", "hit_RT_SE", "variability_of_SE", "D.prime", "beta",
      "N_perservations", "PC_perservations", "hit_RT_block_change",
      "hit_RT_SE_block_change", "hit_RT_ISI_change", "hit_RT_SE_ISI_change",
      "RAW_DS_total", "RAW_DSB", "Raw_DSF", "Raw_SS_total", "Raw_SSB",
      "Raw_SSF", "Standard_DS_total", "Standard_DSB", "Standard_DSF",
      "Standard_SS_total", "Standard_SSB", "Standard_SSB.1", "PS", "SSDS_WJ",
      "SSVM_WJ", "FSIQ")

cpt = c("N_of_omissions", "PC_omissions", "N_commissions", "PC_commissions",
      "hit_RT", "hit_RT_SE", "variability_of_SE", "D.prime", "beta",
      "N_perservations", "PC_perservations", "hit_RT_block_change",
      "hit_RT_SE_block_change", "hit_RT_ISI_change", "hit_RT_SE_ISI_change");    
wiscraw = c("RAW_DS_total", "RAW_DSB", "Raw_DSF", "Raw_SS_total", "Raw_SSB",
      "Raw_SSF")
wiscstd = c("Standard_DS_total", "Standard_DSB", "Standard_DSF",
            "Standard_SS_total", "Standard_SSB", "Standard_SSB.1") 
wj = c("PS", "SSDS_WJ", "SSVM_WJ")
iq = c("FSIQ")

df$mask.id = df$MRN

for (m in c('cpt', 'wiscraw', 'wiscstd', 'wj', 'iq')) {
    eval(parse(text=sprintf('data = df[, c("mask.id", %s)]', m)))
    keep_me = rowSums(is.na(data)) == 0
    data = data[keep_me, ]
    cnames = sapply(colnames(data)[2:ncol(data)], function(y) sprintf('v_%s', y))
    colnames(data)[2:ncol(data)] = cnames
    fname = sprintf('~/data/baseline_prediction/cog_%s_09142018.RData.gz', m)
    save(data, file=fname, compress=T)
    print(sprintf('%s has n=%d', m, nrow(data)))
}

# combine everything with a decent amount of subjects
load('~/data/baseline_prediction/cog_iq_09142018.RData.gz')
all_data = data
for (m in c('cpt', 'wiscraw', 'wj')) {
    fname = sprintf('~/data/baseline_prediction/cog_%s_09142018.RData.gz', m)
    load(fname)
    all_data = merge(all_data, data, by='mask.id')
}
print(sprintf('combined has n=%d', nrow(all_data)))
data = all_data
save(data, file='~/data/baseline_prediction/cog_all_09142018.RData.gz', compress=T)

clin = read.csv('~/data/baseline_prediction/long_clin_0913.csv')
colnames(clin)[1]='mask.id'
write.csv(clin, file='~/data/baseline_prediction/cog_gf_09142018.csv', row.names=F)
```


