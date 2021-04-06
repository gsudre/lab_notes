# 2021-03-31 20:22:34

Quick note on how I prepared the data for Gang, so he can play with his Bayesian
model.

I started compiled DTI data in all_DTI_tracts_11122020.csv, used only the
baseline scan for each subject, and the oldest person in each famID, to
hopefully get the best scan that way. This way we also remove any confounders of
same person or same family in the dataset. The data is already coded by SID, and
I removed other variables that I didn't think to be needed. That left 522
subjects in the end.

Then, we merge with our clinical variables:

```r
r$> source('~/research_code/lab_mgmt/merge_on_closest_date.R')                                 

r$> dti = read.csv('~/data/data_for_gang.csv')                                                 

r$> clin = read.csv('/Volumes/NCR/clinical_data_compiled//augmented_anon_clinical_02222021.csv'
    )     

    r$> clin = clin[clin$age_clin != 'child', ]                                                    

r$> clin$age_clin = as.numeric(clin$age_clin)                                                  

r$> m = mergeOnClosestAge(dti, clin, unique(dti$SID), x.id='SID', x.age = 'age_scan', y.id='SID
    ', y.age='age_clin') 

write.csv(m, row.names=F, file='~/data/data_for_gang_clin.csv')
```

I then only kept the DX_dsm variable, and renamed it as group.

# 2021-04-02 15:50:53

Gang complained the data was too complex. So, I'll reduce it to only the DTI-TK
atlas, and split it into the format he requested:

```r
df = read.csv('~/data/data_for_gang_clin.csv') 
df = df[, 1:41]
df2 = df[1, 1:8]
df2$FA = NA
df2$AD = NA
df2$RD = NA
df2$tract = NA
df2 = df2[, c(1, 12, 2:11)]
tracts = sapply(colnames(df)[9:ncol(df)],
                function(x) { res = strsplit(x, '_')[[1]];
                              paste0(res[2:length(res)], collapse='_') })
tracts = unique(tracts)
cnt = 1
for (s in unique(df$SID)) {
    subj_row = which(df$SID == s)
    for (t in tracts) {
        df2[cnt, 'tract'] = t
        df2[cnt, 'SID'] = s
        df2[cnt, 3:9] = df[subj_row, 2:8]
        for (m in c('FA', 'AD', 'RD')) {
            df2[cnt, m] = df[subj_row, sprintf('%s_%s', m, t)]
        }
        cnt = cnt + 1
    }
}
write.csv(df2, row.names=F, file='~/data/data_for_gang_simple.csv')
```