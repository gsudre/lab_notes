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

# 2021-04-07 11:18:47

We had a nice chat with Gang yesterday. He'll send me his code to test some
easier cases, like the effects of removing outliers (he doesn't think it's
needed), and checking the effects of using a different atlas, or doing only kids
or only adults. I should also check how it handles predicting continuous variables.

He'll take care of harder cases, like longitudinal data and familiar
relationships. I'll start by sending him some longitudinal data. Again, I'll
keep it to only one per family. I'll choose the person with the most scans. If
it's a tie, the one with highest average age with hopes that yields better
scans.

```r
df = read.csv('/Volumes/NCR/brain_data_compiled/all_DTI_tracts_03222021.csv')
df = df[, 1:55]
# removing 1090 for now because I need to reprocess her scan
df = df[df$maskid != 1090, ]
df2 = df[1, c(1:3, 8:9, 18:19, 22)]
df2$FA = NA
df2$AD = NA
df2$RD = NA
df2$tract = NA
df2 = df2[, c(1, 2, 12, 3:11)]
tracts = sapply(colnames(df)[23:ncol(df)],
                function(x) { res = strsplit(x, '_')[[1]];
                              paste0(res[2:length(res)], collapse='_') })
tracts = unique(tracts)
# variables to be repeated
var_repeat = colnames(df2)[c(1, 4:9)]
cnt = 1
for (f in unique(df$famID)) {
    fam_rows = which(df$famID == f)
    ftable = table(df[fam_rows,'SID'])
    max_scans = max(ftable)
    # if one subject with the most scans
    if (sum(ftable == max_scans) == 1) {
        s = names(which.max(ftable))
    } else { # resolve tie by overall oldest scans
        max_subjs = names(ftable[ftable == max_scans])
        ages = c()
        for (s in max_subjs) {
            subj_rows = which(df$SID == as.numeric(s))
            ages = c(ages, mean(df[subj_rows, 'age_scan']))
        }
        s = max_subjs[which.max(ages)]
    }

    # now that we have the subject for the family, grab all their scans
    subj_rows = which(df$SID == s)
    for (subj_row in subj_rows) {
        for (t in tracts) {
            df2[cnt, 'tract'] = t
            df2[cnt, 'SID'] = s
            df2[cnt, var_repeat] = df[subj_row, var_repeat]
            for (m in c('FA', 'AD', 'RD')) {
                df2[cnt, m] = df[subj_row, sprintf('%s_%s', m, t)]
            }
            cnt = cnt + 1
        }
    }
}

# now I'll add the everADHD field as the single group metric
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('/Volumes/NCR/clinical_data_compiled/augmented_anon_clinical_02222021.csv')
clin = clin[clin$age_clin != 'child', ]
clin$age_clin = as.numeric(clin$age_clin)
m = mergeOnClosestAge(df2, clin[, c('SID', 'age_clin', 'everADHD_dsm')],
                      unique(df2$SID), x.id='SID',
                      x.age = 'age_scan', y.id='SID', y.age='age_clin')
# shuffling some variables around
m = m[, c(1, 2, 4:9, 14, 3, 10:12)]
colnames(m)[9] = 'group'
m$group = ifelse(m$group == 'yes', 'ADHD', 'NV')
write.csv(m, row.names=F, file='~/data/long_data_for_gang.csv')
```

And just so I'm not looking this up all the time, here are the upload
instructions:

![](images/2021-04-07-12-03-34.png)


# 2021-05-04 14:31:01

Let's send Gang some data with family IDs in it as well:

```r
df = read.csv('/Volumes/NCR/brain_data_compiled/all_DTI_tracts_03222021.csv')
df = df[, 1:55]
# removing 1090 for now because I need to reprocess her scan
df = df[df$maskid != 1090, ]
df2 = df[1, c(1:3, 6, 8:9, 18:19, 22)]
df2$FA = NA
df2$AD = NA
df2$RD = NA
df2$tract = NA
df2 = df2[, c(4, 2, 6, 1, 3, 5, 7:13)]
tracts = sapply(colnames(df)[23:ncol(df)],
                function(x) { res = strsplit(x, '_')[[1]];
                              paste0(res[2:length(res)], collapse='_') })
tracts = unique(tracts)
# variables to be repeated
var_repeat = colnames(df2)[c(1:9)]
cnt = 1
for (subj_row in 1:nrow(df)) {
    for (t in tracts) {
        df2[cnt, 'tract'] = t
        df2[cnt, 'SID'] = df[subj_row, 'SID']
        df2[cnt, var_repeat] = df[subj_row, var_repeat]
        for (m in c('FA', 'AD', 'RD')) {
            df2[cnt, m] = df[subj_row, sprintf('%s_%s', m, t)]
        }
        cnt = cnt + 1
    }
}

# now I'll add the everADHD field as the single group metric
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
clin = read.csv('/Volumes/NCR/clinical_data_compiled/augmented_anon_clinical_02222021.csv')
clin = clin[clin$age_clin != 'child', ]
clin$age_clin = as.numeric(clin$age_clin)
m = mergeOnClosestAge(df2, clin[, c('SID', 'age_clin', 'everADHD_dsm')],
                      unique(df2$SID), x.id='SID',
                      x.age = 'age_scan', y.id='SID', y.age='age_clin')
# shuffling some variables around
m = m[, c(1:13, 15)]
colnames(m)[14] = 'group'
m$group = ifelse(m$group == 'yes', 'ADHD', 'NV')
write.csv(m, row.names=F, file='~/data/long_data_for_gang_with_FAMID.csv')
```