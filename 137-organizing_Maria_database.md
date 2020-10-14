# 2020-10-13 09:40:23

I want to create a big file with all subejcts in the USA Trios study. The
subjects should not be part of MTA and Paisas, so let's sart with just that one
as we only have phenotypes for it at the moment.

There are demographic data in development_history.txt, which I'll need to verify
against what Ariel sent. 

For now, let's make sure the study code matches the person and family code, then
we can check what kind of data each person has.

```
r$> df = read.table('~/data/all_families/wender.txt', sep = "\t", header=1)                                                          

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))                                                                      
[1] 0

r$> df = read.table('~/data/all_families/vanderbilt.txt', sep = "\t", header=1)

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))
[1] 0

r$> df = read.table('~/data/all_families/swan.txt', sep = "\t", header=1)

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))
[1] 0

r$> df = read.table('~/data/all_families/scid.txt', sep = "\t", header=1)

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))
[1] 0

r$> df = read.table('~/data/all_families/restrospective_dica.txt', sep = "\t", header=1)                                             

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))                                                                      
[1] 0

r$> df = read.table('~/data/all_families/dica.txt', sep = "\t", header=1)                                                            

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))                                                                      
[1] 0

r$> df = read.table('~/data/all_families/caars.txt', sep = "\t", header=1)                                                           

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))                                                                      
[1] 0

r$> df = read.table('~/data/all_families/development_history.txt', sep = "\t", header=2)                                             

r$> sum(df$Study != sprintf('%04d-%02d', df$Family, df$Member))                                                                      
[1] 0
```

Note that I had to remove one header cell from the developmen file in Excel to
make it work. I just did it in my local copy though.

So, it looks like the code in the database was automatically generated and it
all works. 

I'm going to start with the devlopmental history spreadsheet because it has all
the demographic info, and then add columns to that. Note that not everyone in
that spreadsheet has dev_history either, but at least the records are there.

```r
df = read.table('~/data/all_families/development_history.txt', sep = "\t",
                header=2)
data = df[, 1:11]
data = data[!duplicated(data$Study), ]
data$has_devHistory = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_devHistory'] = T

df = read.table('~/data/all_families/caars.txt', sep = "\t", header=1)
data$has_caars = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_caars'] = T
# also make sure everyone in CAARs is already in the main data table
sum(! imgood %in% data$Study) 

df = read.table('~/data/all_families/dica.txt', sep = "\t", header=1)
data$has_dica = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_dica'] = T
sum(! imgood %in% data$Study) 

df = read.table('~/data/all_families/restrospective_dica.txt', sep = "\t", header=1)
data$has_retrospectiveDica = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_retrospectiveDica'] = T
sum(! imgood %in% data$Study) 

df = read.table('~/data/all_families/scid.txt', sep = "\t", header=1)
data$has_scid = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_scid'] = T
sum(! imgood %in% data$Study) 

df = read.table('~/data/all_families/swan.txt', sep = "\t", header=1)
data$has_swan = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_swan'] = T
sum(! imgood %in% data$Study) 

df = read.table('~/data/all_families/vanderbilt.txt', sep = "\t", header=1)
data$has_vanderbilt = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_vanderbilt'] = T
sum(! imgood %in% data$Study) 

df = read.table('~/data/all_families/wender.txt', sep = "\t", header=1)
data$has_wender = F
imgood = df[df$Q1 != '', 'Study']
data[data$Study %in% imgood, 'has_wender'] = T
sum(! imgood %in% data$Study)

write.csv(data, file='~/data/USA_trios_demo.csv', row.names=F)
```

Now we should check the genotyped data. The tricky thing here is that even
though one of the projects islabeled Colombia and the other USA, all samples in
the Colombia sheet are marked as USA, and it's a mix in the USA set. So, let me
try to do this by BSN...

```
r$> df = read.xls('/Volumes/NCR/from_Maria//Broad_family_samples/PGC-ADHD_MuenkeAcosta_USAtrios_Summary Report_GSA-MD_FINAL_sendout_2
    -20-20.xlsx', 'phenofile')                                                                                                       

r$> df$BSN = gsub(x=df$Collaborator_Sample_ID, pattern='Muenke1-', replacement = '')                                                 
r$> table(df$BSN %in% data$BSN)                                                                                                      

FALSE  TRUE 
   14   826 
```

So, the first issue here is that out of the 840 samples, the BSN of 14 of them
is not in the phenotype file. That might not necessarily be bad, as it could be
that those people didn't get phenotyped. But we'd need to at least figure out
their identity and some demographic data. Let's keep on checking though.

```
r$> data2 = merge(data, df, by='BSN', all.x=F, all.y=F)
r$> sum(data2$SEX != data2$Gender)                                                                                                   
[1] 10

r$> data2[data2$ADHD=='Y','ADHD'] = 'ADHD'                                                                                           

r$> data2[data2$ADHD=='N','ADHD'] = 'Control'  
r$> sum(data2$ADHD != data2$Primary_Disease)                                                                                         
[1] 47
```

Then, if we merge the data based on the BSNs we found, we get 10 people with
wrong Sex, and 47 with wrong DX. Some of those could actually be wrong coding,
but it's worrisome because I'm not completely confident that the sample number
sent to Broad is the BSN.

I'm going to go ahead and merge in the genotype information, but add a few flags
to mismatches. Also, I looked at the samples in the Colombia Broad file, and
their BSNs are indeed quite low and would in that case match the actual
Colombian samples in the Max spreadsheet in the Sharepoint folder.

```r
df = read.xls('/Volumes/NCR/from_Maria//Broad_family_samples/PGC-ADHD_MuenkeAcosta_USAtrios_Summary Report_GSA-MD_FINAL_sendout_2
    -20-20.xlsx', 'phenofile')
df$BSN = gsub(x=df$Collaborator_Sample_ID, pattern='Muenke1-', replacement = '')
data2 = merge(data, df[, c('BSN', 'Gender', 'Primary_Disease')], all.x=T)
r$> df[which(! df$BSN %in% data2$BSN), 'BSN']                                                                                        
 [1] "9645" "9971" "9573" "9646" "9931" "9028" "9930" "9030" "7754" "9647" "9968" "8392" "8335" "9969"
```

So, I still need to add those 14 BSNs by hand to my file. And then maybe some
other coming from Ben's file that might have BSNs but no phenotype or who have
been previously genotyped. Let's do some formatting in Excel then, getting some
information from USA_SAMPLE_TRIOS.xlsx.

I'm gonna go ahead and redo the merge with the genotype data in Excel. But then
I found a few BSNs assigned to 2 different people. Gonna try using the USA Trios
spreadsheets to disambiguate. In looking at those spreadsheets, and also Max's
spreadsheet, most of the errors are in the phenotype (development_history) file.
But I'm correcting them in the Excel USA_trios_demo.xlsx in sharepoint. All the
DX mismatches were errors in the Broad phenotype file, as both Max's and the USa
Trios spreadsheet agreed with what was in the phenotype file. Now we just need
to add the information on the samples we have an we're good to go.

