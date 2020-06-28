# 2020-06-28 10:49:57

Let's see how using comorbidities and substance abuse change our results. 

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
rownames(data) = data$submitted_name  # just to ensure compatibility later
rm_me = rownames(data) %in% c('68080')
data = data[!rm_me, ]
data = data[, !grepl(colnames(data), pattern='^ENS')]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
data = data[!duplicated(data$hbcc_brain_id),]
data = merge(data, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                            'substance', 'substance_group')],
             by='hbcc_brain_id', all.x=T, all.y=F)
# using the kmeans 3 clusters using 10 PCs result
imWNH = which(data$C1 > 0 & data$C2 < -.065)
dataWNH = data[imWNH, ]
```

At this point we have 60 subjects (most of them with both brain regions - 113
data points), and only 36 WNH.

```
> table(data$Diagnosis, data$comorbid_group)
         
          no yes
  Case    15  12
  Control 33   0
> table(dataWNH$Diagnosis, dataWNH$comorbid_group)
         
          no yes
  Case    13   8
  Control 15   0
```

The numbers above are in units of subjects. Because the comorbidity rate is much
higher in cases, it'd make sense to remove the people with comorbidity from the
analysis, instead of adding it as a covariate. We'd be losing a good chunk of
data, but I'm not sure how the models would run otherwise.

```
> table(data$Diagnosis, data$substance_group)
         
           0  1  2
  Case    14  7  6
  Control 33  0  0
> table(dataWNH$Diagnosis, dataWNH$substance_group)
         
           0  1  2
  Case    10  6  5
  Control 15  0  0
```

Same thing for substance abuse. I might need to remove them entirely. Let's see
what we get if I keep the data squeaky clean:

```r
datac = data[data$substance_group==0 & data$comorbid_group=='no',]
imWNH = which(datac$C1 > 0 & datac$C2 < -.065)
datacWNH = datac[imWNH, ]
```

```
> summary(datac$Diagnosis)
   Case Control 
      8      33 
> summary(datacWNH$Diagnosis)
   Case Control 
      6      15 
```

It's worth running it to see what happens. Or I can run it separately for each
cleaning. 