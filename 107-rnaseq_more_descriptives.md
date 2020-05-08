# 2020-05-08 09:04:52

Let's do some more descriptives in the data. First, which of the different
variables should be used when predicting ADHD. In other words, what covaries
with ADHD that's not grex?

```r
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_data_04292020.rds')
data = data[-c(which(rownames(data)=='57')), ] # removing ACC outlier
data = data[data$Region==myregion, ]
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[, 1:33]
more = more[!duplicated(more$hbcc_brain_id), ]
data2 = merge(data, more, by='hbcc_brain_id', all.x=F, all.y=F)
```

```
> colnames(more)[1:40]
 [1] "hbcc_brain_id"           "dg_name"                
 [3] "original_brain_number"   "Region"                 
 [5] "bainbank"                "submitted_name"         
 [7] "project"                 "machine_type"           
 [9] "pcnt_optical_duplicates" "clusters"               
[11] "read_length"             "experiment_type"        
[13] "run_date"                "batch"                  
[15] "date_sent_to_NISC"       "Pair"                   
[17] "Age"                     "Sex"                    
[19] "Race.x"                  "Ethnicity.x"            
[21] "RINe"                    "Diagnosis"              
[23] "RINe.1"                  "PMI"                    
[25] "pH"                      "Manner.of.Death"        
[27] "Source"                  "PMINTERVAL"             
[29] "ADHD_type"               "comorbid"               
[31] "comorbid_group"          "substance"              
[33] "substance_group"         "grex1"          
```

Based on exploring the data, it makes sense to evaluate these:

```
"bainbank", "pcnt_optical_duplicates" "clusters" run_date, batch, date_Sent_to_NISC, Age, Sex, Race.x, Ethnicity.x, RINe, PMI, pH, Manner.of.Death, 
comorbid_group, substance_group
```

```
> t = table(data2$Diagnosis.x, data2$bainbank.x)
> print(t)
         
          nimh_hbcc pitt umbn
  Case           10    5   10
  Control        14    4   12
> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 0.50926, df = 2, p-value = 0.7752

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 0.50926, df = NA, p-value = 0.8161
> idx = data2$Diagnosis.x =='Case'
> t.test(data2[idx,'pcnt_optical_duplicates.x'], data2[!idx, 'pcnt_optical_duplicates.x'])

	Welch Two Sample t-test

data:  data2[idx, "pcnt_optical_duplicates.x"] and data2[!idx, "pcnt_optical_duplicates.x"]
t = 1.5346, df = 52.512, p-value = 0.1309
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.2032867  1.5263534
sample estimates:
mean of x mean of y 
 3.727200  3.065667 

> wilcox.test(data2[idx,'pcnt_optical_duplicates.x'], data2[!idx, 'pcnt_optical_duplicates.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "pcnt_optical_duplicates.x"] and data2[!idx, "pcnt_optical_duplicates.x"]
W = 444.5, p-value = 0.2434
alternative hypothesis: true location shift is not equal to 0
> t.test(data2[idx,'clusters.x'], data2[!idx, 'clusters.x'])

	Welch Two Sample t-test

data:  data2[idx, "clusters.x"] and data2[!idx, "clusters.x"]
t = 0.64446, df = 45.014, p-value = 0.5226
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -11908095  23114402
sample estimates:
mean of x mean of y 
112646312 107043159 

> wilcox.test(data2[idx,'clusters.x'], data2[!idx, 'clusters.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "clusters.x"] and data2[!idx, "clusters.x"]
W = 355, p-value = 0.7417
alternative hypothesis: true location shift is not equal to 0

> data2$rd2 = factor(data2$run_date.x)
> t = table(data2$Diagnosis.x, data2$rd2)
> print(t)
         
          1-May-19 16-Aug-19 25-Jun-19 3-May-19
  Case           5         4        11        5
  Control        5         9         7        9
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 3.5294, df = NA, p-value = 0.3418

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 3.5294, df = 3, p-value = 0.317

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect
> t = table(data2$Diagnosis.x, data2$date_sent_to_NISC.x)
> print(t)
         
          sent_to_NISC_02212019 sent_to_NISC_05222019 sent_to_NISC_07232019
  Case                       13                    11                     1
  Control                    21                     7                     2
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 2.6721, df = NA, p-value = 0.3498

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 2.6721, df = 2, p-value = 0.2629

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect
> wilcox.test(data2[idx,'Age.x'], data2[!idx, 'Age.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "Age.x"] and data2[!idx, "Age.x"]
W = 332, p-value = 0.4724
alternative hypothesis: true location shift is not equal to 0

> t.test(data2[idx,'Age.x'], data2[!idx, 'Age.x'])

	Welch Two Sample t-test

data:  data2[idx, "Age.x"] and data2[!idx, "Age.x"]
t = -0.60125, df = 49.108, p-value = 0.5504
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -5.736798  3.094410
sample estimates:
mean of x mean of y 
 21.35026  22.67145 

> t = table(data2$Diagnosis.x, data2$Sex.x)
> print(t)
         
           F  M
  Case     2 23
  Control  8 22
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 3.1941, df = NA, p-value = 0.08346

> chisq.test(t)

	Pearson's Chi-squared test with Yates' continuity correction

data:  t
X-squared = 2.0625, df = 1, p-value = 0.151

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect

> t.test(data2[idx,'RINe.x'], data2[!idx, 'RINe.x'])

	Welch Two Sample t-test

data:  data2[idx, "RINe.x"] and data2[!idx, "RINe.x"]
t = -1.414, df = 48.862, p-value = 0.1637
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6941148  0.1207815
sample estimates:
mean of x mean of y 
 5.180000  5.466667 

> wilcox.test(data2[idx,'RINe.x'], data2[!idx, 'RINe.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "RINe.x"] and data2[!idx, "RINe.x"]
W = 281, p-value = 0.1133
alternative hypothesis: true location shift is not equal to 0

> t.test(data2[idx,'PMI.x'], data2[!idx, 'PMI.x'])

	Welch Two Sample t-test

data:  data2[idx, "PMI.x"] and data2[!idx, "PMI.x"]
t = 1.5131, df = 40.967, p-value = 0.1379
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.853186 12.925186
sample estimates:
mean of x mean of y 
   26.596    21.060 

> wilcox.test(data2[idx,'PMI.x'], data2[!idx, 'PMI.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "PMI.x"] and data2[!idx, "PMI.x"]
W = 461, p-value = 0.1482
alternative hypothesis: true location shift is not equal to 0

```

Philip said to use (refactored) run_date instead of batch.

There's some stuff missing in race and ethnicity, so maybe that will be better
handled by using the PCs anyways.


# TODO
* do these results change for the Caudate?
* maybe subgroup in the ADHD_type variables?