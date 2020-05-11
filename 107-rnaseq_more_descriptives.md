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

# 2020-05-11 07:55:56

```
> data3 = data2[data2$Manner.of.Death.x!='unknown' & data2$Manner.of.Death.x!='Suicide (probable)',]
> data3$mod = factor(data3$Manner.of.Death.x)
> t = table(data3$Diagnosis.x, data3$mod)
> print(t)
         
          Accident Homicide natural Suicide
  Case          10        1       7       6
  Control        7        9      10       3
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 8.0588, df = NA, p-value = 0.05147

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 8.0588, df = 3, p-value = 0.04481

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect
> t = table(data2$Diagnosis.x, data2$comorbid_groupmod)
data2$comorbid_group
> t = table(data2$Diagnosis.x, data2$comorbid_group)
> print(t)
         
          no yes
  Case    15  10
  Control 30   0
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 14.667, df = NA, p-value = 0.0004998

> chisq.test(t)

	Pearson's Chi-squared test with Yates' continuity correction

data:  t
X-squared = 12.101, df = 1, p-value = 0.0005039

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect
> t = table(data2$Diagnosis.x, data2$substance_group)
> print(t)
         
           0  1  2
  Case    13  7  5
  Control 30  0  0
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 18.419, df = NA, p-value = 0.0004998

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 18.419, df = 2, p-value = 0.0001001

Warning message:
In chisq.test(t) : Chi-squared approximation may be incorrect
```

So, subtance_group and comorbid_group will definitely have an impact in
Diagnosis, and possibly manner of death.

Looking at Caudate samples might have an impact because it's slightly larger.

## Caudate

```
> t = table(data2$Diagnosis.x, data2$bainbank.x)
> print(t)
         
          nimh_hbcc pitt umbn
  Case           11    5    9
  Control        14    5   14
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 0.35017, df = NA, p-value = 0.8286

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 0.35017, df = 2, p-value = 0.8394
> idx = data2$Diagnosis.x =='Case'
> t.test(data2[idx,'pcnt_optical_duplicates.x'], data2[!idx, 'pcnt_optical_duplicates.x'])

	Welch Two Sample t-test

data:  data2[idx, "pcnt_optical_duplicates.x"] and data2[!idx, "pcnt_optical_duplicates.x"]
t = 0.61078, df = 55.997, p-value = 0.5438
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4592760  0.8621851
sample estimates:
mean of x mean of y 
 3.086000  2.884545 

> wilcox.test(data2[idx,'pcnt_optical_duplicates.x'], data2[!idx, 'pcnt_optical_duplicates.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "pcnt_optical_duplicates.x"] and data2[!idx, "pcnt_optical_duplicates.x"]
W = 460.5, p-value = 0.4558
alternative hypothesis: true location shift is not equal to 0
> t.test(data2[idx,'clusters.x'], data2[!idx, 'clusters.x'])

	Welch Two Sample t-test

data:  data2[idx, "clusters.x"] and data2[!idx, "clusters.x"]
t = 0.99749, df = 53.931, p-value = 0.323
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -8858511 26400569
sample estimates:
mean of x mean of y 
119579674 110808645 

> wilcox.test(data2[idx,'clusters.x'], data2[!idx, 'clusters.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "clusters.x"] and data2[!idx, "clusters.x"]
W = 490, p-value = 0.2267
alternative hypothesis: true location shift is not equal to 0
> data2$rd2 = factor(data2$run_date.x)
> t = table(data2$Diagnosis.x, data2$rd2)
> print(t)
         
          1-May-19 16-Aug-19 25-Jun-19 28-Jun-19 3-May-19
  Case           3         3         1        11        7
  Control        4         9         1         9       10
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 2.8225, df = NA, p-value = 0.6437

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 2.8225, df = 4, p-value = 0.588
> t = table(data2$Diagnosis.x, data2$date_sent_to_NISC.x)
> print(t)
         
          sent_to_NISC_02212019 sent_to_NISC_05222019 sent_to_NISC_07232019
  Case                       12                    12                     1
  Control                    22                    10                     1
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 2.0587, df = NA, p-value = 0.4078

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 2.0587, df = 2, p-value = 0.3572
> wilcox.test(data2[idx,'Age.x'], data2[!idx, 'Age.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "Age.x"] and data2[!idx, "Age.x"]
W = 350.5, p-value = 0.3341
alternative hypothesis: true location shift is not equal to 0

> t.test(data2[idx,'Age.x'], data2[!idx, 'Age.x'])

	Welch Two Sample t-test

data:  data2[idx, "Age.x"] and data2[!idx, "Age.x"]
t = -0.96587, df = 51.821, p-value = 0.3386
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -6.327863  2.215813
sample estimates:
mean of x mean of y 
 20.89104  22.94707 
> t = table(data2$Diagnosis.x, data2$Sex.x)
> print(t)
         
           F  M
  Case     4 21
  Control  9 24
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 1.0394, df = NA, p-value = 0.3408

> chisq.test(t)

	Pearson's Chi-squared test with Yates' continuity correction

data:  t
X-squared = 0.49224, df = 1, p-value = 0.4829
> t.test(data2[idx,'RINe.x'], data2[!idx, 'RINe.x'])

	Welch Two Sample t-test

data:  data2[idx, "RINe.x"] and data2[!idx, "RINe.x"]
t = 0.96981, df = 55.993, p-value = 0.3363
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1529306  0.4399609
sample estimates:
mean of x mean of y 
 5.992000  5.848485 

> wilcox.test(data2[idx,'RINe.x'], data2[!idx, 'RINe.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "RINe.x"] and data2[!idx, "RINe.x"]
W = 453, p-value = 0.5291
alternative hypothesis: true location shift is not equal to 0
> t.test(data2[idx,'PMI.x'], data2[!idx, 'PMI.x'])

	Welch Two Sample t-test

data:  data2[idx, "PMI.x"] and data2[!idx, "PMI.x"]
t = 2.1196, df = 39.594, p-value = 0.04036
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
  0.3551991 15.0222555
sample estimates:
mean of x mean of y 
 28.11600  20.42727 

> wilcox.test(data2[idx,'PMI.x'], data2[!idx, 'PMI.x'], exact=F)

	Wilcoxon rank sum test with continuity correction

data:  data2[idx, "PMI.x"] and data2[!idx, "PMI.x"]
W = 549, p-value = 0.03265
alternative hypothesis: true location shift is not equal to 0
> data3 = data2[data2$Manner.of.Death.x!='unknown' & data2$Manner.of.Death.x!='Suicide (probable)',]
> data3$mod = factor(data3$Manner.of.Death.x)
> t = table(data3$Diagnosis.x, data3$mod)
> print(t)
         
          Accident Homicide natural Suicide
  Case           9        1       7       7
  Control        8        9      11       4
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 7.1694, df = NA, p-value = 0.07196

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 7.1694, df = 3, p-value = 0.06669
> t = table(data2$Diagnosis.x, data2$comorbid_group)
> print(t)
         
          no yes
  Case    14  11
  Control 33   0
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 17.918, df = NA, p-value = 0.0004998

> chisq.test(t)

	Pearson's Chi-squared test with Yates' continuity correction

data:  t
X-squared = 15.17, df = 1, p-value = 9.827e-05
> t = table(data2$Diagnosis.x, data2$substance_group)
> print(t)
         
           0  1  2
  Case    13  6  6
  Control 33  0  0
> chisq.test(t, simulate.p.value=T)

	Pearson's Chi-squared test with simulated p-value (based on 2000
	replicates)

data:  t
X-squared = 19.972, df = NA, p-value = 0.0004998

> chisq.test(t)

	Pearson's Chi-squared test

data:  t
X-squared = 19.972, df = 2, p-value = 4.604e-05
```

The Caudate results were also emphatic about substance_group and comorbid_group.
It also highlighted PMI, which wasn't there for ACC. Manner of death is not that
significant anymore, but worth considering later.

## Are there siginficnace DX differences in removed genes?

As we remove genes because X number of subjects didn't have the gene expressed,
is X differently proportioned towards ADHD?

```r
library(caret)
data = readRDS('~/data/rnaseq_derek/complete_data_04292020.rds')
data = data[data$Region=='ACC',]
grex_names = colnames(data)[grepl(colnames(data), pattern='^ENS')]
pp = preProcess(data[, grex_names],
                method=c('zv', 'nzv', 'range'), rangeBounds=c(0,1))
a = predict(pp, data[, grex_names])
n0 = colSums(a==0)
nobad = n0 == 1
some_bad_grex = names(n0)[!nobad]

# how many of the genes with at least 1 bad subject, this subject is the bad one
data$bad_grex = rowSums(a[, some_bad_grex]==0)
idx = data$Diagnosis =='Case'
wilcox.test(data[idx, 'bad_grex'], data[!idx, 'bad_grex'])
```

Not really significant. What if we make it more stringent what we consider to be
a bad gene?

```r
idx = data$Diagnosis =='Case'
pvals = c()
for (i in 2:30) {
	nobad = n0 < i
	some_bad_grex = names(n0)[!nobad]
	data$bad_grex = rowSums(a[, some_bad_grex]==0)
	pvals = c(pvals,
			  wilcox.test(data[idx, 'bad_grex'], data[!idx, 'bad_grex'])$p.value)
}
```

```
> pvals
 [1] 0.7604294 0.7169008 0.6803029 0.6706853 0.6444429 0.6706853 0.6236354
 [8] 0.6236354 0.6236354 0.5978993 0.5921902 0.5921902 0.5921649 0.5752084
[15] 0.5640274 0.5257175 0.4887694 0.4991846 0.4815368 0.4683246 0.5044323
[22] 0.5097089 0.4939660 0.4288569 0.4336811 0.4052167 0.3690341 0.3308329
[29] 0.3474464
```

Nothing really there. In other words, the ADHD group doesn't have a different
number of removed grex than NVs.

But what if we do this in a per gene basis? Say we
consider a gene to be bad after we remove X subjects. Is X equally distributed
between Case and Control?

```r
n0 = colSums(a==0)
prop_sig = c()
for (bad_subjs in 2:(max(n0)-1)) {
	print(bad_subjs)
	imbad = n0 > bad_subjs
	bad_grex = names(n0)[imbad]
	pvals = c()
	for (g in bad_grex) {
		tmp = rep('bad', length=nrow(data))
		tmp[a[, g] > 0] = 'good'
		# is the proportion of good and bad samples different between ADHD and NV?
		pvals = c(pvals, chisq.test(table(data$Diagnosis, tmp))$p.value)
	}
	# record the percentage of genes, out of all bad ones, with a different 
	# proportion between ADHD and NV p < .05
	prop_sig = c(prop_sig, sum(pvals<.05)/length(pvals))
	print(prop_sig)
}
```

```
 [1] 0.017691063 0.019127774 0.020294561 0.020695567 0.021669031 0.021645022
 [7] 0.021340162 0.022059920 0.022155877 0.022008936 0.022505234 0.022597832
[13] 0.022656855 0.022974359 0.022173913 0.021839080 0.021057786 0.021099245
[19] 0.020827548 0.020703934 0.020641473 0.019849418 0.020817844 0.020824990
[25] 0.020319931 0.021942110 0.022773279 0.023127753 0.024009604 0.022875817
[31] 0.022873481 0.023052464 0.022827041 0.022094140 0.021074816 0.022538553
[37] 0.020000000 0.021084337 0.018518519 0.017307692 0.017316017 0.009828010
[43] 0.008746356 0.006944444 0.004366812 0.006329114 0.000000000 0.000000000
[49] 0.000000000
```

And let's repeat the same analysis for the Caudate:

```
 [1] 0.02069385 0.02204214 0.02326377 0.02423806 0.02544304 0.02439024
 [7] 0.02433628 0.02527233 0.02527021 0.02534268 0.02600867 0.02587413
[13] 0.02595251 0.02692233 0.02668568 0.02695764 0.02708804 0.02631579
[19] 0.02634913 0.02777778 0.02761795 0.02887865 0.02926669 0.02917399
[25] 0.03070508 0.03062474 0.03138815 0.02939782 0.02962206 0.02734375
[31] 0.02766085 0.02803129 0.02642857 0.02852737 0.03022670 0.03266788
[37] 0.03174603 0.03459459 0.03444181 0.03434610 0.03106509 0.02662230
[43] 0.02651515 0.01766004 0.01842105 0.02208202 0.01145038 0.01492537
[49] 0.02068966 0.01960784 0.01886792
```

Nothing really there...

# TODO
* maybe subgroup in the ADHD_type variables?