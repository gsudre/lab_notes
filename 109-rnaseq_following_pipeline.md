# 2020-05-15 07:33:48

Let's continue our observations from 108 and follow a pipeline. Starting with
this one:

https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

```r
library(edgeR)
myregion = 'ACC'
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
counts = t(data[data$Region==myregion, grex_vars])
d0 <- DGEList(counts)
d0 <- calcNormFactors(d0)
cutoff <- 1  # subjective... play with this
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
group = data[data$Region==myregion, 'Diagnosis']
plotMDS(d, col = as.numeric(group))
```

![](images/2020-05-15-07-44-50.png)

```r
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
```

![](images/2020-05-15-07-46-05.png)

This is supposedly a good fit, base don the tutorial.

```r
fit <- lmFit(y, mm)
contr <- makeContrasts(groupCase - groupControl, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
```

This seems to finish the analysis, but as expected we don't have anything useful
after corretion with FDR.

We could check for pH, RIN and batch variables. What happens if I just chuck
them into the model?

```r
RIN = data[data$Region==myregion, 'RINe']
pH = data[data$Region==myregion, 'pH']
site = data[data$Region==myregion, 'bainbank']
batch = factor(data[data$Region==myregion, 'run_date'])
```

pH has many NAs, and I cannot add site and batch at the same time without
altering the variables, as there as cross-overs in the categories with zero
samples. It's just a matter of recoding, or doing the batch correction twice,
but for now we're just playing.

```r
mm <- model.matrix(~0 + group + batch + RIN)
y <- voom(d, mm, plot = T)
```

![](images/2020-05-15-08-05-27.png)

```r
fit <- lmFit(y, mm)
contr <- makeContrasts(groupCase - groupControl, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```

```r
fit <- lmFit(y, mm)
contr <- makeContrasts(groupCase - groupControl, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```

Still nothing. Maybe there really isn't anything if we use FDR. But we do have
plenty of other things to try, and we're just starting...

## Effects of combat

```r
library(sva)
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[data$Region=='ACC', ]
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
batch = as.numeric(data$run_date)
group = as.numeric(data$Diagnosis)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)
# now I'll further adjust it for brain bank
batch = as.numeric(data$bainbank)
adjusted_counts2 <- ComBat_seq(adjusted_counts, batch=batch, group=group)

d0 <- DGEList(adjusted_counts2)
d0 <- calcNormFactors(d0)
cutoff <- 1  # subjective... play with this
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
group = data[, 'Diagnosis']
plotMDS(d, col = as.numeric(group))
```

![](images/2020-05-15-09-08-03.png)

I forgot to remove 57...

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ] # removing ACC outlier
data = data[data$Region=='ACC', ]
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
batch = as.numeric(data$run_date)
group = as.numeric(data$Diagnosis)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)
# now I'll further adjust it for brain bank
batch = as.numeric(data$bainbank)
adjusted_counts2 <- ComBat_seq(adjusted_counts, batch=batch, group=group)

d0 <- DGEList(adjusted_counts2)
d0 <- calcNormFactors(d0)
cutoff <- 1  # subjective... play with this
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
group = data[, 'Diagnosis']
plotMDS(d, col = as.numeric(group))
```

![](images/2020-05-15-09-11-20.png)

This is better. Maybe it even correlates with population...

```r
mm <- model.matrix(~0 + group)
y <- voom(d, mm, plot = T)
```

![](images/2020-05-15-09-12-50.png)

```r
fit <- lmFit(y, mm)
contr <- makeContrasts(groupCase - groupControl, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```

![](images/2020-05-15-09-15-31.png)

This is looking much more promising... I could try being more stringent in the
cleaning, or adding a few more covariates, like RIN and PCs?

```r
RIN = data[, 'RINe']
mm <- model.matrix(~0 + group + RIN)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(groupCase - groupControl, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```

![](images/2020-05-15-09-18-11.png)

One hit... maybe adding PCs, or even keeping it to WNH only? Apparently we
cannot run more than 1 continuous variable in the model matrix? It's giving me
errors. Let's focus on the WNH only then.

```r
imWNH = data$C1 > 0 & data$C2 < -.075
data = data[which(imWNH),]
count_matrix = t(data[, grex_vars])
batch = as.numeric(data$run_date)
group = as.numeric(data$Diagnosis)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)
# now I'll further adjust it for brain bank
batch = as.numeric(data$bainbank)
adjusted_counts2 <- ComBat_seq(adjusted_counts, batch=batch, group=group)

d0 <- DGEList(adjusted_counts2)
d0 <- calcNormFactors(d0)
cutoff <- 1  # subjective... play with this
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]
group = data[, 'Diagnosis']
plotMDS(d, col = as.numeric(group))
```

![](images/2020-05-15-09-34-35.png)

There's some interesting separation there.

```r
mm <- model.matrix(~0 + Diagnosis + RINe, data=data)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(DiagnosisCase - DiagnosisControl,
                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```

![](images/2020-05-15-09-37-15.png)

Now we're talking!

```
> length(which(top.table$adj.P.Val < 0.05))
[1] 33
> rownames(top.table[top.table$adj.P.Val < 0.05, ])
 [1] "ENSG00000198692.10" "ENSG00000099725.14" "ENSG00000067048.17"
 [4] "ENSG00000206159.11" "ENSG00000067646.12" "ENSG00000129824.16"
 [7] "ENSG00000154620.6"  "ENSG00000099715.14" "ENSG00000228411.1" 
[10] "ENSG00000165246.14" "ENSG00000183878.15" "ENSG00000176728.9" 
[13] "ENSG00000131002.12" "ENSG00000226555.1"  "ENSG00000012817.15"
[16] "ENSG00000259917.1"  "ENSG00000196436.8"  "ENSG00000241859.7" 
[19] "ENSG00000114374.13" "ENSG00000251022.6"  "ENSG00000215580.11"
[22] "ENSG00000196584.3"  "ENSG00000260372.7"  "ENSG00000124782.20"
[25] "ENSG00000250483.1"  "ENSG00000260197.1"  "ENSG00000288049.1" 
[28] "ENSG00000229236.3"  "ENSG00000271741.1"  "ENSG00000217896.2" 
[31] "ENSG00000215414.4"  "ENSG00000223773.7"  "ENSG00000103995.14"
```

What if we select WNH after COMBAT?

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ] # removing ACC outlier
data = data[data$Region=='ACC', ]
grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
batch = as.numeric(data$run_date)
group = as.numeric(data$Diagnosis)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)
# now I'll further adjust it for brain bank
batch = as.numeric(data$bainbank)
adjusted_counts2 <- ComBat_seq(adjusted_counts, batch=batch, group=group)

imWNH = which(data$C1 > 0 & data$C2 < -.075)
d0 <- DGEList(adjusted_counts2)
d0 <- calcNormFactors(d0[, imWNH])
cutoff <- 1  # subjective... play with this
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

mm <- model.matrix(~0 + Diagnosis + RINe, data=data[imWNH, ])
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(DiagnosisCase - DiagnosisControl,
                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 20)
```

![](images/2020-05-15-10-44-50.png)

Results not as good...

Is it COMBAT or RINe doing the benefit?

```r
data = readRDS('~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
data = data[-c(which(rownames(data)=='57')), ] # removing ACC outlier
data = data[data$Region=='ACC', ]
imWNH = data$C1 > 0 & data$C2 < -.075
data = data[which(imWNH),]

grex_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
count_matrix = t(data[, grex_vars])
batch = as.numeric(data$run_date)
group = as.numeric(data$Diagnosis)
adjusted_counts <- ComBat_seq(count_matrix, batch=batch, group=group)
# now I'll further adjust it for brain bank
batch = as.numeric(data$bainbank)
adjusted_counts2 <- ComBat_seq(adjusted_counts, batch=batch, group=group)

d0 <- DGEList(adjusted_counts2)
d0 <- calcNormFactors(d0)
cutoff <- 1  # subjective... play with this
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,]

mm <- model.matrix(~0 + Diagnosis + RINe, data=data)
y <- voom(d, mm, plot = F)
fit <- lmFit(y, mm)
contr <- makeContrasts(DiagnosisCase - DiagnosisControl,
                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
length(which(top.table$adj.P.Val < 0.05))
```

I only get 20 without RIN:

![](images/2020-05-15-10-52-05.png)

Maybe we should use the intercept? Well, nothing wrong with use RINe as a
covariate...


# TODO
* try caudate
* how do these results look like in the entire population?
* look at gene functions
* play with other gene removal thresholds
* try other pipelines from note 108