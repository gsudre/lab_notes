# 2020-10-02 11:57:43

I'm going to reproduce Alex's methylation analysis here, but using the updated
packages and also using a more principled way for selecting the PCs. His
original script is in data/Methylation_pipeline-Sep22.rmd

```r
library(doParallel)
library(data.table)
library(ChAMP)
library(ChAMPdata)
library(ewastools)
library(MethylToSNP)
library(qqman)
library(EnhancedVolcano)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(ggfortify)
library(ggbiplot)
library(stringr)
library(scales)
library(writexl)
library(limma)

registerDoParallel(cores = 3)

saveExcel<- function(v) {
    name <- deparse(substitute(v))
    writexl::write_xlsx(cbind(index=rownames(v), as.data.frame(v)),
                        paste0(name, ".xlsx"))
}

addalpha <- function(colors, alpha=1.0) {
  r <- col2rgb(colors, alpha=T)
  # Apply alpha
  r[4,] <- alpha*255
  r <- r/255.0
  return(rgb(r[1,], r[2,], r[3,], r[4,]))
}

rootDir <- "~/data/methylation_post_mortem"
rawDataDir <- paste(rootDir, '/Shaw_2019', sep="")
samplesFile <- paste(rootDir, '/samples.csv', sep="")
```

I'm going to stop here for now, but the next step is to transform Alex's samples
file to .csv to make sure it works and loads the data properyl. Stopping because
I want to check that the PCA analysis works on my RNAseq data as well before
running this to send a list of genes and probes to David.

# 2020-10-07 13:27:42

Let's pick it up from where we left of.

```r
# disable scientific notation
options(scipen = 999)
samples <- read.csv(samplesFile)
samples$Sample_Group <- as.factor(as.character(samples$Diagnosis))
samples$Basename <- paste0(rawDataDir, '/',  samples$Sentrix_ID, '/',
                           samples$Sentrix_ID,'_',samples$Sentrix_Position)
raw = read_idats(samples$Basename, quiet = FALSE)
raw$sample_names <- samples$Sample_Name
```

Loading the samples took a few minutes, but it looks like it worked. Moving on
to QC:

```r
# checks control metrics according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf
ctrls = control_metrics(raw)
stripchart(ctrls$`Bisulfite Conversion I Green`, method="jitter", pch=4,
           xlab='Bisulfite Conversion I Green',xlim=c(0,10))
abline(v=1,col=2,lty=3)

stripchart(ctrls$`Bisulfite Conversion I Red`, method="jitter", pch=4,
           xlab='Bisulfite Conversion I Red',xlim=c(0,10))
abline(v=1,col=2,lty=3)

stripchart(ctrls$`Bisulfite Conversion II`, method="jitter", pch=4,
           xlab='Bisulfite Conversion II',xlim=c(0,10))
abline(v=1,col=2,lty=3)

failed.samples <- samples[sample_failure(ctrls),
                          c("Sample_Group", "Sentrix_ID")]
```

![](images/2020-10-07-20-52-26.png)
![](images/2020-10-07-20-52-44.png)
![](images/2020-10-07-20-53-03.png)

QC looks good, no failed samples. Checking sex:

```r
pheno = data.table(samples)
pheno$sex <- tolower(pheno$Sex)
pheno[, c("X","Y") := check_sex(raw)]
pheno[, predicted_sex := predict_sex(X, Y, which(sex=="m"), which(sex=="f"))]

tmp = pheno[predicted_sex==sex]
plot(Y ~ X, data=tmp, pch=ifelse(tmp$sex=="f",1,4), asp=1,
     xlab="Normalized X chromosome intensities",
     ylab="Normalized Y chromosome intensities")
tmp = pheno[predicted_sex!=sex]
```
![](images/2020-10-07-20-59-23.png)

There were no mismatches between predicted and declared sex.

# 2020-10-08 06:08:39

Let's continue with preprocessing:

```r
#  Preprocess with normalize450k / ewastools that includes bg and dye corrections
meth = raw %>% detectionP %>% mask(0.01) %>% correct_dye_bias() %>% dont_normalize
# split filtered results between markers with rs in their name and without
colnames(meth) <- samples$Sample_Name
meth.rs <- meth[grep("rs", rownames(meth)),]
meth.nors <-meth[grep("rs", rownames(meth), invert=TRUE),]
# Filter probes: sex chr, SNP (general, EUR, AFR), CpG, Multi-hit
# apply typical champ filters, including common general SNPs
filtered <- champ.filter(beta=meth.nors, pd=samples, filterDetP = FALSE,
                         filterXY = TRUE, autoimpute=FALSE, filterNoCG = TRUE,
                         filterMultiHit = TRUE, filterBeads = FALSE,
                         fixOutlier=FALSE, arraytype="EPIC")
meth.filtered <- na.omit(filtered$beta)
dim(meth.filtered)
```

We are now down to 647K markers, from the initial 865K. I'll skip checking the
QC control samples because Alex has already done it, and he didn't include them
in his samples spreadsheet. I'd have to reconstruct the sheet from the sample
sheet of each box to redo it. Not worth it if he has already confirmed that the
QC controls match...

So, let's get an idea about the data we're working with:

```r
table(samples$Sex)
table(samples$Diagnosis)
table(samples$Brainbank)
table(samples$Manner.of.Death)
table(samples$Kit)
```

```
F  M 
21 94 

   Case Control 
     51      64 

nimh_hbcc      pitt      umbn 
       50        19        46 

          Accident           Homicide            natural            Suicide Suicide (probable)            unknown
                32                 20                 38                 21                  2                  2

Qiagen   Zymo
    59     56
```

I'll do the same modifications I did for the RNAseq data:

```r
samples$Sentrix_ID <- as.factor(samples$Sentrix_ID)
samples$Sample_Group <- as.factor(samples$Sample_Group)
samples$Row <- as.factor(substr(samples$Sentrix_Position, 1, 3))
samples$Scanner <- as.factor(samples$Scanner)
samples$Region <- as.factor(samples$Region)
samples$Sex <- as.factor(samples$Sex)
samples$Kit <- as.factor(samples$Kit)
samples$Sample_Plate <- as.factor(samples$Sample_Plate)

more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[!duplicated(more$hbcc_brain_id),]
samples = merge(samples, more[, c('hbcc_brain_id', 'comorbid', 'comorbid_group',
                                  'substance', 'substance_group', 'C1', 'C2',
                                  'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9',
                                  'C10', 'POP_CODE')],
             by='hbcc_brain_id', all.x=T, all.y=F)
samples$POP_CODE = as.character(samples$POP_CODE)
samples[which(samples$POP_CODE=='WNH'), 'POP_CODE'] = 'W'
samples[which(samples$POP_CODE=='WH'), 'POP_CODE'] = 'W'
samples$POP_CODE = factor(samples$POP_CODE)
samples[which(samples$Manner.of.Death=='Suicide (probable)'),
        'Manner.of.Death'] = 'Suicide'
samples[which(samples$Manner.of.Death=='unknown'),
        'Manner.of.Death'] = 'natural'
samples$MoD = factor(samples$Manner.of.Death)
```

Before we run the PCA analysis, let's split the data between ACC and Caudate.

```r
is.caudate <- samples$Region == 'Caudate'
is.outlier <- samples$Sample_Name == "1908_ACC"

samples.caudate <- samples[is.caudate, ]
samples.acc <- samples[!is.caudate & !is.outlier, ]

meth.caudate <- meth.filtered[, is.caudate]
meth.acc <- meth.filtered[, !is.caudate & !is.outlier]

# Stratified PCA
set.seed(42)
meth.caudate.pca <- prcomp(t(meth.caudate), scale.=TRUE)
meth.acc.pca <- prcomp(t(meth.acc), scale.=TRUE)
```

Like in the RNAseq analysis, let's figure out how many PCs to use for ACC and
Caudate.

```r
library(nFactors)
eigs <- meth.acc.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
pcs.acc = data.frame(meth.acc.pca$x[, keep_me])

std_dev <- meth.acc.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, xlab = "Principal Component",
             ylab = "Proportion of Variance Explained",
             type = "b", main='ACC')


eigs <- meth.caudate.pca$sdev^2
nS = nScree(x=eigs)
keep_me = 1:nS$Components$nkaiser
pcs.caudate = data.frame(meth.caudate.pca$x[, keep_me])

std_dev <- meth.caudate.pca$sdev
pr_var <- std_dev^2
prop_varex <- pr_var/sum(pr_var)
plot(prop_varex, xlab = "Principal Component",
             ylab = "Proportion of Variance Explained",
             type = "b", main='Caudate')
```

We got 9 for ACC and 13 for Caudate.

![](images/2020-10-08-07-07-56.png)
![](images/2020-10-08-07-09-19.png)

Now let's make the same plots as the RNAseq analysis, but most importantly we
need to see which PCs are correlated with the data at a certain threshold.

```r
num_vars = c('pH', 'Age', 'RINe', 'PMI', 'RIN',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(pcs.acc)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(samples.acc[, x], pcs.acc[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

library(corrplot)
corrplot(t(num_corrs), method='color', tl.cex=1, cl.cex=1)
```

![](images/2020-10-08-07-24-01.png)

```r
categ_vars = c('MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex',
               "Sentrix_ID", "Row", "Scanner", "Kit", "Sample_Group",
               "Brainbank")
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(pcs.acc[, y], samples.acc[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}
corrplot(t(categ_corrs), method='color', tl.cex=1, cl.cex=1, is.corr=F)
```

![](images/2020-10-08-07-23-39.png)

```
r$> which(num_pvals < .01, arr.ind = T)                                                                                                                                        
    row col
pH    1   2
pH    1   4
C2    7   5
C3    8   5
C4    9   5
C7   12   5
RIN   5   6

r$> which(categ_pvals < .01, arr.ind = T)                                                                                                                                      
           row col
Sentrix_ID   6   2
Brainbank   11   2
Kit          9   5
Brainbank   11   5
Brainbank   11   9
```

So, for ACC we will remove PCs 2, 4, 5, 6, and 9. And the minimum p-value for
Sample_Group was .19, so we're good there. Let's check the Caudate:

```r
num_vars = c('pH', 'Age', 'RINe', 'PMI', 'RIN',
             'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'C10')
pc_vars = colnames(pcs.caudate)
num_corrs = matrix(nrow=length(num_vars), ncol=length(pc_vars),
                   dimnames=list(num_vars, pc_vars))
num_pvals = num_corrs
for (x in num_vars) {
    for (y in pc_vars) {
        res = cor.test(samples.caudate[, x], pcs.caudate[, y])
        num_corrs[x, y] = res$estimate
        num_pvals[x, y] = res$p.value
    }
}

library(corrplot)
corrplot(t(num_corrs), method='color', tl.cex=1, cl.cex=1)
```

![](images/2020-10-08-07-28-24.png)

```r
categ_vars = c('MoD', 'substance_group',
               'comorbid_group', 'POP_CODE', 'Sex',
               "Sentrix_ID", "Row", "Scanner", "Kit", "Sample_Group",
               "Brainbank")
categ_corrs = matrix(nrow=length(categ_vars), ncol=length(pc_vars),
                   dimnames=list(categ_vars, pc_vars))
categ_pvals = categ_corrs
for (x in categ_vars) {
    for (y in pc_vars) {
        res = kruskal.test(pcs.caudate[, y], samples.caudate[, x])
        categ_corrs[x, y] = res$statistic
        categ_pvals[x, y] = res$p.value
    }
}
corrplot(t(categ_corrs), method='color', tl.cex=1, cl.cex=1, is.corr=F)
```

![](images/2020-10-08-07-29-19.png)

```
r$> which(num_pvals < .01, arr.ind = T)                                                                                                                                        
    row col
pH    1   1
C9   14   2
C3    8   5
PMI   4   6
C2    7  12
C4    9  12

r$> which(categ_pvals < .01, arr.ind = T)                                                                                                                                      
           row col
Sentrix_ID   6   1
Brainbank   11   1
Brainbank   11   5
Brainbank   11   6
Scanner      8  12
```

And for Caudate we will remove PCs 1, 2, 5, 6, and 12. And the minimum p-value for
Sample_Group was .1, so we're also good there.