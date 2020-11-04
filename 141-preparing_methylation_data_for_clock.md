# 2020-11-02 17:09:26

Let's format the methylation data so Gauri can run it in the epigenetic clock
from Shireby2020. I'll use my note on reproducing Alex's analysis as a model.

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

# disable scientific notation
options(scipen = 999)
samples <- read.csv(samplesFile)
samples$Sample_Group <- as.factor(as.character(samples$Diagnosis))
samples$Basename <- paste0(rawDataDir, '/',  samples$Sentrix_ID, '/',
                           samples$Sentrix_ID,'_',samples$Sentrix_Position)
raw = read_idats(samples$Basename, quiet = FALSE)
raw$sample_names <- samples$Sample_Name

# checks control metrics according to https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf
ctrls = control_metrics(raw)
failed.samples <- samples[sample_failure(ctrls),
                          c("Sample_Group", "Sentrix_ID")]

pheno = data.table(samples)
pheno$sex <- tolower(pheno$Sex)
pheno[, c("X","Y") := check_sex(raw)]
pheno[, predicted_sex := predict_sex(X, Y, which(sex=="m"), which(sex=="f"))]

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

Before we start using the PCs for cleaning, let me go ahead and export the
uncleaned data...

In this process I found out that there was an error in the methylation analysis!
When I merged the samples with the pop and PC data it resorted the matrix,
making it lose the connection with the methyl matrix! Will start a new note to
fix that!

# 2020-11-03 15:10:19

In note 142 I already exported some data for the clock, but there were some
probes missing, so I want to make sure they weren't removed in the cleaning:

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

rootDir <- "~/data/methylation_post_mortem"
rawDataDir <- paste(rootDir, '/Shaw_2019', sep="")
samplesFile <- paste(rootDir, '/samples.csv', sep="")

# disable scientific notation
options(scipen = 999)
samples <- read.csv(samplesFile)
samples$Sample_Group <- as.factor(as.character(samples$Diagnosis))
samples$Basename <- paste0(rawDataDir, '/',  samples$Sentrix_ID, '/',
                           samples$Sentrix_ID,'_',samples$Sentrix_Position)
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

raw = read_idats(samples$Basename, quiet = FALSE)
raw$sample_names <- samples$Sample_Name

#  Preprocess with normalize450k / ewastools that includes bg and dye corrections
meth = raw %>% detectionP %>% mask(0.01) %>% correct_dye_bias() %>% dont_normalize
# split filtered results between markers with rs in their name and without
colnames(meth) <- samples$Sample_Name
meth.rs <- meth[grep("rs", rownames(meth)),]
meth.nors <-meth[grep("rs", rownames(meth), invert=TRUE),]
```

Then we check if we have all probes:

```r
clock_markers = read.table('~/data/methylation_post_mortem/CorticalClockCoefs.tx
    t', head=1)

clock_data = meth.nors[rownames(meth.nors) %in% clock_markers$probe, ]

is.caudate <- samples$Region == 'Caudate'
is.outlier <- samples$Sample_Name == "1908_ACC"

samples.caudate <- samples[is.caudate, ]
samples.acc <- samples[!is.caudate & !is.outlier, ]

clock.caudate <- clock_data[, is.caudate]
clock.acc <- clock_data[, !is.caudate & !is.outlier]

new_names = gsub(x=colnames(clock.acc), pattern='_ACC', replacement='')
colnames(clock.acc) = new_names
new_names = gsub(x=colnames(clock.caudate), pattern='_Caudate', replacement='')
colnames(clock.caudate) = new_names
saveRDS(clock.acc, file='~/data/methylation_post_mortem/acc_methylClockOnly_raw_11032020.rds')
saveRDS(clock.caudate, file='~/data/methylation_post_mortem/caudate_methylClockOnly_raw_11032020.rds')

ages = samples$Age
names(ages) = samples$Sample_Name
new_names = gsub(x=names(ages), pattern='_ACC', replacement='')
new_names = gsub(x=new_names, pattern='_Caudate', replacement='')
names(ages) = new_names
ages = ages[!duplicated(new_names)]
df = data.frame(IID=names(ages), age=ages)
write.csv(df, file='~/data/methylation_post_mortem/IID_age.csv', row.names=F)
```

Gauri ran the new sets, so let's see if there is any difference between controls
and cases. First, is there a difference just on age?

```
r$> preds = read.csv('~/data/clock//ACC_Caudate_allprobes_GS_11032020/ACC_allprobes_
    CorticalClock_GS_11032020/CorticalPred_acc.csv')      
r$> idx = samples.acc$Sample_Group == 'Case'                                        

r$> t.test(samples.acc$Age[idx], samples.acc$Age[!idx])                             

        Welch Two Sample t-test

data:  samples.acc$Age[idx] and samples.acc$Age[!idx]
t = -0.68821, df = 50.461, p-value = 0.4945
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -5.866391  2.871714
sample estimates:
mean of x mean of y 
 20.83026  22.32759 


r$> idx = samples.caudate$Sample_Group == 'Case'                                    

r$> t.test(samples.caudate$Age[idx], samples.caudate$Age[!idx])                     

        Welch Two Sample t-test

data:  samples.caudate$Age[idx] and samples.caudate$Age[!idx]
t = -0.9081, df = 53.232, p-value = 0.3679
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -6.213085  2.340186
sample estimates:
mean of x mean of y 
 21.01062  22.94707 
```

OK, so that's good. Now, let's import the predicted ages and see what we get:

```
r$> t.test(preds$brainpred[idx], preds$brainpred[!idx])                             

        Welch Two Sample t-test

data:  preds$brainpred[idx] and preds$brainpred[!idx]
t = -0.29585, df = 52.107, p-value = 0.7685
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -5.736057  4.261950
sample estimates:
mean of x mean of y 
 21.97844  22.71550 


r$> d = preds$Age - preds$brainpred                                                 

r$> t.test(d[idx], d[!idx])                                                         

        Welch Two Sample t-test

data:  d[idx] and d[!idx]
t = -0.92695, df = 53, p-value = 0.3582
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -2.4053947  0.8848248
sample estimates:
 mean of x  mean of y 
-1.1481867 -0.3879018 


r$> d = abs(preds$Age - preds$brainpred)                                            

r$> t.test(d[idx], d[!idx])                                                         

        Welch Two Sample t-test

data:  d[idx] and d[!idx]
t = -0.29228, df = 52.999, p-value = 0.7712
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -1.2271893  0.9150216
sample estimates:
mean of x mean of y 
 2.352704  2.508788 
```

Nothing for Caudate either... well, valiant effort.
