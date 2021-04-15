# 2021-04-15 14:07:22

Let me try some very basic DESeq2 to see if I can get some single gene results.
I think our main results will remain there, because they're quite robust
anyways.

As usual, I'm using this:

http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

But I do want to start from scratch:

```r
# I'm running all these in RStudio because it makes graphics much easier
data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name

library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]
df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))

# make sure we only have 
library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = df,
                              design = ~ Diagnosis)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

library(pcaExplorer)
pcaExplorer(dds = dds)
```

I'm using the log2 of the data for dst. I get this:

![](images/2021-04-15-15-23-25.png)

But I get the same thing for other normalizations.

Let's do this for ACC only then:

```r
idx = df$Region=='ACC'
dds <- DESeqDataSetFromMatrix(countData = data[, idx],
                              colData = df[idx, ],
                              design = ~ Diagnosis)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds = DESeq(dds)
library(pcaExplorer)
pcaExplorer(dds = dds)
```

![](images/2021-04-15-15-31-11.png)

There's still something funky there. Let's see if we can find out.

![](images/2021-04-15-15-33-15.png)

PC1 could be brainbank, but there might be a better one...

![](images/2021-04-15-15-39-45.png)

This is also worrisome. Some samples have way more reads than others... 

![](images/2021-04-15-15-41-33.png)

There you go. It's the same as batch. But it doesn't account for the difference
in the PCA plot, especially PC2:

![](images/2021-04-15-15-44-00.png)

Would SV help here?

First, let's check other methods:

```r
library(Glimma)
glimmaMDS(dds)
```

This just uses the same plot we've been looking at using PCA, but now it's MDS.
But it doesn't look as nice, so I'd go with pcaExplorer anyways. It does make me
see that the second PC is related to Sex:

![](images/2021-04-15-15-48-59.png)

Assuming PCs and MDS align.

The other thing we can do with glima is the MDA plot, but I won't do that now.
Let me see if regionReport shows anything else we should look at:

```r
dir.create('~/tmp/DESeq2-example', showWarnings = FALSE, recursive = TRUE)
library('ggplot2')
library('regionReport')
report <- DESeq2Report(dds, project = 'DESeq2 HTML report',
    intgroup = c('Diagnosis', 'Region'), outdir = '~/tmp/DESeq2-example',
    output = 'index', theme = theme_bw())
```

# TODO
 * make sure we only have one sample for each brain
 * make sure all variables are of correct type
 * see if SVA would help get rid of those main effects in the PCA plot














data = t(data)
data = cbind(as.numeric(sn), data)
colnames(data)[1] = 'submitted_name'
df = read.csv('~/data/rnaseq_derek/UPDATED_file_for_derek_add_cause_of_death.csv')
df = df[!duplicated(df$submitted_name),]

m = merge(df, data, by='submitted_name', all.x=F, all.y=T)
pop_code = read.csv('~/data/rnaseq_derek/file_pop.csv')
m2 = merge(m, pop_code, by='hbcc_brain_id')
pcs = read.table('~/data/rnaseq_derek/HM3_b37mds.mds', header=1)
myids = sapply(1:nrow(pcs), function(x) as.numeric(gsub('BR', '',
                                                        strsplit(as.character(pcs[x,'IID']), '_')[[1]][1])))
pcs$numids = myids
m3 = merge(m2, pcs, by.x='hbcc_brain_id', by.y='numids', all.x=T, all.y=F)
```

Choosing which of the repeated sample in Caudate to keep:

```r
cdata = m3[m3$Region=='Caudate',]
grex_names = colnames(m3)[grepl(colnames(m3), pattern='^ENS')]
n0 = rowSums(cdata[, grex_names]==0)
idx = which(cdata$hbcc_brain_id==2877)
idx

[1] 25 26 27 28 29 30

n0[idx]
   49    50    51    52    53    54 
26733 26301 26410 26070 26502 24469 
```

So, row 54 has the fewest genes with zero transcription counts. We'll take that.

```r
cdata = cdata[-c(25:29), ]
m4 = rbind(cdata, m3[m3$Region=='ACC',])
saveRDS(m4, file='~/data/rnaseq_derek/complete_rawCountData_05132020.rds')
```
