# 2021-04-21 14:15:21

Let's do the PCA figures:

```r
data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = round(data)
sub_name = gsub(x=colnames(data), pattern='X', replacement='')
colnames(data) = sub_name
# # this is a repeat for Caudate hbcc 2877, but has more genes with zeros than
# # its other replicate
# data = data[, ! colnames(data) %in% c('66552')]
# # outliers based on PCA plots
# outliers = c('68080','68096', '68108', '68084', '68082')
# data = data[, ! colnames(data) %in% outliers]

library(gdata)
df = read.xls('~/data/post_mortem/POST_MORTEM_META_DATA_JAN_2021.xlsx')
data = data[, colnames(data) %in% df$submitted_name]
df = df[df$submitted_name %in% colnames(data), ]
df = df[order(df$submitted_name), ]
data = data[, order(df$submitted_name)]

# cleaning up some variables
df$Individual = factor(df$hbcc_brain_id)
df[df$Manner.of.Death=='Suicide (probable)', 'Manner.of.Death'] = 'Suicide'
df[df$Manner.of.Death=='unknown', 'Manner.of.Death'] = 'natural'
df$MoD = factor(df$Manner.of.Death)
df$Sex = factor(df$Sex)
df$batch = factor(df$batch)
df$run_date = factor(gsub(df$run_date, pattern='-', replacement=''))
df$Diagnosis = factor(df$Diagnosis, levels=c('Control', 'Case'))
df$Region = factor(df$Region, levels=c('Caudate', 'ACC'))
df$substance_group = factor(df$substance_group)
df$comorbid_group = factor(df$comorbid_group_update)
df$evidence_level = factor(df$evidence_level)
df$brainbank = factor(df$bainbank)
# replace the one subject missing population PCs by the median of their
# self-declared race and ethnicity
idx = (df$Race.x=='White' & df$Ethnicity.x=='Non-Hispanic' & !is.na(df$C1))
pop_pcs = c('C1', 'C2', 'C3', 'C4', 'C5')
med_pop = apply(df[idx, pop_pcs], 2, median)
df[which(is.na(df$C1)), pop_pcs] = med_pop
df$BBB = factor(sapply(1:nrow(df),
                        function(x) sprintf('%s_%s',
                                    as.character(df[x,'brainbank']),
                                    as.character(df[x, 'batch']))))
df$BBB2 = NA                                                                        
df[df$brainbank=='nimh_hbcc', 'BBB2'] = 1                                           
df[df$batch==3, 'BBB2'] = 2                                                         
df[df$batch==4, 'BBB2'] = 3      
df$BBB2 = factor(df$BBB2)                                                   
                    
library(GenomicFeatures)
txdb <- loadDb('~/data/post_mortem/Homo_sapies.GRCh38.97.sqlite')
txdf <- select(txdb, keys(txdb, "GENEID"), columns=c('GENEID','TXCHROM'),
            "GENEID")
bt = read.csv('~/data/post_mortem/Homo_sapiens.GRCh38.97_biotypes.csv')
bt_slim = bt[, c('gene_id', 'gene_biotype')]
bt_slim = bt_slim[!duplicated(bt_slim),]
txdf = merge(txdf, bt_slim, by.x='GENEID', by.y='gene_id')
tx_meta = data.frame(GENEID = substr(rownames(data), 1, 15))
tx_meta = merge(tx_meta, txdf, by='GENEID', sort=F)
imautosome = which(tx_meta$TXCHROM != 'X' &
                tx_meta$TXCHROM != 'Y' &
                tx_meta$TXCHROM != 'MT')
data = data[imautosome, ]
tx_meta = tx_meta[imautosome, ]

library("DESeq2")
fm_str = '~ BBB2 + Diagnosis'
dds <- DESeqDataSetFromMatrix(countData = data,
                                colData = df,
                                design = as.formula(fm_str))
library(pcaExplorer)
pcaExplorer(dds)
```

Grabbed the PCa figures I needed. Let's grab the data for Table 2:

```
r$> load('~/data/post_mortem/basic_DGE_04202021.RData')                                 

r$> dim(dds.ACC)                                                                        
[1] 23675    53

r$> dim(dds.Caudate)                                                                    
[1] 24059    56

r$> df = rbind(colData(dds.ACC), colData(dds.Caudate))

r$> chisq.test(table(df$Region, df$Diagnosis))                                          

        Pearson's Chi-squared test with Yates' continuity correction

data:  table(df$Region, df$Diagnosis)
X-squared = 0.0038412, df = 1, p-value = 0.9506


r$> dim(df)                                                                             
[1] 109  61

r$> idx = df$Region=='ACC'                                                              

r$> mean(df[idx,'RINe'])                                                                
[1] 5.330189

r$> sd(df[idx,'RINe'])                                                                  
[1] 0.7556416

r$> mean(df[idx,'pcnt_optical_duplicates'])                                             
[1] 3.33

r$> sd(df[idx,'pcnt_optical_duplicates'])                                               
[1] 1.644002

r$> mean(df[!idx,'RINe'])                                                               
[1] 5.910714

r$> sd(df[!idx,'RINe'])                                                                 
[1] 0.5895078

r$> mean(df[!idx,'pcnt_optical_duplicates'])                                            
[1] 2.960179

r$> sd(df[!idx,'pcnt_optical_duplicates'])                                              
[1] 1.307745

r$> mean(df[!idx,'clusters'])                                                           
[1] 113772537

r$> sd(df[!idx,'clusters'])                                                             
[1] 33776915

r$> mean(df[idx,'clusters'])                                                            
[1] 108093178

r$> sd(df[idx,'clusters'])                                                              
[1] 30783939

r$> x='RINe'; t.test(df[idx & df$Diagnosis!='Control', x], df[idx&df$Diagnosis=='Control', x])
                                                                                              
l
         Welch Two Sample t-test

data:  df[idx & df$Diagnosis != "Control", x] and df[idx & df$Diagnosis == "Control", x]
t = -1.2432, df = 46.513, p-value = 0.22
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.6828642  0.1613125
sample estimates:
mean of x mean of y 
 5.187500  5.448276 


r$> x='RINe'; t.test(df[!idx&df$Diagnosis!='Control',x], df[!idx&df$Diagnosis=='Control',x])  

        Welch Two Sample t-test

data:  df[!idx & df$Diagnosis != "Control", x] and df[!idx & df$Diagnosis == "Control", x]
t = 0.87704, df = 53.999, p-value = 0.3844
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.1728015  0.4415515
sample estimates:
mean of x mean of y 
 5.987500  5.853125 

r$> x='pcnt_optical_duplicates'; t.test(df[!idx&df$Diagnosis!='Control',x], df[!idx&df$Diagnos
    is=='Control',x])                                                                         

        Welch Two Sample t-test

data:  df[!idx & df$Diagnosis != "Control", x] and df[!idx & df$Diagnosis == "Control", x]
t = 0.65519, df = 53.993, p-value = 0.5151
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.4604934  0.9075767
sample estimates:
mean of x mean of y 
 3.087917  2.864375 


r$> x='pcnt_optical_duplicates'; t.test(df[idx&df$Diagnosis!='Control',x], df[idx&df$Diagnosis
    =='Control',x])                                                                           

        Welch Two Sample t-test

data:  df[idx & df$Diagnosis != "Control", x] and df[idx & df$Diagnosis == "Control", x]
t = 1.5047, df = 50.385, p-value = 0.1387
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.2239922  1.5626991
sample estimates:
mean of x mean of y 
 3.696250  3.026897 
```

## Developmental sets

```r
r = 'Caudate'

df = read.table(sprintf('~/data/post_mortem/Project_WG26_DGE_%s_BBB2_%s_developmental_10K/enrichment_results_WG26_DGE_%s_BBB2_%s_developmental_10K.txt',
                          r, tolower(r), r, tolower(r)),
                  header=1, sep='\t')[, 1:6]
df[df$link=='overlap', 'link'] = 'pan-developmental'
df$link = as.character(df$link)
df$link = factor(df$link, levels=c("pan-developmental", "prenatal",
                                   "infant (0-2 yrs)", "child (3-11 yrs)",
                                   "adolescent (12-19 yrs)", "adult (>19 yrs)"))
df$Direction = ifelse(df$normalizedEnrichmentScore > 0, 'up', 'down')
df$str = sapply(df$pValue, function(x) sprintf('p = %.1e', x))
df$star_pos = abs(df$normalizedEnrichmentScore) + .1
df[df$pValue ==0, 'str'] = 'p < 1e-5'

star_me = as.character(df[df$FDR < .05, 'link'])
stars.df <- df[df$link %in% star_me, c('link', 'star_pos')]

quartz()
library(ggplot2)
ggplot(data=df, aes(x=link, y=abs(normalizedEnrichmentScore), fill=Direction)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=str), vjust=1.6, size=5, color='white')+
  theme_minimal() + ylab('Absolute Normalized Enrichment Score') + xlab('') +
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=14)) +
  geom_text(data = stars.df, aes(y = star_pos, fill=NA), label = "*",
            size = 10) + ggtitle(r)
```

![](images/2021-04-21-17-42-06.png)


## ACC Mollecular Functions

```r
r = 'ACC'
db = 'geneontology_Molecular_Function_noRedundant'

df = read.table(sprintf('~/data/post_mortem/Project_WG26_DGE_%s_BBB2_%s_10K/enrichment_results_WG26_DGE_%s_BBB2_%s_10K.txt',
                          r, db, r, db),
                  header=1, sep='\t')[, 1:7]
df = df[df$FDR < .05, c('description', 'normalizedEnrichmentScore')]
df$Direction = ifelse(df$normalizedEnrichmentScore > 0, 'up', 'down')
df$star_pos = abs(df$normalizedEnrichmentScore)
my_order = order(df$star_pos)
df$description = as.character(df$description)
# organizing based on separate groups (for now, don't do anything)
mylevels = c("respiratory burst",
             "pri-miRNA transcription by RNA polymerase II",
             "microtubule organizing center organization",
             "regulation of synapse structure or activity",
             "regulation of trans-synaptic signaling",
             "chemical synaptic transmission, postsynaptic",
             "synapse organization",
             "vascular endothelial growth factor receptor signaling pathway",
             "gamma-aminobutyric acid signaling pathway",
             "serotonin receptor signaling pathway"
             )
mylevels = df$description
df$description = factor(df$description, levels=mylevels)

library(ggplot2)
ggplot(data=df, aes(x=description, y=abs(normalizedEnrichmentScore), fill=Direction)) +
  geom_bar(stat="identity")+ coord_flip() +
  theme_minimal() + ylab('Absolute Normalized Enrichment Score') + xlab('') +
  theme(axis.text.x = element_text(size=16),
        axis.title.y = element_text(size=16),
        axis.text.y = element_text(size=12),
        legend.text = element_text(size=14))
```

![](images/2021-04-21-17-59-19.png)

# Specific genes of a given pathway

```r
r = 'ACC'
db = 'geneontology_Molecular_Function_noRedundant'

df = read.table(sprintf('~/data/post_mortem/Project_WG26_DGE_%s_BBB2_%s_10K/enrichment_results_WG26_DGE_%s_BBB2_%s_10K.txt',
                          r, db, r, db),
                  header=1, sep='\t')
genes = df[df$description=='GABA receptor activity', 'userId']
gene_list = strsplit(genes, ';')[[1]]

load('~/data/post_mortem/basic_DGE_04202021.RData')
res_str = sprintf('dds = dds.%s', r)
eval(parse(text=res_str))

vsd <- vst(dds, blind=FALSE)
norm.cts <- assay(vsd)
mat <- limma::removeBatchEffect(norm.cts, vsd$BBB2)

gnames = data.frame(full=rownames(counts(dds)),
                    nov=substr(rownames(counts(dds)), 1, 15))
mart = readRDS('~/data/rnaseq_derek/mart_rnaseq.rds')
gnames = merge(gnames, mart, by.x='nov', by.y='ensembl_gene_id')
keep_me = gnames$nov %in% gene_list
gene_ids = gnames[keep_me, ]

resid_expr = reshape2::melt(mat[gene_ids$full,])
colnames(resid_expr) = c('gene', 'submitted_name', 'normCount')
junk = colData(vsd)[, c('Diagnosis', 'submitted_name')]
resid_expr = merge(resid_expr, junk, by='submitted_name')
resid_expr = merge(resid_expr, gene_ids, by.x='gene', by.y='full')

# plotting each of the significant genes
library(ggpubr)
library(ggbeeswarm)
quartz()
myplots = list()
clrs = c("green3", "red")
for (g in 1:nrow(gene_ids)) {
    cat(gene_ids[g, 'nov'], '\n')
    d = as.data.frame(resid_expr[resid_expr$nov == gene_list[g],])
    p = (ggplot(d, aes(x=Diagnosis, y=normCount, color = Diagnosis,
                    fill = Diagnosis)) + 
        scale_y_log10() +
        geom_boxplot(alpha = 0.4, outlier.shape = NA, width = 0.8,
                    lwd = 0.5) +
        stat_summary(fun = mean, geom = "point", color = "black",
                    shape = 5, size = 3,
                    position=position_dodge(width = 0.8)) +
        scale_color_manual(values = clrs) +
        scale_fill_manual(values = clrs) +
        geom_quasirandom(color = "black", size = 1, dodge.width = 0.8) +
        theme_bw() + #theme(legend.position = "none") + 
        ggtitle(gene_ids[g, 'hgnc_symbol']))
    myplots[[g]] = p
}
p = ggarrange(plotlist=myplots)
print(p)
```

![](images/2021-04-21-18-12-07.png)
