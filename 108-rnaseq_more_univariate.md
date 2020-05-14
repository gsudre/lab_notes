# 2020-05-13 07:01:37

Going back to the items in note 106, let's try to do the univariate analysis
again, but removing the bad genes at different thresholds.

Philip is interested in the DX*Region interaction, so the most efficient thing
here will be to run it for all genes, and then start removing genes as we go.
I'll add the population variables, along with the ones that were significant
from note 107.

```r
library(nlme)
pthresh = .1

data = readRDS('~/data/rnaseq_derek/complete_data_04292020.rds')
data = data[-c(which(rownames(data)=='57')), ] # removing ACC outlier
more = readRDS('~/data/rnaseq_derek/data_from_philip_POP_and_PCs.rds')
more = more[, 1:33]
more = more[!duplicated(more$hbcc_brain_id), ]
data = merge(data, more, by='hbcc_brain_id', all.x=F, all.y=F)

data$hbcc_brain_id = as.factor(data$hbcc_brain_id)

# dependent
dep_vars = colnames(data)[grepl(colnames(data), pattern='^ENS')]
# variables to be tested/screened
test_vars = c(# brain-related
              "bainbank.x", 'PMI.x',
              # technical
              'RINe.x', 'run_date.x',
              #clinical
              'comorbid_group', 'substance_group',
              # others
              'Sex.x', 'Age.x',
              sapply(1:10, function(x) sprintf('C%d', x)))
# spit out the results
out_fname = sprintf('~/data/rnaseq_derek/univar_pLT%.02f_DXRegion.csv',
                    pthresh)

hold = c()
for (dp in 1:length(dep_vars)) {
    if (dp %% 50 == 0) {
        print(sprintf('%d of %d (%s)', dp, length(dep_vars), out_fname))
    }

    dep_var = dep_vars[dp]
    fm_str = paste(dep_var, ' ~ Diagnosis.x*Region.x +',
                   paste(test_vars, collapse='+'), sep="")
    fit = try(lme(as.formula(fm_str), ~1|hbcc_brain_id, data=data, na.action=na.omit))
    if (length(fit) > 1) {
        res = summary(fit)$tTable
        # filtering variables
        sig_vars = c()
        for (v in 1:length(test_vars)) {
            # rows in results table that correspond to the screened variable
            var_rows = which(grepl(rownames(res),
                            pattern=sprintf('^%s', test_vars[v])))
            for (r in var_rows) {
                if (res[r, 'p-value'] < pthresh) {
                    sig_vars = c(sig_vars, test_vars[v])
                }
            }
        }
        # factors might get added several times, so here we clean it up
        sig_vars = unique(sig_vars)
        if (length(sig_vars) > 0) {
            clean_fm_str = paste(dep_var, ' ~ Diagnosis.x*Region.x + ',
                                 paste(sig_vars, collapse='+'), sep="")
        } else {
            clean_fm_str = paste(dep_var, ' ~ Diagnosis.x*Region.x', sep="")
        }
        # new model
        clean_fit = try(lme(as.formula(clean_fm_str), ~1|hbcc_brain_id, data=data,
                        na.action=na.omit))
        if (length(clean_fit) > 1) {
            res = data.frame(summary(clean_fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = clean_fm_str
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        } else {
            res = data.frame(summary(fit)$tTable)
            # remove intercept
            res = res[2:nrow(res),]
            res$dep_var = dep_var
            res$formula = NA
            res$orig_formula = fm_str
            res$predictor = rownames(res)
        }
        hold = rbind(hold, res)
    }
}
write.csv(hold, file=out_fname, row.names=F)
```

At a first look I'm noticing a whole bunch of genes where the brain bank is
highly significant. It could be that one of them having so few is the culprit
but it'd still be nice to correct for that. Time for combat? Let's at least see
if combat correct for that, so we dont' need to keep adding it as a covariate in
the models.

## Raw counts

I asked Derek for the raw counts so I can use Combat. Let's give it a try:

```r
df = read.csv('~/data/rnaseq_derek/UPDATED_file_for_derek_add_cause_of_death.csv')
df = df[!duplicated(df$submitted_name),]
data = read.table('~/data/rnaseq_derek/adhd_rnaseq_counts.txt', header=1)
rownames(data) = data[,1]
data[,1] = NULL
data = t(data)
sn = gsub(x=rownames(data), pattern='X', replacement='')
data = cbind(as.numeric(sn), data)
colnames(data)[1] = 'submitted_name'
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

Now we can try correcting for batch effects.

Reading about sva, combat, and combat-seq
(https://www.bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf),
I'm not sure whether I should try removing just the batch or the site as well.
Maybe do batch first, and then try to estimate if there are any leftovers?

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
```

This is an incredible resource:

http://www.nathalievialaneix.eu/doc/pdf/tutorial-rnaseq.pdf

these might also become useful:

* https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html
* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4937821/
* https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4934518/
  
# TODO:
 * does the distribution shape matter? what if we inverse transform it? should
   we do it within region?