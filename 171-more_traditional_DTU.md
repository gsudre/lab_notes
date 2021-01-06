# 2020-12-31 15:19:27

Let's try a more traditional DTU analysys, like the one here:

http://bioconductor.org/packages/release/workflows/vignettes/rnaseqDTU/inst/doc/rnaseqDTU.html

I'm having trouble importing David's data... not quite sure what format it's in.
Why don't I just go ahead and re-run the quantification step? I can use salmon,
like in the workflow above, or use rsem and kallisto like this:

https://ycl6.gitbook.io/guide-to-rna-seq-analysis/raw-read-processing/mapping

Or maybe spend some more time trying to do the tximport using none type? I did
some research and id4, the item right before "protein_coding" type of column, is
the transcript length in base pairs unit. ID1 is the transcript id. 

# 2021-01-04 09:22:57

Can we find some example RSEM TPM output?

https://github.com/bli25broad/RSEM_tutorial

So, it's just a matter of extracting some of that information. Let's do a single
sample.

# 2021-01-05 11:13:29

None of that was working, so I asked David and Derek for the data. Let's run
this for RSEM first, then we can try Kallisto. Finally, if none of that works,
we can try Salmon, using the raw fastq files, as Derek saud he didn't do any
trimming.

```r
# from David
df = read.delim('~/data/isoforms/shaw_adhd.rsem_output.tpm.tsv')
# from Derek
cnts = read.delim('~/data/isoforms/isoform_counts.txt')
a = lapply(df[,1], function(x) strsplit(as.character(x), split="\\|"))
meta_iso = t(data.frame(a))
colnames(meta_iso) = c('id1', 'ensembleID', 'id2', 'id3', 'iso_name',
                        'hgnc_symbol','id4', 'read_type')
data_iso = df[, 2:ncol(df)]

subjs = colnames(data_iso)
fnames = c()
for (s in subjs) {
    print(s)
    sdf = cbind(meta_iso, data_iso[, s])
    colnames(sdf)[ncol(sdf)] = 'tpm'
    if (s %in% colnames(cnts)) {
        m = merge(sdf, cnts[, c('sample_id', s)], by.x='id1', by.y='sample_id')
        colnames(m)[ncol(m)] = 'count'
        fname = sprintf('~/tmp/%s.tsv', s)
        write.table(m, file=fname, sep='\t', row.names=F, quote=F)
        fnames = c(fnames, fname)
    } else {
        cat('not found\n')
    }
}
txdf = data.frame(TXNAME=meta_iso[,'id1'], GENEID=meta_iso[,'ensembleID'])
txdf = txdf[order(txdf$GENEID),] 
rownames(txdf) = NULL
rsem = tximport(fnames, type='none', geneIdCol = 'ensembleID', txIdCol='id1',
             abundanceCol = 'tpm', lengthCol='id4', countsCol='count',
             txOut=T, importer=read.delim, countsFromAbundance="dtuScaledTPM",
             tx2gene=txdf)
save(rsem, subjs, file='~/data/isoforms/tximport_rsem_DTU.RData')
```

Now, let's go back to the analysis. I'm actually using this as the basis:

https://ycl6.gitbook.io/guide-to-rna-seq-analysis/differential-expression-analysis/differential-transcript-usage/dtu-using-drimseq

