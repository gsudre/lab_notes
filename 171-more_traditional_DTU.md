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

```r
df = read.delim('~/data/isoforms/shaw_adhd.rsem_output.tpm.tsv')
a = lapply(df[,1], function(x) strsplit(as.character(x), split="\\|"))
meta_iso = t(data.frame(a))
colnames(meta_iso) = c('id1', 'ensembleID', 'id2', 'id3', 'iso_name',
                        'hgnc_symbol','id4', 'read_type')
data_iso = df[, 2:ncol(df)]