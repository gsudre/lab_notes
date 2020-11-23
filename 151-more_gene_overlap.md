# 2020-11-20 15:04:46

Let's play a bit more with GeneOverlap just using our rnaseq results and our
ready-made gene sets:

```r
load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
library(GeneOverlap)
library(WebGestaltR)
db_file = '~/data/post_mortem/acc_developmental.gmt'
gmt = readGmt(db_file) # already in gene symbols
# and convert it to lists
mylist = list()
for (s in unique(gmt$geneSet)) {
    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
}
t = .001
ra = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
rc = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
gom.obj <- newGOM(list(rnaseq_acc=ra, rnaseq_caudate=rc),
                  mylist, spec='hg19.gene')
getMatrix(gom.obj, name='pval')
```

Have to play a bit more with the thresholds, and potentially add lists with
overlapping genes in them?

# 2020-11-23 06:17:43

I actually don't think it's too bad to use hg19 for this, because we need to
consider the chances of being in either list too.

```r
library(GeneOverlap)
library(WebGestaltR)

load('~/data/rnaseq_derek/rnaseq_results_11122020.rData')
db_file = '~/data/post_mortem/acc_developmental.gmt'
gmt = readGmt(db_file) # already in gene symbols
# and convert it to lists
mylist = list()
for (s in unique(gmt$geneSet)) {
    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
}
for (t in c(.05, .01, .005, .001)) {
    ra = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
    rc = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
    cat(sprintf('\nt=%.3f, ra=%d, rc=%d\n', t, length(ra), length(rc)))
    gom.obj <- newGOM(list(rnaseq_acc=ra, rnaseq_caudate=rc),
                    mylist, spec='hg19.gene')
    print(getMatrix(gom.obj, name='pval'))
}
```

```
t=0.050, ra=1335, rc=1063
               dev1_c0.9 dev2_c0.9   dev3_c0.9  dev4_c0.9  dev5_c0.9 overlap_c0.9
rnaseq_acc     0.6236385 0.8759384 0.020775497 1.00000000 0.09771868    0.9532548
rnaseq_caudate 0.1064240 1.0000000 0.007331972 0.08262743 0.51812821    0.7375957

t=0.010, ra=325, rc=220
               dev1_c0.9 dev2_c0.9 dev3_c0.9 dev4_c0.9 dev5_c0.9 overlap_c0.9
rnaseq_acc     0.8636162         1 0.4154584 1.0000000 0.6844767    0.9824292
rnaseq_caudate 0.5972892         1 0.3041968 0.0461671 0.5412226    0.8290920

t=0.005, ra=165, rc=117
               dev1_c0.9 dev2_c0.9 dev3_c0.9 dev4_c0.9 dev5_c0.9 overlap_c0.9
rnaseq_acc     0.9570919         1         1 1.0000000 0.4421944    0.9857973
rnaseq_caudate 0.3790669         1         1 0.1673051 0.3386872    0.9342330

t=0.001, ra=38, rc=26
               dev1_c0.9 dev2_c0.9 dev3_c0.9 dev4_c0.9 dev5_c0.9 overlap_c0.9
rnaseq_acc             1         1         1         1         1    1.0000000
rnaseq_caudate         1         1         1         1         1    0.6257793
```

Results are not good. We might need to go for enrichment here. Let's just check
adhd_genes though:

```r
db_file = '~/data/post_mortem/adhd_genes.gmt'
gmt = readGmt(db_file) # already in gene symbols
# and convert it to lists
mylist = list()
for (s in c('GWAS1', 'GWAS', 'TWAS1', 'TWAS2', 'TWAS', 'CNV1', 'CNV2')) {
    mylist[[s]] = unique(gmt$gene[gmt$geneSet==s])
}
for (t in c(.05, .01, .005, .001)) {
    ra = rnaseq_acc[rnaseq_acc$P.Value < t, 'hgnc_symbol']
    rc = rnaseq_caudate[rnaseq_caudate$P.Value < t, 'hgnc_symbol']
    cat(sprintf('\nt=%.3f, ra=%d, rc=%d\n', t, length(ra), length(rc)))
    gom.obj <- newGOM(list(rnaseq_acc=ra, rnaseq_caudate=rc),
                    mylist, spec='hg19.gene')
    print(getMatrix(gom.obj, name='pval'))
}
```

```
t=0.050, ra=1335, rc=1063
                   GWAS1      GWAS     TWAS1      TWAS
rnaseq_acc     1.0000000 0.9231465 1.0000000 0.4623828
rnaseq_caudate 0.6770602 0.9506344 0.5426837 0.5996965

t=0.010, ra=325, rc=220
               GWAS1      GWAS TWAS1 TWAS
rnaseq_acc         1 0.9059446     1    1
rnaseq_caudate     1 1.0000000     1    1

t=0.005, ra=165, rc=117
               GWAS1      GWAS TWAS1 TWAS
rnaseq_acc         1 0.6976711     1    1
rnaseq_caudate     1 1.0000000     1    1

t=0.001, ra=38, rc=26
               GWAS1 GWAS TWAS1 TWAS
rnaseq_acc         1    1     1    1
rnaseq_caudate     1    1     1    1
```

Again, nothing to be proud off.