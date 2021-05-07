# 2021-05-05 19:36:44

Let's look into running MAGMA for other GWAS as well. I'll re-run the ADHD
analysis using the same code as note 202, but this time not in the temporary
folder:

```bash
# BW
conda activate base

cd ~/data/post_mortem
mkdir MAGMA

cd MAGMA
module load plink
echo "g1000_afr.bed g1000_afr.bim g1000_afr.fam" > merge_list.txt;
# they suggest using MAF when subsampling, so let's use it for concatenating too
# I got the g1000 files form their website: https://ctg.cncr.nl/software/magma
plink --bfile g1000_eur --merge-list merge_list.txt --maf 1e-5 \
    --flip-scan --make-bed --out g1000_BW
plink --bfile g1000_eur --exclude g1000_BW-merge.missnp \
    --make-bed --out d1
plink --bfile g1000_afr --exclude g1000_BW-merge.missnp \
    --make-bed --out d2
echo "d2.bed d2.bim d2.fam" > merge_list.txt;
plink --bfile d1 --merge-list merge_list.txt --maf 1e-5 \
    --make-bed --out g1000_BW

module load MAGMA
magma --annotate --snp-loc g1000_BW.bim \
    --gene-loc /usr/local/apps/MAGMA/gene_location/NCBI37.3/NCBI37.3.gene.loc \
    --out annot_BW
magma --bfile g1000_BW --pval ~/pgc2017/adhd_jun2017 N=55374 \
    --gene-annot annot_BW.genes.annot --out genes_BW
for r in 'ACC' 'Caudate'; do
    magma --gene-results genes_BW.genes.raw \
        --gene-covar ../MAGMA_bigger_log10_dge_${r}.tab \
        --out MAGMA_bigger_log10_gc_dge_${r};
done;

magma --annotate --snp-loc g1000_eur.bim \
    --gene-loc /usr/local/apps/MAGMA/gene_location/NCBI37.3/NCBI37.3.gene.loc \
    --out annot_WNH
magma --bfile g1000_eur --pval ~/pgc2017/adhd_eur_jun2017 N=53293 \
    --gene-annot annot_WNH.genes.annot --out genes_WNH
for r in 'ACC' 'Caudate'; do
    magma --gene-results genes_WNH.genes.raw \
        --gene-covar ../MAGMA_bigger_WNH_log10_dge_${r}.tab \
        --out MAGMA_bigger_WNH_log10_gc_dge_${r};
done;
```

Then, we just run the other GWASes we find, keeping in mind that they use
different populations:

```bash
# AAD
# keep only rsid for each SNP
cut -d" " -f 2 pgc_alcdep.discovery.aug2018_release.txt | cut -d":" -f 1 > rsids.txt;
awk '{ print $1,$3,$4,$5,$6,$7,$8}' pgc_alcdep.discovery.aug2018_release.txt > tmp.txt;
paste rsids.txt tmp.txt > aad.txt;
magma --bfile g1000_BW --pval aad.txt N=52848 \
    --gene-annot annot_BW.genes.annot --out AAD_genes_BW
for r in 'ACC' 'Caudate'; do
    magma --gene-results AAD_genes_BW.genes.raw \
        --gene-covar ../MAGMA_bigger_log10_dge_${r}.tab \
        --out MAGMA_AAD_bigger_log10_gc_dge_${r};
done;

# AAD, WNH
cut -d" " -f 2 pgc_alcdep.eur_discovery.aug2018_release.txt | cut -d":" -f 1 > rsids.txt;
awk '{ print $1,$3,$4,$5,$6,$7,$8}' pgc_alcdep.eur_discovery.aug2018_release.txt > tmp.txt;
paste rsids.txt tmp.txt > aad_WNH.txt;
magma --bfile g1000_eur --pval aad_WNH.txt N=46568 \
    --gene-annot annot_WNH.genes.annot --out AAD_genes_WNH
for r in 'ACC' 'Caudate'; do
    magma --gene-results AAD_genes_WNH.genes.raw \
        --gene-covar ../MAGMA_bigger_log10_dge_${r}.tab \
        --out MAGMA_AAD_WNH_bigger_log10_gc_dge_${r};
done;

# ASD is eur only
magma --bfile g1000_eur --pval iPSYCH-PGC_ASD_Nov2017 N=46351 \
    --gene-annot annot_WNH.genes.annot --out ASD_genes
for r in 'ACC' 'Caudate'; do
    magma --gene-results ASD_genes.genes.raw \
        --gene-covar ../MAGMA_bigger_log10_dge_${r}.tab \
        --out MAGMA_ASD_bigger_log10_gc_dge_${r};
done;
for r in 'ACC' 'Caudate'; do
    magma --gene-results ASD_genes.genes.raw \
        --gene-covar ../MAGMA_bigger_log10_dge_${r}.tab \
        --out MAGMA_ASD_WNH_bigger_log10_gc_dge_${r};
done;
```

I'll give a break on these because they take a while to find the right GWAS.
Maybe I'll ask one of the IRTAs to do it. I'll work on the other figures for
now.

