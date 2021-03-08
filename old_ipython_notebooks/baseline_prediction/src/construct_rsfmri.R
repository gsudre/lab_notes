gf_fname = '~/data/baseline_prediction/fmri/good_scans_01312017.csv'
gf = read.csv(gf_fname)
nscans = dim(gf)[1]
data = matrix(nrow=nscans, ncol=2290870)
for (i in 1:nscans) {
  cat(sprintf('scan %d / %d\n', i, nscans))
  a = read.table(sprintf("%s/%04d_maskAve.1D", data_dir, gf$Mask.ID[i]))
  b = cor(a)
  data[i, ] = b[upper.tri(b)]
}
maskids = gf$Mask.ID
save(maskids, data, compress=T,
     file='~/data/baseline_prediction/fmri/rois_spheres/spheres_corr.RData')