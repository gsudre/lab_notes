# assumes ldata has been processed and we set root_dir and root_fname

# open first results file to figure out model names
var_names = load(sprintf('%s/%s_split%02d.RData', root_dir, root_fname, 1))

# figure out model names and create appropriate collection files
run_models = c()
for (v in var_names) {
  if (grepl('Fit$', v)) {
    run_models = c(run_models, sub('Fit', '', v))
  }
}

rndForestAll = c()
lrAll = c()
lsvmAll = c()
rsvmAll = c()
xgbAll = c()
gbmAll = c()
save_list = c('train_idx')
for (m in run_models) {
  save_list = c(save_list, sprintf('%sFit', m))
}

files = list.files(path=root_dir, pattern=sprintf('%s_split*', root_fname))
# for each results file
for (f in files) {
  load(sprintf('%s/%s', root_dir, f))
  Xtrain = as.data.frame(ldata[train_idx, ])
  colnames(Xtrain) = colnames(ldata)
  ytrain = groups[train_idx]
  Xtest = as.data.frame(ldata[-train_idx, ])
  colnames(Xtest) = colnames(ldata)
  ytest = groups[-train_idx]
  
  for (m in run_models) {
    eval(parse(text=sprintf('%sRes = eval_model(%sFit, Xtest, ytest)', m, m)))
    eval(parse(text=sprintf('%sAll = rbind(%sAll, %sRes)', m, m, m)))
  }
}