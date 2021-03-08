# Collect the results form the different tables in Paul Taylor's output

subj_dir = '/Volumes/Labs/Shaw/dti_robust_tsa/paul_taylor/'
subjs = read.table(sprintf('%s/paul_parsed.txt', subj_dir))
data = c()
for (s in subjs[[1]]) {
  print(s)
  subj_data = c()
  files = c('xad', 'xae', 'xaf', 'xag', 'xah', 'xai', 'xaj', 'xak', 'xal', 'xam', 'xan', 'xao',
            'xap', 'xaq', 'xar', 'xas', 'xat')
  titles = read.table(sprintf('%s/%s/parsed_tables/xac', subj_dir, s))[1,]
  for (f in files) {
    fpath = sprintf('%s/%s/parsed_tables/%s', subj_dir, s, f)
    title.line <- readLines(fpath, n=1)
    dtype = substring(title.line, 3)
    a = read.table(fpath)
    b = a[upper.tri(a, diag=TRUE)]
    subj_data = c(subj_data, b)
  }
  data = rbind(data, subj_data)
}
