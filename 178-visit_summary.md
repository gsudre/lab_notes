# 2021-01-14 11:06:53

Let's create a script that summarizes what each person did in their last visit.
For that, download all applicable custom data forms from Labmatrix. We won't
process anything that has the _old prefix, so that we leave room to clean up the
Custom data forms if needed.

```r
forms_dir = '~/data/labmatrix'

library(gdata)
subjects_file = sprintf('%s/subjects.xlsx', forms_dir)
subjs = read.xls(subjects_file)
res = data.frame(SID = subjs$Individual.ID...Subjects,
                 MRN = subjs$Medical.Record...MRN...Subjects,
                 Name = subjs$Full.Name...Subjects,
                 nuclearFamID = subjs$Nuclear.ID...FamilyIDs,
                 extendedFamID = subjs$Extended.ID...FamilyIDs,
                 DOB = subjs$Date.of.Birth...Subjects,
                 lastAge = NA)
res = res[!is.na(res$MRN), ]
res = res[res$MRN > 0, ]
res = res[res$DOB != '', ]

forms = list.files(path = forms_dir, pattern = 'xlsx$')
forms = forms[!grepl(pattern='subject', forms)]
forms = forms[!grepl(pattern='_old', forms)]

# load each form
for_names = c()
for (f in forms) {
    cat(f, '\n')
    form_name = gsub(f, pattern='.xlsx', replacement='')
    tmp = read.xls(sprintf('%s/%s', forms_dir, f))
    dates = rep('', length=nrow(res))
    for (s in 1:nrow(res)) {
        # grab the subject's data in this form and compute the last one
        sdata = tmp[tmp$subject.id == res$SID[s], ]
        if (nrow(sdata) > 0) {
            last_date = max(as.Date(sdata$record.date.collected,
                                    format = '%m/%d/%Y'))
            dates[s] = as.character(last_date)
        }
    }
    res = cbind(res, dates)
    colnames(res)[ncol(res)] = form_name
    form_names = c(form_names, form_name)
}

# compute last age for each subject
for (s in 1:nrow(res)) {
    last_date = which.max(sapply(res[s, form_names],
                           function(x) as.Date(x, format='%Y-%m-%d')))
    age = (as.Date(res[s, form_names[last_date]], format="%Y-%m-%d") -
           as.Date(res[s, 'DOB'], format = '%m/%d/%Y'))
    res[s, 'lastAge'] = as.numeric(age/365.25)
}
write.csv(res, file=sprintf('%s/last_visit.csv', forms_dir), row.names=F)
```

# TODO
 * make sure the final subject number conforms to last CIER