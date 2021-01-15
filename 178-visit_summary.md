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

# storing all data for all subjects, for future checks
all_res = c()

# load each form
form_names = c()
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

            junk = data.frame(date = sdata$record.date.collected)
            junk$form = form_name
            junk$SID = res$SID[s]
            all_res = rbind(all_res, junk)
        }
    }
    res = cbind(res, dates)
    colnames(res)[ncol(res)] = form_name
    form_names = c(form_names, form_name)
}

# compute last age for each subject
for (s in 1:nrow(res)) {
    if (! all(res[s, form_names] == '')) {
        last_date = which.max(sapply(res[s, form_names],
                            function(x) as.Date(x, format='%Y-%m-%d')))
        age = (as.Date(res[s, form_names[last_date]], format="%Y-%m-%d") -
            as.Date(res[s, 'DOB'], format = '%m/%d/%Y'))
        res[s, 'lastAge'] = as.numeric(age/365.25)
    }
}

# find date and value for first IQ
cat('Finding last IQs\n')
tmp = read.xls(sprintf('%s/IQ.xlsx', forms_dir))
vals = rep('', length=nrow(res))
for (s in 1:nrow(res)) {
    # grab the subject's data in this form and compute the last one
    sdata = tmp[tmp$subject.id == res$SID[s], ]
    if (nrow(sdata) > 0) {
        last_date = which.max(as.Date(sdata$record.date.collected,
                                       format = '%m/%d/%Y'))
        vals[s] = as.character(sdata[last_date, 'FSIQ'])
    }
}
res = cbind(res, vals)
colnames(res)[ncol(res)] = 'lastIQval'

# now that everything is computed, let's clean it up

# remove anyone that we couldn't compute dates
res_clean = res[!is.na(res$lastAge), ]
# only keep people with at least one entry after 2018
keep_me = c()
for (s in 1:nrow(res_clean)) {
    years = as.numeric(substr(res_clean[s, form_names], 1, 4))
    if (sum(years >= 2018, na.rm=T) > 0) {
        keep_me = c(keep_me, s)
    }
}
res_clean = res_clean[keep_me, ]
write.csv(res_clean, file=sprintf('%s/last_visit.csv', forms_dir), row.names=F)

# make sure all consent dates have something on them
subjects_file = sprintf('%s/manual/all_visits.xlsx', forms_dir)
df = read.xls(subjects_file)
date_vars = colnames(df)[grepl(colnames(df), pattern='^Visit')]
for (s in 1:nrow(df)) {
    sdates = df[s, date_vars]
    sdates = sdates[!is.na(sdates)]
    sdates = sdates[sdates != '']
    sdates = format(as.Date(sdates, format="%Y-%m-%d"), '%m/%d/%Y')
    # try to find these subject dates in all_res
    for (d in sdates) {
        idx = which(all_res$SID == df[s, 'SID'] & all_res$date == d)
        if (length(idx) == 0) {
            cat('Cannot find', df[s, 'SID'], 'on', d, '\n')
        }
    }
}
```

# TODO
 * What is adult self?
 * Who has been genotyped 
 * Who has been sequenced
 * is subject in all_visits? is all_visits subj in res?
 * Where is the biospecimen acquired >= 2018
 * Who has been genotyped 
 * Who has been sequenced


Who is in the study, when they signed a consent, and what they did in that day, FSIQ, date when obtained (make sure they match date of visit(, match to participant study membership file
How it was created, and how it was checked 
Biomaterials (make sure date matches what’s in the spreadsheet)