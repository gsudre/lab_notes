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
forms = forms[!grepl(pattern='special', forms)]

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

# create column for whether the subject has been genotyped or sequenced
cat('Genotyped or sequenced?\n')
for (md in c('SNP', 'WES')) {
    tmp = read.xls(sprintf('%s/special_%s.xlsx', forms_dir, md))
    vals = rep('', length=nrow(res))
    for (s in 1:nrow(res)) {
        vals[s] = res$SID[s] %in% tmp$subject.id
    }
    res = cbind(res, vals)
    colnames(res)[ncol(res)] = sprintf('has%sdata', md)
}

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

# make sure all consent dates have something on them
subjects_file = sprintf('%s/manual/all_visits_20210115.xlsx', forms_dir)
df = read.xls(subjects_file)
date_vars = colnames(df)[grepl(colnames(df), pattern='^Visit')]
missing_file = '~/tmp/missing.txt'
cat('Checking if we can find the consent dates...\n')
for (s in 1:nrow(df)) {
    sdates = df[s, date_vars]
    sdates = sdates[!is.na(sdates)]
    sdates = sdates[sdates != '']
    sdates = format(as.Date(sdates, format="%Y-%m-%d"), '%m/%d/%Y')
    # try to find these subject dates in all_res
    for (d in sdates) {
        idx = which(all_res$SID == df[s, 'SID'] & all_res$date == d)
        if (length(idx) == 0) {
            cat('Cannot find', df[s, 'SID'], 'on', d, '\n', file=missing_file,
                append=T)
            sdata = all_res[all_res$SID==df[s, 'SID'],]
            sdata = sdata[!duplicated(sdata),]
            sdata = sdata[order(as.Date(sdata$date, format='%m/%d/%Y')),]
            for (i in 1:nrow(sdata)) {
                cat('\t', paste0(sdata[i, ]), '\n', file=missing_file, append=T)
            }
        }
    }
}

# asnwering the question  'where are the biospecimens associated with consents obtained from 2018-2020'

# first, are all biosamples in that range associated with a consent
closest = rep('', length=nrow(res_clean))
time_diff = rep('', length=nrow(res_clean))
for (s in 1:nrow(res_clean)) {
    sdates = df[df$SID==res_clean$SID[s], date_vars]
    sdates = sdates[!is.na(sdates)]
    sdates = sdates[sdates != '']
    # subject consent dates
    sdates = as.Date(sdates, format="%Y-%m-%d")
    diffs = as.Date(res_clean[s, 'biosample'], format='%Y-%m-%d') - sdates
    consents_before = which(diffs >= 0)
    sdates = sdates[consents_before]
    diffs = diffs[consents_before]
    if (length(diffs) > 0) {
        mydiff = which.min(diffs)
        time_diff[s] = diffs[mydiff]
        closest[s] = format(sdates[mydiff], '%m/%d/%Y')
    }
}
res_clean = cbind(res_clean, closest)
colnames(res_clean)[ncol(res_clean)] = 'biosampleConsent'
res_clean = cbind(res_clean, time_diff)
colnames(res_clean)[ncol(res_clean)] = 'biosampleConsentTimeDiff'

write.csv(res_clean, file=sprintf('%s/last_visit.csv', forms_dir), row.names=F)


# spitting out long form
long = c()
for (s in 1:nrow(df)) {
    sdates = df[s, date_vars]
    sdates = sdates[!is.na(sdates)]
    sdates = sdates[sdates != '']
    sdates = format(as.Date(sdates, format="%Y-%m-%d"), '%m/%d/%Y')
    # try to find these subject dates in all_res
    for (d in sdates) {
        long = rbind(long, c(df[s, 'SID'], df[s, 'MRN'], d))
    }
}

# let's make sure every consent obtained in 2018 or later has a biosample
colnames(long) = c('SID', 'MRN', 'date')
consents = as.data.frame(long)
consents = consents[as.Date(consents$date, , format='%m/%d/%Y') > 
                    as.Date('01/01/2018', format='%m/%d/%Y'), ]
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
tmp = read.xls(sprintf('%s/biosample.xlsx', forms_dir))
tmp$Notes = NULL
m = mergeOnClosestDate(consents, tmp, unique(consents$SID), x.id = 'SID',
                       y.id='subject.id', x.date='date',
                       y.date = 'record.date.collected')
write.csv(m, file='~/tmp/2018.csv', row.names=F)

# let's make sure every biosample obtained in 2018 or later has a consent
colnames(long) = c('SID', 'MRN', 'date')
consents = as.data.frame(long)
tmp = read.xls(sprintf('%s/biosample.xlsx', forms_dir))
tmp$Notes = NULL
tmp = tmp[as.Date(tmp$record.date.collected, format='%m/%d/%Y') > 
          as.Date('01/01/2018', format='%m/%d/%Y'), ]
source('~/research_code/lab_mgmt/merge_on_closest_date.R')
m = mergeOnClosestDate(tmp, consents, unique(tmp$subject.id), y.id = 'SID',
                       x.id='subject.id', y.date='date',
                       x.date = 'record.date.collected')
write.csv(m, file='~/tmp/2018b.csv', row.names=F)

```

# TODO
 * What is adult self?
 * Who has been genotyped 
 * Who has been sequenced
 * is subject in all_visits? is all_visits subj in res?
 * Where is the biospecimen acquired >= 2018

Who is in the study, when they signed a consent, and what they did in that day, FSIQ, date when obtained (make sure they match date of visit(, match to participant study membership file
How it was created, and how it was checked 
Biomaterials (make sure date matches what’s in the spreadsheet)

Philip,

Sorry for the delay. There was a lot to do, and I'm afraid still some to be
done. I was able to match a file BTRIS support had sent me last year, with the
dates for all consents (as seen in the Protocols tab in CRIS), to those dates in
the wide form spreadsheet. There were many missing that I had to add by hand
(over 120), and I had to correct some of the dates that were there as well (40
or so). I sent you an updated wide list in a new spreadsheet (conset_dates tab).
There are still a few dates in the wide list that we need to chat about: they
are the six with red comments in columns N:R in tab dates_comparison.

Now that I have a more reliable set of dates I can move on to add the info you
asked about biospecimens associated with consents from 2018 and on. I'll have
time to work on this over the weekend/holiday, so I'll send you a new
spreadsheet soon.

