# 2020-11-09 20:27:15

Since we had this long break because of COVID, all Freesurfer and DTI data
finished being processed. It makes sense to create a big csv with all the data
in it.

## MPRAGE

We start from Labmatrix to get the best mprage QC, SID, and age_scan. We also
make sure to only get the scans that were processed.

Then:

```bash
# ncrshell01
bash ~/research_code/lab_mgmt/get_freesurfer_roi_data.sh ~/tmp/ids10.txt
```

Then we organize everything in Excel, and to create the summary metrics we do:

