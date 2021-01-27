# 2021-01-27 13:12:28

Philip asked me to retrieve weight, respiration and heart rate from scans.
Weight is easy, and I could't find height in the DICOM headers. 

## weight

```bash
# ncrshell01
cd /mnt/NCR/sudregp/MR_data_by_maskid
/bin/ls -1 > ~/tmp/maskids.txt;
rm -rf ~/tmp/weights.csv;
for m in `cat ~/tmp/maskids.txt`; do
    echo $m;
    val=`dicom_hdr $m/2*/*01/*0001.dcm | grep Weight`;
    echo $m,$val >> ~/tmp/weights.csv;
done
```

It looks like I can derive RR and HR from ECG:

https://gist.github.com/raphaelvallat/55624e2eb93064ae57098dd96f259611

Or even :

https://python-heart-rate-analysis-toolkit.readthedocs.io/en/latest/quickstart.html#getting-heart-rate-over-time

In the end, Philip needs a spreasheet containing respiration rate and heart rate
as well.