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

# 2021-02-02 07:30:59

Let's try the raw script:

```python
import numpy as np
from scipy.interpolate import splrep, splev
from mne.filter import filter_data, resample
from scipy.signal import detrend, find_peaks

# Load and preprocess data
fname = '/Volumes/NCR/MR_data_by_maskid/2650/E29010/ECG_epiRTnih_scan_0010_20190809_185538.1D'
fid = open(fname)
ecg = [float(line.rstrip()) for line in fid]
fid.close()
sf_ori = 16
sf = 100

dsf = sf / sf_ori
ecg = resample(ecg, dsf)
ecg = filter_data(ecg, sf, 2, 30, verbose=0)

# Select only a 20 sec window
window = 20
start = 155
ecg = ecg[int(start*sf):int((start+window)*sf)]

# R-R peaks detection
rr, _ = find_peaks(ecg, distance=40, height=0.5)

# R-R interval in ms
rr = (rr / sf) * 1000
rri = np.diff(rr)

# Interpolate and compute HR
def interp_cubic_spline(rri, sf_up=4):
    """
    Interpolate R-R intervals using cubic spline.
    Taken from the `hrv` python package by Rhenan Bartels.
    
    Parameters
    ----------
    rri : np.array
        R-R peak interval (in ms)
    sf_up : float
        Upsampling frequency.
    
    Returns
    -------
    rri_interp : np.array
        Upsampled/interpolated R-R peak interval array
    """
    rri_time = np.cumsum(rri) / 1000.0
    time_rri = rri_time - rri_time[0]
    time_rri_interp = np.arange(0, time_rri[-1], 1 / float(sf_up))
    tck = splrep(time_rri, rri, s=0)
    rri_interp = splev(time_rri_interp, tck, der=0)
    return rri_interp

sf_up = 4
rri_interp = interp_cubic_spline(rri, sf_up) 
hr = 1000 * (60 / rri_interp)
print('Mean HR: %.2f bpm' % np.mean(hr))

# Detrend and normalize
edr = detrend(hr)
edr = (edr - edr.mean()) / edr.std()

# Find respiratory peaks
resp_peaks, _ = find_peaks(edr, height=0, distance=sf_up)

# Convert to seconds
resp_peaks = resp_peaks
resp_peaks_diff = np.diff(resp_peaks) / sf_up

# Extract the mean respiratory rate over the selected window
mresprate = resp_peaks.size / window
print('Mean respiratory rate: %.2f Hz' % mresprate)
print('Mean respiratory period: %.2f seconds' % (1 / mresprate))
print('Respiration RMS: %.2f seconds' % np.sqrt(np.mean(resp_peaks_diff**2)))
print('Respiration STD: %.2f seconds' % np.std(resp_peaks_diff))
```

This code seems to be working. Now I'll need to extract only the resting scans,
then get their sampling rate from the scan file:

```
RT Physio: ENABLED 1.0
RT Physio: sampling 16 ms
```

And finally spit out everything to a CSV so we can check for outliers, which
signals failed, etc. Probably also good to spit out the length of each signal,
so we know which measurements are reliable or not.