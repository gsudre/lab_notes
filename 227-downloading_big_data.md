# 2021-05-12 05:59:35

I was chatting with Luke yesterday, and for ABCd the approach will be to
download their data already in BIDS format from the ABCD-BIDS repository
(https://collection3165.readthedocs.io/en/stable/). Their pipeline for the
derivatives is not necessarily clear, and it doesn't provide the voxel-wise data
we might eventually need in MNI space. The only data they give in MNI space is
not filtered. It's also questionable whether their filtering is the best.

Let's download the source BIDS files and use our own pipeline. This way the
pipeline will be constant across datasets, and we can always download more ABCD
data (from new releases) ahead of ABCD-BIDS, if needed.

So, the datasets I'll need are:

inputs.anat.T1w
inputs.anat.T2w
inputs.dwi.dwi
inputs.func.task-rest

I'm dowloading them in ncrshell01, directly to NCR drive, to NCR/ABCD/BIDS. I'll
download a first batch of 1000 subjects, but the tracker is in Teams. Then we
can proceed with fmriprep and xcpengine.

```bash
# ncrshell01
cd /mnt/NCR/sudregp/abcdbids/nda-abcd-s3-downloader
conda activate NDA
./download.py -i ../datastructure_manifest.txt -o /mnt/NCR/sudregp/ABCD/BIDS/ -s download_batch1.txt -d sources.txt -p 4
```

The download instructions mainly come from here:
https://github.com/ABCD-STUDY/nda-abcd-s3-downloader

