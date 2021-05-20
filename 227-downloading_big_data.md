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

# 2021-05-18 16:28:16

For fmri, we can just use the files Adam's team have already processed in
/data/ABCD_MBDU/abcd_bids/derivatives/fmriprep/fmriprep_20.2.0. But they didn't
download the DTI data, so I'll need to do that.

```bash
# ncrshell01
cd /mnt/NCR/sudregp/abcdbids/nda-abcd-s3-downloader
conda activate NDA

./download.py -i ../datastructure_manifest.txt -o /mnt/NCR/sudregp/ABCD/BIDS/ \
    -s download_batch1.txt -d dwi_subsets.txt -p 4

# some subjects don't have dwi data, so let's just figure out who they are
# to update the spreadsheet, and remove them from the folder to save space
rm -rf no_dwi.txt;
for s in `cat download_batch1.txt`; do
    if [[ ! -d /mnt/NCR/sudregp/ABCD/BIDS/${s}/ses-baselineYear1Arm1/dwi ]]; then
        echo $s >> no_dwi.txt;
    fi;
done

for s in `cat no_dwi.txt`; do
    rm -rf /mnt/NCR/sudregp/ABCD/BIDS/${s};
done

# remove the extra fmaps for BOLD that we won't need
rm /mnt/NCR/sudregp/ABCD/BIDS/*/ses-baselineYear1Arm1/fmap/*run*
```

# useful links
 * https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0226715
 * https://onlinelibrary.wiley.com/doi/full/10.1002/hbm.24691