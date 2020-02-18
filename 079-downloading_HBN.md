# 2020-02-14 18:31:36

Started downloading the MRI data (and eventually EEG) for the Healthy Brain
Network consortium. 

http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/index.html

They store everything into S3 buckets, so I'll need something like CyberDuck to
download it directly. I'm putting it all in NCR/HBN. 

All the phenotype files were already downloaded.

I also looked at the derivatives folder, but it looks it's only freesurfer and
ants. That could be helpful, but for now I won't download it. We can always
calculate it later through our own pre-processing.

Still need to decide if we're downloading the Staten Island cohort or not, as
it's 1.5T:
http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/MRI%20Protocol.html

But CyberDuck is crapping out after a few downloads... it just stalls. Let's try
to script it out.

```bash
brew install duck

duck --username anonymous --list s3:/fcp-indi/data/Archives/HBN/MRI/Site-RU | grep "sub-" > ~/tmp/files.txt
cd /Volumes/NCR/HBN/MRI/Site-RU/
for f in `cat ~/tmp/files.txt`; do
    wget https://fcp-indi.s3.amazonaws.com/data/Archives/HBN/MRI/Site-RU/$f;
done
```

And just do that for each site.

# 2020-02-18 09:47:36

Philip said we don't need the SI data (1.5T), and the CUNY data is only 25
participants. Since I don't know the scanner there, I'll also leave it off. So,
byt the end of this we'll have 1227 paerticipants from RU and 905 from CBIC.

For the phenotyping data, the download is a bit different:

http://fcon_1000.projects.nitrc.org/indi/cmi_healthy_brain_network/Pheno_Access.html

There also seems ot be a NKI Enhanced Rockland sample:

http://fcon_1000.projects.nitrc.org/indi/enhanced/neurodata.html

I'll grab that (probably the BIDS converted one), when I'm done with HBN. But
focus on the kids only, which should be something close to 300:

https://s3.amazonaws.com/fcp-indi/data/Projects/RocklandSample/RawDataBIDS/participants.tsv

