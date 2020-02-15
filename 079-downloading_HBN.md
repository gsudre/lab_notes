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
cd /Vo
for f in `cat ~/tmp/files.txt`; do
    wget https://fcp-indi.s3.amazonaws.com/data/Archives/HBN/MRI/Site-RU/$f;
done
```

And just do that for each site.