# 2019-10-18 09:27:36

Here are the steps we're using to clean a new SNP dataset. We basically use the
ENIGMA protocol. 

http://enigma.ini.usc.edu/wp-content/uploads/2012/07/ENIGMA2_1KGP_cookbook_v3.pdf

But because we have multiple datasets, we need to merge them first.

# 2019-10-21 15:56:03

First, we need to install GenomeStudio all over again:

wget
ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/genomestudio/genomestudio-software-v2-0-4-5-installer.zip

Then, install it on virtualbox.

Also get the PLINK export tool:

https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/genomestudio/plink-input-report-plugin-v2-1-4.zip

NOT NEEDED....

I was doing it wrong. So, in the citrix Genome Studio v2.0
(citrix.nhgri.nih.gov), load sample sheet and have the same folder as the data
source. Then, for the manifest, use the Illumina folder that also has the
cluster file!

That should load the project. Then, we can just export everything, and save the
results into the correct folders in the shaw drive.

Note that I ended up using my own virtual machine, after installing the software
in centaur/home$/app_bk, because the PLINK export module isn't set up properly
in citrix, and also not in sciware. I cannot see it at all in Citrix, and I get
a permissions error in Sciware.