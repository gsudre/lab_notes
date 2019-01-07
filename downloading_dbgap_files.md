# 2019-01-07 16:41:47

Quick note on how I did it. 

* Go to the dbGap project page and create a new download
  (https://dbgap.ncbi.nlm.nih.gov/aa/wga.cgi?page=pi_requests) 
* In the Downloads page, download the repository key.
* Use the Aspera command line to download the files, because it's easier to
  monitor remotely. I'll change depending on the package being downloaded, but
  the gist of it is:

```
ascp -QTr -l 300M -k 1 -i ~/Applications/Aspera\ CLI/etc/asperaweb_id_dsa.openssh -W A1FB65FC09DF84E48A4C49861EDD08D28BF293C754DF64AB2B532A7AD331B6F5E92DD738A25CF4CDB932D9ECF829526F1D dbtest@gap-upload.ncbi.nlm.nih.gov:data/instant/shawp/66134 /Volumes/Shaw/pgc_from_dgbap/
```

* The downloaded files are encrypted, so we need to decrypt them. So, go to the
  bin directory in the SRA toolkit folder you downloaded from
  https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc&f=std, and
  run ./vdb-config -i
* Import the key downloaded above.
* Go to the workspace directory created when the key was imported.
* **From inside that working directory**, run:
  
```
~/sratoolkit.2.9.2-mac64/bin/vdb-decrypt /Volumes/Shaw/pgc_from_dgbap/66133
```

which will decrypt all files in place.

These pages were useful:

https://www.ncbi.nlm.nih.gov/books/NBK63512/