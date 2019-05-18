# 2019-05-15 13:34:09

These are the data for Martine. We couldn't send the matching DTI data to the
structural data Philip had sent her before because that data was mostly 1.5T.
So, let's choose based on DTI and then pick the matching Freesurfer to those. 

Philip said we can send anyone under 35 y.o., but only cross sectional and one
per family. I think the easiest thing will be to get all our DTI and family IDs,
and take the oldest per family, which will hopefully also be the best quality.
Here, I'll consider only one per extended family as well. So, I start with this
search in Labmatrix:

![](images/2019-05-15-13-41-03.png)

Note that I decided to go with the TORTOISE pre-processing for this because
we've used more often in previous work. I'm still playing a bit with the FDT
pipeline, especially the failed scans, so we should hold off in using them for
now.

The first thing I did was to figure out who had bene processed already. Esteban
is still working on some of them, which means these are not necessarily the
oldest or the best we'll have. But at the moment, that's true.

So, we end up with a list of 512 mask ids (have_tortoise tab in
dti_for_enigma.xlsx). Now, let's extract the data ENIGMA wants based on the
protocol in their webpage, and then attach the Freesurfer data to that.

First, I'm following the steps from here:
http://enigma.ini.usc.edu/wp-content/uploads/DTI_Protocols/ENIGMA_TBSS_protocol_USC.pdf

Note that I'm not running DTI-TK, because one of the steps in the ENIGMA steps
will be to transform everything into their templates. This way we can go from
subject native space directly to that.

```bash
#caterpie
maskidDir=/mnt/shaw/MR_data_by_maskid/
cd /mnt/shaw/dti_robust_tsa/enigma
for m in `cat ids512.txt`; do
    tensorDir=${maskidDir}/${m}/edti_proc/edti_DMC_DR_R1_SAVE_DTITK/;
    if [ ! -d $tensorDir ]; then
        tensorDir=${maskidDir}/${m}/edti_proc/edti_DMC_R1_SAVE_DTITK/;
        tensorFile=edti_DMC_R1_tensor.nii;
    else
        tensorFile=edti_DMC_DR_R1_tensor.nii;
    fi;
    if [ ! -d $tensorDir ]; then
        echo $tensorDir "does not exist. Skipping..."
    else
        TVtool -in ${tensorDir}/${tensorFile} -fa -out ./${m}_fa.nii;
    fi;
done
```

# 2019-05-16 10:17:31

It turns out that 3 IDs weren't really processed. So, I'm replacing those ids.
They were:

```
2135
2080
2232
2251
```

```bash
#caterpie
cd /mnt/shaw/dti_robust_tsa/enigma
source activate python2
tbss_1_preproc *.nii
tbss_2_reg -t enigmaDTI/ENIGMA_DTI_FA.nii.gz
tbss_3_postreg -S
```

# TODO

* Update excel spreadsheets with IDs that changed