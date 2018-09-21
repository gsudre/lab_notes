# 2018-09-20 16:00:14

Ayaka asked to run FATCAT in a few IDs for for project. This is how we did it
for DLPFC:

```bash
# caterpie
bin=/usr/local/neuro/TORTOISE_V3.1.2/DIFFPREPV312/bin/bin/ConvertOldListfileToNew;
for m in `cat ~/tmp/no_fatcat.txt`; do     
    echo "Working on $m";
    mydir=/mnt/shaw/data_by_maskID/${m};
    mkdir ${mydir}/tortoise_converted_v2;
    if [ -e ${mydir}/edti_proc/edti_DMC_DR_R1.list ]; then
        ${bin} ${mydir}/edti_proc/edti_DMC_DR_R1.list ${mydir}/edti_proc/edti_DMCstructural.nii ${mydir}/tortoise_converted_v2/${m}.list;
    else
        ${bin} ${mydir}/edti_proc/edti_DMC_R1.list ${mydir}/edti_proc/edti_DMCstructural.nii ${mydir}/tortoise_converted_v2/${m}.list;
    fi;
    cp ${mydir}/edti_proc/edti_DMCtemplate.nii ${mydir}/edti_proc/edti_DMCstructural.nii ${mydir}/tortoise_converted_v2/;
done
```

# 2018-09-21 10:15:00

Some of the IDs Ayaka had sent didn't have good DTI (i.e. were bad according to
Labmatrix and didn't have the necessary directories ot be converted).

NEXT STEPS

Then, we copy the export directory of all subjects to biowulf, and run fatcat:

while read m; do echo "bash ~/research_code/dti/run_fatcat_fs_exported.sh /scratch/sudregp/tortoise_exported_v2/ ${m}" >> swarm.fatcat; done < ~/tmp/400.txt