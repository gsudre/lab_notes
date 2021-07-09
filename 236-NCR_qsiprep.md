# 2021-07-01 06:43:31

For our data, we'll have a few extra wrinkles when creating the BIDS directory
for DTI. Talking with Marine (Teams, 5/19/21, 1:50PM), the idea would be to use
the re-acquired volumes when possible, replacing the old ones. Second best would
be using everything. 

So, why not create up to 3 sessions for each mask id? Session 1 is the vanilla
60 or 80, which we'll process similarly to PNC. Session 2 includes the 99 run.
Session 3 has the 99 run replacing the original volumes. Session 3 would only be
created after qsiprep, by removing the re-acquired volumes from the eddy output.
Not sure how much difference this will make, or even if it's feasible, because
it'd only work if all volumes are processed independently.

It was getting a bit too complex for dcm2bids, so I'll need to do it manually:
_
```bash
# ncrshell01
net_dir=/mnt/NCR/sudregp/MR_data_by_maskid;
for m in `cat subjs.txt`; do
    echo $m;
    cd /mnt/shaw/sudregp/NCR_BIDS/sub-${m};

    # find name of date folders
    ls -1 $net_dir/${m}/ | grep -e ^20 > ~/tmp/date_dirs;

    # for each date folder, check for dti scans
    mkdir -p dwi;
    while read d; do
        if grep -q cdiflist08 $net_dir/${m}/${d}/*README*; then
            nruns=3;
            tail -n +2 $net_dir/${m}/E*/cdiflist08 | split \
                --numeric-suffixes=1 -l 20;
        elif grep -q cdiflist09 $net_dir/${m}/${d}/*README*; then
            nruns=4;
            tail -n +2 $net_dir/${m}/E*/cdiflist09 | split \
                --numeric-suffixes=1 -l 20;
        else
            echo "No cdi sequence for $m";
        fi;
        # we need to convert the runs in order
        for s in `seq 1 $nruns`; do
            grep cdiflist $net_dir/${m}/${d}/*README* | \
                grep -i _g0${s} > ~/tmp/scan;
            awk '{for(i=1;i<=NF;i++) {if ($i ~ /Series/) print $i}}' \
                ~/tmp/scan | sed "s/Series://g" > ~/tmp/scan_clean;
            mr_dir=`cat ~/tmp/scan_clean | sed "s/,//g"`;
            echo "Converting ${net_dir}/${m}/${d}/${mr_dir}/";
            dcm2niix_afni -o dwi/ -z y -f sub-${m}_run-${s}_dwi \
                ${net_dir}/${m}/${d}/${mr_dir}/;
            # replace gradients with scaled versions
            cp dwi/sub-${m}_run-${s}_dwi.bval \
                dwi/sub-${m}_run-${s}_dwi.bval.orig;
            cp dwi/sub-${m}_run-${s}_dwi.bvec \
                dwi/sub-${m}_run-${s}_dwi.bvec.orig;
            cp dwi/sub-${m}_run-${s}_dwi.nii.gz \
                dwi/sub-${m}_run-${s}_dwi.nii.gz.orig;
            
            1dDW_Grad_o_Mat++           \
                -in_col_vec        x0${s}          \
                -out_row_vec       x0${s}_rowvec
            cd dwi;
            # disable reorient because qsiprep couldn't merge it
            fat_proc_convert_dcm_dwis \
                -innii sub-${m}_run-${s}_dwi.nii.gz \
                -inbval sub-${m}_run-${s}_dwi.bval.orig \
                -inbvec ../x0${s}_rowvec -prefix sub-${m}_run-${s}_dwi \
                -flip_z -no_qc_view -reorig_reorient_off;
            cp sub-${m}_run-${s}_dwi_bval.dat sub-${m}_run-${s}_dwi.bval;
            cp sub-${m}_run-${s}_dwi_rvec.dat sub-${m}_run-${s}_dwi.bvec;
            # cleaning up
            shopt -s extglob;
            rm -v !(*.bvec|*.bval|*.nii.gz|*.json);
            shopt -u extglob;
            cd ..;

            # # convert gradients to the correct magnitude and format
            # 1dDW_Grad_o_Mat++           \
            #     -in_col_vec        x0${s}          \
            #     -in_bvals          dwi/sub-${m}_run-${s}_dwi.bval.orig  \
            #     -out_row_vec       dwi/sub-${m}_run-${s}_dwi.bvec       \
            #     -out_row_bval_sep  dwi/sub-${m}_run-${s}_dwi.bval        \
            #     -unit_mag_out -flip_z


            # replace PhaseEncodingAxis by PhaseEncodingDirection in JSON to
            # conform with BIDS
            sed -i -e "s/PhaseEncodingAxis/PhaseEncodingDirection/" \
                dwi/sub-${m}_run-${s}_dwi.json;
        done
        rm -f x??;
    done < ~/tmp/date_dirs;
done;
```

# 2021-07-09 13:03:59

I wasn't able to get qsiprep running in the individual sessions, so it'll have
to be done after concatenation. It's not ideal to match PNC, but at least it'd
be easier to deal with the 99 runs. They'd likely not have enough data to run
denoising on their own. Now I can play with processing without them, adding
them, or replacing them.

First, let's check that the qsiprep results look OK.

```bash
cd /lscratch/$SLURM_JOBID
cp -r /lscratch/18740433/qsiprep_output/qsiprep/sub-1233 .
cd sub-1233/dwi
conda activate mrtrix3
dwi2tensor -fslgrad sub-1233_run-1_space-T1w_desc-preproc_dwi.bvec \
    sub-1233_run-1_space-T1w_desc-preproc_dwi.bval \
    -nthreads $SLURM_CPUS_PER_TASK \
    sub-1233_run-1_space-T1w_desc-preproc_dwi.nii.gz tensors.nii.gz
tensor2metric tensors.nii.gz -fa dti_FA.nii.gz -ad dti_AD.nii.gz \
    -rd dti_RD.nii.gz -vector dti_V1.nii.gz -force
mrcalc dti_V1.nii.gz -abs fac.nii.gz
module load afni
fat_proc_decmap -in_fa dti_FA.nii.gz -in_v1 fac.nii.gz -prefix DEC \
    -mask sub-1233_run-1_space-T1w_desc-brain_mask.nii.gz
# then look at the pngs in the QC folder
```

OK, now we need to come up with a game plan for all DTI we have.

```bash
# ncrshell
net_dir=/mnt/NCR/sudregp/MR_data_by_maskid;
rm -rf ~/tmp/dti_count.csv;
for m in `cat ~/tmp/ncr.txt`; do
    n8=`grep cdiflist08 $net_dir/${m}/2*/*README* | wc -l`;
    n9=`grep cdiflist09 $net_dir/${m}/2*/*README* | wc -l`;
    n99=`grep cdiflist99 $net_dir/${m}/2*/*README* | wc -l`;
    ng=`grep edti_g0 $net_dir/${m}/2*/*README* | wc -l`;
    echo $m,$n8,$n9,$n99,$ng >> ~/tmp/dti_count.csv;
done
```

```bash
# ncrshell01
net_dir=/mnt/NCR/sudregp/MR_data_by_maskid;
for m in `cat ~/tmp/batch2.txt`; do
    echo $m;
    cd /mnt/shaw/sudregp/NCR_BIDS/sub-${m};

    # find name of date folders
    ls -1 $net_dir/${m}/ | grep -e ^20 > ~/tmp/${m}_date_dirs;

    # for each date folder, check for dti scans
    mkdir -p dwi;
    while read d; do
        if grep -q cdiflist08 $net_dir/${m}/${d}/*README*; then
            nruns=3;
            tail -n +2 ~/research_code/cdiflist08 | split \
                --numeric-suffixes=1 -l 20;
        elif grep -q cdiflist09 $net_dir/${m}/${d}/*README*; then
            nruns=4;
            tail -n +2 ~/research_code/cdiflist09 | split \
                --numeric-suffixes=1 -l 20;
        else
            echo "No cdi sequence for $m";
        fi;
        # we need to convert the runs in order
        nii='-innii'
        bval='-inbval'
        bvec='-inbvec'
        for s in `seq 1 $nruns`; do
            grep cdiflist $net_dir/${m}/${d}/*README* | \
                grep -i _g0${s} > ~/tmp/${m}_scan;
            awk '{for(i=1;i<=NF;i++) {if ($i ~ /Series/) print $i}}' \
                ~/tmp/${m}_scan | sed "s/Series://g" > ~/tmp/${m}_scan_clean;
            mr_dir=`cat ~/tmp/${m}_scan_clean | sed "s/,//g"`;
            echo "Converting ${net_dir}/${m}/${d}/${mr_dir}/";
            dcm2niix_afni -o dwi/ -z y -f sub-${m}_run-${s}_dwi \
                ${net_dir}/${m}/${d}/${mr_dir}/;
            1dDW_Grad_o_Mat++           \
                -in_col_vec        x0${s}          \
                -out_row_vec       x0${s}_rowvec
            # replace PhaseEncodingAxis by PhaseEncodingDirection in JSON to
            # conform with BIDS
            sed -i -e "s/PhaseEncodingAxis/PhaseEncodingDirection/" \
                dwi/sub-${m}_run-${s}_dwi.json;
            # construct string
            nii=$nii' 'sub-${m}_run-${s}_dwi.nii.gz;
            bval=$bval' 'sub-${m}_run-${s}_dwi.bval;
            bvec=$bvec' '../x0${s}_rowvec;
        done
        cd dwi;
        fat_proc_convert_dcm_dwis $nii $bval $bvec -flip_z -no_qc_view \
            -prefix sub-${m}_dwi;
        # cleaning up
        mv sub-${m}_run-1_dwi.json sub-${m}_dwi.json;
        mv sub-${m}_dwi_rvec.dat sub-${m}_dwi.bvec;
        mv sub-${m}_dwi_bval.dat sub-${m}_dwi.bval;
        rm *run-* *dat *txt;
        cd ..;
        rm -f x??;
    done < ~/tmp/${m}_date_dirs;
done;
```


# TODO
 * use checkFlip script on PNC dataset!
 * check flip for other sequencs too
 * test that directions make sense (/scratch/sudregp/ncr_qsiprep)
 * add 99
 * add new sequence check that each sequence has the correct number of volumes
 * ignore sequences that don't have everything



