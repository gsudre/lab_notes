# 2021-06-16 10:35:53

Here's how we're running MRIQC for ABCD. First, get a list of how many have 4 in
the BIDS directory:

```bash
rm -rf rest_count.txt;
for s in `cat subjs.txt`; do
    echo $s;
    nrests=`/bin/ls -1 /data/ABCD_MBDU/abcd_bids/bids/${s}/ses-baselineYear1Arm1/func/*rest*gz | wc -l`;
    echo $s,$nrests >> rest_count.txt;
done;
```

I'll start with a batch of 11 of the people with 4 runs.

```bash
# bw
cd /scratch/sudregp

main_bids=/data/ABCD_MBDU/abcd_bids/bids/;
out_mriqc=/scratch/sudregp/mriqc_output;
bname=batch1;

rm -rf swarm.$bname; 
for m in `cat batch1.txt`; do
    echo "bash ~/research_code/mriqc_wrapper.sh $m $main_bids $out_mriqc " >> swarm.$bname;
done
swarm -g 10 -t 4 --logdir trash_mriqc --gres=lscratch:10 --time 45:00 \
    -f swarm.$bname --partition quick,norm --job-name $bname -m mriqc/0.16.1
```

I might need to make this take a bit longer. For now, let's see which subjects
failed:

```bash
cd /scratch/sudregp/mriqc_output
rm -rf ../redo.txt
for m in `cat ../batch1.txt`; do
    nres=`ls sub-${m}/ses-*/anat/*json 2>/dev/null | wc -l`;
    if [[ $nres != 1 ]]; then
        echo $m >> ../redo.txt;
    fi;
done
```

```bash
# bw
cd /scratch/sudregp

main_bids=/data/ABCD_MBDU/abcd_bids/bids/;
out_mriqc=/scratch/sudregp/mriqc_output;
bname=redo;

rm -rf swarm.$bname; 
for m in `cat ${bname}.txt`; do
    echo "bash ~/research_code/mriqc_wrapper.sh $m $main_bids $out_mriqc " >> swarm.$bname;
done
swarm -g 40 -t 4 --logdir trash_mriqc --gres=lscratch:20 --time 4:00:00 \
    -f swarm.$bname --partition quick,norm --job-name $bname -m mriqc/0.16.1
```

```bash
rm -rf /scratch/sudregp/redo.txt
for m in `cat /scratch/sudregp/batch?.txt`; do
    nres=`ls /scratch/sudregp/mriqc_output/sub-${m}/ses-*/anat/*json 2>/dev/null | wc -l`;
    if [[ $nres != 1 ]]; then
        echo $m >> /scratch/sudregp/redo.txt;
    fi;
done
```

Time to see why the redos are not working. At this point I've run them 2 times,
with all the memory they should need. Are they not in BIDS?















While I'm re-running the redos, let's make sure the data collection is still
working:








```bash
#interactive
docker run -it --rm -v /Volumes/Shaw/NCR_BIDS/:/data:ro \
    -v /Volumes/Shaw/mriqc_output/:/out poldracklab/mriqc:latest /data /out \
    group --no-sub


# directories for mriqc to play with
mkdir ${TMPDIR}/mriqc_output
mkdir ${TMPDIR}/mriqc_work

mriqc ${TMPDIR}/BIDS/ ${TMPDIR}/mriqc_output participant \
    --participant_label ${s} -m T1w --no-sub --n_procs $SLURM_CPUS_PER_TASK \
    --mem_gb $SLURM_MEM_PER_NODE -w ${TMPDIR}/mriqc_work;



**Need to resubmit anything that didn't finish**