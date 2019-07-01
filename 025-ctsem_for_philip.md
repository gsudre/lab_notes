# 2019-07-01 13:05:10

I created a script for Philip to run the CTSEM analysis in
~/research_code/ctsem_loop.R. The cript runs alright with a for loop for SX and
different phenotypes, but then I recoded it so it would take a few variables as
arguments, and this way we could run it in parallel in the cluster.

To do it that way, I'm going to set it up in an interactive node, because it'll
make it easier for other people to run it (instead of swarming). But of course,
it could be done in a swarm fashion later.

```bash
#bw in
TMPDIR=/lscratch/$SLURM_JOBID;
cd $TMPDIR;
params=./params.txt;
cp ~/tmp/wider_devel_names.csv $TMPDIR/
for sx in SX_inatt SX_HI SX_total; do
    for p in fa_lcst ixi_ADC_left_cst; do
        echo $sx >> $params;
        echo $p >> $params;
        echo $TMPDIR/res_TI1_${sx}_${p}.csv >> $params;
    done;
    for p in {5..46}; do
        echo $sx >> $params;
        echo Y${p} >> $params;
        echo $TMPDIR/res_TI1_${sx}_Y${p}.csv >> $params;
    done;
done;

cat $params | parallel --max-args=3 -j $SLURM_CPUS_PER_TASK \
    Rscript ~/research_code/ctsem_loop.R $TMPDIR/wider_devel_names.csv {1} {2} {3};
```

Then, we need to compile all resultas into a single file:

```bash
echo '"lbound","estimate","ubound","StdError","var","sx","phen"' > out_TI1.csv;
for f in `ls $TMPDIR/res_TI1_*csv`; do
    tail -n +2 $f >> out_TI1.csv;
done
```

That's taking longer than what I expected. Let's swarm it and hope it doesn't
cause much problem in the cluster:

```bash
#bw
TMPDIR=~/data/philip;
cd $TMPDIR;
params=swarm.ctsem;
cp ~/tmp/wider_devel_names.csv $TMPDIR/
for sx in SX_inatt SX_HI SX_total; do
    for p in fa_lcst ixi_ADC_left_cst; do
        echo Rscript ~/research_code/ctsem_loop.R $TMPDIR/wider_devel_names.csv $sx $p $TMPDIR/res_TI1_${sx}_${p}.csv >> $params;
    done;
    for p in {5..46}; do
        echo Rscript ~/research_code/ctsem_loop.R $TMPDIR/wider_devel_names.csv $sx Y${p} $TMPDIR/res_TI1_${sx}_Y${p}.csv >> $params;
    done;
done;
swarm -f ${params} --module R --logdir=trash_ctsem --job-name ctsem \
    --partition quick --time=20:00;
```

