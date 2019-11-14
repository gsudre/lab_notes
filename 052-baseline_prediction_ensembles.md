# 2019-11-12 11:47:33

Let's give it a try with ensembles. First, RandomForest. It takes a while to
run, so I'll reduce it to 30min per swarm.

```bash
# bw
export OMP_NUM_THREADS=32
cd ~/data/baseline_prediction/manual_swarms

jname=rfImp;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/rf_classifier.py;
res=~/data/baseline_prediction/manual_results;
for s in `cat ../random25.txt`; do
    for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
        phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $vars $res $s" >> $swarm_file;
            done;
        done;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 8 --logdir=trash_${jname} \
    --job-name ${jname} --time=30:00 --merge-output --partition quick,norm
```

And for XGB we need a different OMP. But I haven't played much with the
configuration parameters, so it's not taking very long.

```bash
# bw
export OMP_NUM_THREADS=1
cd ~/data/baseline_prediction/manual_swarms

jname=xgbImp;
swarm_file=swarm.${jname};
rm -f $swarm_file;
code=~/research_code/baseline_prediction/xgb_classifier.py;
res=~/data/baseline_prediction/manual_results;
for s in `cat ../random25.txt`; do
    for p in dti_fa dti_ad dti_rd struct_area struct_volume struct_thickness; do
        phen=~/data/baseline_prediction/${p}_OD0.95_11052019.csv;
        for i in Next Last Study; do
            for j in SX_inatt SX_HI; do
                echo "python3 $code $phen ${j}_group${i} $vars $res $s" >> $swarm_file;
            done;
        done;
    done;
done;
swarm --gres=lscratch:10 -f $swarm_file -t 32 -g 20 -b 48 --logdir=trash_${jname} \
    --job-name ${jname} --time=5:00 --merge-output --partition quick,norm
```

# 2019-11-13 14:25:42

Let's plot the RF and XGB results to see if there is anything earth-shattering.

```r
mydir='~/data/baseline_prediction/manual_results/'
par(mfrow=c(2, 3))
for (j in c('SX_inatt', 'SX_HI')) {
    for (i in c('Next', 'Last', 'Study')) {
        target = sprintf('%s_group%s', j, i)
        tmp = c()
        for (phen in c('dti_fa', 'dti_ad', 'dti_rd', 'struct_area',
                        'struct_volume', 'struct_thickness')) {
            for (pipe in c('_XGB', '_RF')) {
                data = read.csv(sprintf('%s/classification_results%s_%s_OD0.95_11052019.csv',
                                mydir, pipe, phen), header=0)
                res_rows = which(grepl(data$V1, pattern=target))
                tmp = rbind(tmp,
                            data.frame(group=sprintf('%s_%s (%d)',
                                                     phen, pipe, length(res_rows)),
                                       val=data[res_rows, 'V3']))
            }
        }
        mytitle = sprintf('%s', target)
        boxplot(as.formula('val ~ group'), data=tmp, main=mytitle, ylim=c(0.2,.9),
                ylab='ROC', las=2, xlab='')
    }
}
```


How about using this for confidence interval of validation set? 

https://stackoverflow.com/questions/19124239/scikit-learn-roc-curve-with-confidence-intervals

If we use Uher's paper as an acceptable way to do this
(https://www.nature.com/articles/s41598-018-23584-z),
they only showed a confidence interval on the validation set. And they were
explicit that they only did the separation once... maybe we should do that? Then
the only question is how to pick the seed.

Let me check if there is a way to pick the seed that does best across
datasets... in a way, it would have to be a somewhat meaningful seed. Say, 0, 42,
1234, or 2019.


# TODO
 * Try balancing the data?
 * Do further tuning of XGB
 * Try adding in rsfmri to see if things get better
 * Try combining the datasets
 * Final effort will be to do a descriptive paper, and then have a second part
   showing that similar methods don't work in a predictive framework 
