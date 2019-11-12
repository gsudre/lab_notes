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

# TODO
 * Try balancing the data?
 * Do further tuning of XGB
 * Try adding in rsfmri to see if things get better
 * Try combining the datasets
 * Final effort will be to do a descriptive paper, and then have a second part
   showing that similar methods don't work in a predictive framework 

I will be providing clerical and managerial services to companies in the
health-care domain. These activities include client relationship, deployment of
already-established products to assure patient care, and hiring staff to perform
IT-related tasks, such as website deployment and maintenance.

The company intends to seek federal government grants. In those situations, any
compensation for my services will be kept separate from such funds.

I am responsible for analyzing and integrating data from different domains, such
as genomics, brain imaging, and behavioral data, to understand
neurodevelopmental disorders of childhood. Currently, our lab focuses on
Attention Deficit and Hyperactivity Disorder.

