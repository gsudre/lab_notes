# 2018-10-01 17:18:30

I started playing with the list of binary SX (clinics_binary_sx_baseline.xlsx). Using the DICA and NV_interview
sheets (also, simplex versions), I couldn't find binary SX for 2 out of the 380
in their baseline date. Then, there was 3 with missing interview, only SX, so
they had to be removed.

I'm waiting to hear from Wendy on 7620214, and I had to correct SX for 4848779
because it was wrong (inverted) in Philip's sheet... DICA is correct based on
binary SX.

# 2018-10-02 06:25:58

I heard back from Wendy that the SX was wrong, so I fixed it in the binary file. It'll likely affect Philip's files and slopes, but we'll have to worry about that later. Same thing about the inverted ID too.

I then played with the file to have only Y/N and corrected the variable names (clinics_binary_sx_baseline_10022018.csv), and save an RData compressed file just to keep with past convention.

Now, it's just a matter of swarming it:

```bash
fname=swarm.automl_clinics
rm $fname;
function print_commands () {
    for nn in nonew_ ''; do
        for g in ADHDonly_ ''; do
            for t in diag_group2 OLS_inatt_slope OLS_HI_slope OLS_total_slope \
                random_HI_slope random_total_slope random_inatt_slope \
                group_HI3 group_total3 group_INATT3; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} /data/NCR_SBRB/baseline_prediction/models/${nn}${g}${t} 42" >> $fname;
            done; 
        done;
        # diag_group2 doesn't apply to ADHD_NOS
        for g in ADHDNOS_ ADHDNOS_group; do
            for t in OLS_inatt_slope OLS_HI_slope OLS_total_slope \
            random_HI_slope random_total_slope random_inatt_slope; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} /data/NCR_SBRB/baseline_prediction/models/${nn}${g}${t} 42" >> $fname;
            done; 
        done;
        g=ADHDNOS_;
        for t in group_HI3 group_total3 group_INATT3; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R $f /data/NCR_SBRB/baseline_prediction/long_clin_0918.csv ${nn}${g}${t} /data/NCR_SBRB/baseline_prediction/models/${nn}${g}${t} 42" >> $fname;
        done;
    done;
}

var=raw;
f=/data/NCR_SBRB/baseline_prediction/clinics_binary_sx_baseline_10022018.RData.gz;
print_commands

sed -i -e "s/^/unset http_proxy; /g" $fname;
wc -l $fname;
swarm -f $fname -g 40 -t 32 --time 4:00:00 --partition quick --logdir trash_clinics --job-name clinics -m R --gres=lscratch:10;
```


