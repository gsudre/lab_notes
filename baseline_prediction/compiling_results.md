# 2018-09-12 14:46:38

Here's how we currently grab the results:

```bash
echo "target,pheno,var,nfeat,model,metric,val" > auto_summary.csv;
for y in diag_group2 random_HI_slope random_total_slope random_inatt_slope \
    OLS_inatt_slope OLS_HI_slope OLS_total_slope \
    group_HI3 group_total3 group_INATT3; do
    for f in `grep -l \"${y} *o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 4 | cut -d"_" -f 1,2 -`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        nfeat=`grep -e "26. " $f | tail -1 | awk '{ print $3 }'`;
        echo $y,$phen,$var,$nfeat,$model,$metric,$acc >> auto_summary.csv;
    done;
done
```

It looks like uni01 is not doing much. For the ones that it actually run (sometimes there were no variables p < .01)), the results were not better than when using uni. That could indicate overfitting, but also there could be features that interact well and hacving so few features (As in uni01) means shooting ourselves in the foot.

# 2018-09-14 12:34:50

I updated the code to use all the improvements/enhancements we have in the limited version, except that it's not limited to algorithms or time. In fact, I think the filtering to the structural data puts the number of variables in a manageable number, so we don't need to run the limited code anymore. And I have alreayd noticed that the raw code breaks because it runs out of memory... let's run it anyways, and let it die.

```bash
rm swarm.automl
for var in raw pca uni pca_uni uni_pca uni01; do
    for m in area thickness volume; do
        for sx in inatt HI total; do
            for sl in OLS random; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/struct_${m}_08312018_260timeDiff12mo.RData.gz ~/data/baseline_prediction/struct_gf_09052018_260timeDiff12mo.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
            done;
        done;
        for sx in INATT HI total; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/struct_${m}_08312018_260timeDiff12mo.RData.gz ~/data/baseline_prediction/struct_gf_09052018_260timeDiff12mo.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
        done;
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/struct_${m}_08312018_260timeDiff12mo.RData.gz ~/data/baseline_prediction/struct_gf_09052018_260timeDiff12mo.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
    done;
done;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl;
swarm -f swarm.automl -g 40 -t 32 --time 1-00:00:00 --logdir trash_bin --job-name struc -m R --gres=lscratch:10
```

We can also run the snp code. I also don't think there's much value in running the PCA code there, so it's just the raw code.

```bash
rm swarm.automl
var=raw;
for m in {4..9}; do
    for sx in inatt HI total; do
        for sl in OLS random; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_snps1e0${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
        done;
    done;
    for sx in INATT HI total; do
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_snps1e0${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
    done;
    echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_snps1e0${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
done;
m='prs'
for sx in inatt HI total; do
    for sl in OLS random; do
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
    done;
done;
for sx in INATT HI total; do
    echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
done;
echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/geno3_${m}_09132018.RData.gz ~/data/baseline_prediction/geno3_gf_09142018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
sed -i -e "s/^/unset http_proxy; /g" swarm.automl;
swarm -f swarm.automl -g 40 -t 32 --time 1-00:00:00 --logdir trash_bin --job-name geno -m R --gres=lscratch:10
```

And since we're running all for comparison, let's attach neuropsych as well:

```bash
for var in raw pca uni pca_uni uni_pca uni01; do
    for m in all cpt wiscraw wiscstd wj iq; do
        for sx in inatt HI total; do
            for sl in OLS random; do
                echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/cog_${m}_09142018.RData.gz ~/data/baseline_prediction/cog_gf_09142018.csv ${sl}_${sx}_slope ~/tmp/${m}_${sl}_${sx}" >> swarm.automl;
            done;
        done;
        for sx in INATT HI total; do
            echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/cog_${m}_09142018.RData.gz ~/data/baseline_prediction/cog_gf_09142018.csv group_${sx}3 ~/tmp/${m}_${sx}" >> swarm.automl;
        done;
        echo "Rscript --vanilla ~/research_code/automl/${var}.R ~/data/baseline_prediction/cog_${m}_09142018.RData.gz ~/data/baseline_prediction/cog_gf_09142018.csv diag_group2 ~/tmp/${m}_diag_group2" >> swarm.automl;
        done;
    done;
done;
```

# 2018-09-18 10:03:19

Then we collect all results again, but with a small change to the script so we
can correctly capture the number of features used:

```bash
echo "target,pheno,var,nfeat,model,metric,val" > auto_summary.csv;
for y in diag_group2 random_HI_slope random_total_slope random_inatt_slope \
    OLS_inatt_slope OLS_HI_slope OLS_total_slope \
    group_HI3 group_total3 group_INATT3; do
    for f in `grep -l \"${y} *o`; do
        phen=`head -n 2 $f | tail -1 | awk '{FS=" "; print $6}' | cut -d"/" -f 4 | cut -d"_" -f 1,2 -`;
        var=`head -n 2 $f | tail -1 | awk '{FS=" "; print $5}' | cut -d"/" -f 4 | sed -e "s/\.R//g"`;
        model=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $2}' | cut -d"_" -f 1`;
        acc=`grep -A 1 model_id $f | tail -1 | awk '{FS=" "; print $3}'`;
        metric=`grep -A 0 model_id $f | awk '{FS=" "; print $2}'`;
        if [[ $var == 'raw' ]] && [[ $phen == 'dti_ALL' ]]; then
            nfeat=36066;
        elif [[ $var == 'raw' ]] && [[ $phen = *'dti_'* ]]; then
            nfeat=12022;
        elif grep -q "Running model on" $f; then  # newer versions of the script
            nfeat=`grep "Running model on" $f | awk '{FS=" "; print $5}'`;
        else # older versions were less verbose
            nfeat=`grep -e "26. *[1-9]" $f | grep "\[1\]" | tail -1 | awk '{ print $3 }'`;
        fi
        echo $y,$phen,$var,$nfeat,$model,$metric,$acc >> auto_summary.csv;
    done;
done
```

So I saved auto_summary_09182018.csv in ~/Dcouments/baseline_prediction.

