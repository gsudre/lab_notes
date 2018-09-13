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

