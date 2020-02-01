# 2020-01-31 17:50:00

Now, how do we choose which of the univariate results to report? Let's try an
overall FDR using all variables. I'll start with HI, which seems to have a
stronger effect:

```
> ps = read.csv('~/data/baseline_prediction/prs_start/all_ps.csv')
> hi_p = ps[ps$sx=='hi',]
> p2 = p.adjust(hi_p[,'pval'], method='fdr')
> hi_p[p2<.05,'var']
 [1] ATR_fa           IFO_fa           ADHD_PRS0.001000 ADHD_PRS0.000500
 [5] FSIQ             externalizing    VMI.beery        PS.wj           
 [9] DS.wj            VM.wj           
46 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> hi_p[p2<.1,'var']
 [1] ATR_fa           CST_fa           IFO_fa           ADHD_PRS0.001000
 [5] ADHD_PRS0.000500 FSIQ             externalizing    base_age        
 [9] VMI.beery        PS.wj            DS.wj            VM.wj           
46 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p = ps[ps$sx=='inatt',]
> p2 = p.adjust(inatt_p[,'pval'], method='fdr')
> inatt_p[p2<.05,'var']
[1] ADHD_PRS0.000500 ADHD_PRS0.001000 FSIQ             base_age        
[5] externalizing    VMI.beery        PS.wj            VM.wj           
46 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p[p2<.1,'var']
 [1] IFO_fa           ADHD_PRS0.000500 ADHD_PRS0.001000 ADHD_PRS0.000100
 [5] FSIQ             base_age         externalizing    VMI.beery       
 [9] DSF.wisc         PS.wj            DS.wj            VM.wj           
46 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
```

These are not too bad, but it could be better. Maybe if we select a different
thickness summary? We could also remove the PRS to see how FDR behaves. I'm also
thinking about removing base_age, for the same reason I didn't run the models
for sex. It's more of a covariate, not necessarily a test. And it'l be in the
final model anyways.

```
> ps = read.csv('~/data/baseline_prediction/prs_start/all_ps2.csv')
> hi_p = ps[ps$sx=='hi',]
> p2 = p.adjust(hi_p[,'pval'], method='fdr')
> hi_p[p2<.05,'var']
 [1] VMI.beery        PS.wj            DS.wj            VM.wj           
 [5] ATR_fa           IFO_fa           ADHD_PRS0.001000 ADHD_PRS0.000500
 [9] FSIQ             externalizing    OFC             
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> hi_p[p2<.1,'var']
 [1] VMI.beery        PS.wj            DS.wj            VM.wj           
 [5] ATR_fa           CST_fa           IFO_fa           ADHD_PRS0.001000
 [9] ADHD_PRS0.000500 base_age         FSIQ             externalizing   
[13] OFC              cingulate       
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p = ps[ps$sx=='inatt',]
> p2 = p.adjust(inatt_p[,'pval'], method='fdr')
> inatt_p[p2<.05,'var']
[1] VMI.beery        PS.wj            VM.wj            ADHD_PRS0.001000
[5] ADHD_PRS0.000500 base_age         FSIQ             externalizing   
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p[p2<.1,'var']
 [1] VMI.beery        DSF.wisc         PS.wj            DS.wj           
 [5] VM.wj            IFO_fa           ADHD_PRS0.000100 ADHD_PRS0.001000
 [9] ADHD_PRS0.000500 base_age         FSIQ             externalizing   
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
```

Yeah, I think this is a bit better because now we can include a few of the
thickness results, if especially if we go down to q < .1.

Let's try it without base_age:

```
> ps = ps[ps$var!='base_age',]
> hi_p = ps[ps$sx=='hi',]
> p2 = p.adjust(hi_p[,'pval'], method='fdr')
> hi_p[p2<.05,'var']
 [1] VMI.beery        PS.wj            DS.wj            VM.wj           
 [5] ATR_fa           IFO_fa           ADHD_PRS0.001000 ADHD_PRS0.000500
 [9] FSIQ             externalizing    OFC             
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> hi_p[p2<.1,'var']
 [1] VMI.beery        PS.wj            DS.wj            VM.wj           
 [5] ATR_fa           CST_fa           IFO_fa           ADHD_PRS0.001000
 [9] ADHD_PRS0.000500 FSIQ             externalizing    OFC             
[13] cingulate       
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p = ps[ps$sx=='inatt',]
> p2 = p.adjust(inatt_p[,'pval'], method='fdr')
> inatt_p[p2<.05,'var']
[1] VMI.beery        PS.wj            VM.wj            ADHD_PRS0.001000
[5] ADHD_PRS0.000500 FSIQ             externalizing   
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p[p2<.1,'var']
 [1] VMI.beery        DSF.wisc         PS.wj            DS.wj           
 [5] VM.wj            IFO_fa           ADHD_PRS0.000100 ADHD_PRS0.001000
 [9] ADHD_PRS0.000500 FSIQ             externalizing   
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
```

How about removing PRS?

```
> ps = ps[!grepl(ps$var, pattern='^ADHD'),]
> hi_p = ps[ps$sx=='hi',]
> p2 = p.adjust(hi_p[,'pval'], method='fdr')
> hi_p[p2<.05,'var']
[1] VMI.beery     PS.wj         DS.wj         VM.wj         ATR_fa       
[6] IFO_fa        FSIQ          externalizing OFC          
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> hi_p[p2<.1,'var']
 [1] VMI.beery     DSF.wisc      PS.wj         DS.wj         VM.wj        
 [6] ATR_fa        CST_fa        IFO_fa        FSIQ          externalizing
[11] OFC           cingulate    
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p = ps[ps$sx=='inatt',]
> p2 = p.adjust(inatt_p[,'pval'], method='fdr')
> inatt_p[p2<.05,'var']
[1] VMI.beery     PS.wj         VM.wj         FSIQ          externalizing
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
> inatt_p[p2<.1,'var']
[1] VMI.beery     DSF.wisc      PS.wj         DS.wj         VM.wj        
[6] IFO_fa        FSIQ          externalizing
43 Levels: ADHD_PRS0.000050 ADHD_PRS0.000100 ... VMI.beery
```

Now it's just a matter of selecting which one we want to go with to put in the
big model...

# TODO
* construct model with everything
* ML
* are we still using rsfmri?