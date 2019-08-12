# 2019-08-12 16:44:41

It's possible that my MELODIC results are not so good because my mask is taking
into considerationt oo much stuff that's not grey matter. So, let's repeat the
results using a better grey matter mask, and see how they look. I'll also redo
the yeo results just for kicks, but stick to FD1.0 for now.

And, since we're at it, why not compute everything for Z and noZ?

```bash
cd ~/data/heritability_change/xcp-36p_despike
melodic -i fd1_epi.txt -o groupmelodic_gray.ica -v --nobet -m gray_matter_mask.nii --tr=2.5 --report --Oall -a concat;
```


