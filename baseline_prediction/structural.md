# 2018-09-12 14:39:51

We're having some issues running even the limited version of the autoML scripts.
Not sure what's going on, but it doesn't even look like we're getting to run
automl, as we're not stopping after some time. Let's run a few tests...

# 2018-09-13 11:52:26

It looks like it takes a long time to just prepare the data. But using all the
variables in the structural (voxels) we already go up to 115Gb, and then the
models start running. I'm sure it will break as soon as they run longer... ye,
it did.

That's fine, as our best models were all using univariate, so we could go with
that, or even PCA. We just need to make sure the code is working fine for that.


