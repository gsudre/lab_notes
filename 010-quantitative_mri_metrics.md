# 2019-03-22 17:21:40

I wanted to play a little bit with some quantitative MRI metrics, both
anatomical and functional. There's always MRIQC, but I'm not a big fan of how it
needs Docker and many different tools. So, I found this:

http://preprocessed-connectomes-project.org/quality-assessment-protocol/

which is actually the backbones metrics of MRIQC. Also, it only requires AFNI,
which is a plus.

Then, even though I have already done much of the grunt work in collecting the
Freesurfer metrics, there is a nice script from here that does everything too,
including collection and QC metrics:

https://github.com/PennBBL/conte/wiki/FreeSurfer-v5.3-QA-on-CFN

So, I just need to tweak it a bit to do it for our data.

Finally, I found this tool that actually predicts visual QCs based on a trained
supervised model, or we could even train it with our own raters (to the extent
we trust them):

https://github.com/Qoala-T/QC

So, that's an option as well. I'd prefer to use the multi-dimensional outlier
detection approach fist, but they're certainly not exclusive.

## QAP

OK, since I'm about to dive into rsFMRI again, let's run QAP first so we can get
some extra QC metrics for functional MRI.

I'll do it in ncrshell01 because its access to the shared drive is much faster.

```bash
# Had to install using python 2.7, because QP didn't seem compatible with 
# python 3! The install kept breaking when trying to install using pip
source activate QAP
```

# 2019-03-25 13:36:29

Nevermind, even after successfully installing it, I saw several issues with
other dependencies, that the old code didn't work. Well, I'll have to use mriqc
then. Let's put the resting state data into BIDS format, then. Maybe just the
processed data, right? Makes no sense to convert files that AFNI failed...

But for that I can just go back to the 003 note.
