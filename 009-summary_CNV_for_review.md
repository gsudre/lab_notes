# 2019-03-19 16:27:38

It's been a while since I went through these results, but I'm using the Slack
data to refresh my memory on what we did. From the get-go, I just noticed that
GATK 4.0 released their own tool to do CNV discovery in WES germline data, which
was not available at the time I was doing this analysis. It's still somewhat
experimental in the sense that there is no official documentation for it, just a
few notes here and there. Still, it might be interesting to try it out, as much
of our pipeline uses GATK anyways, and possibly compare to the results in the
other tools we're exploring.

So, the crux of it is that there really isn't a standard method to do this
analysis out there. So, we were trying to evaluate the outputs of different
tools to get a comprehensive view of the land.

## Background

Based on my previous notes, CNV methods use the whole exome sequence, aligned to
a pre-determined reference genome, to make inferences about the copy number of
the different parts of the exome.

It made sense to start with the GATK Best Practices processing pipeline because
of its wide-acceptance in the field, and we were already using it for other
denovo detection methods. CNV tools pick up the WES outputs after alignment and
marking of duplicate reads. 

One common limitation in selecting tools for CNV analysis was the ability to
find denovo variants (DNVs) without the need of a case/control structure. Much
of the work in the field is done for cancer research and depends on that design,
but it’s not our case. Another common requirement was regular maintenance of the
software, and favorable mention in survey papers comparing different tools. 

The tool also needed to be shown to work on WES, as many tools were developed
for WGS and their assumptions and methods do not apply for WES. Also, I chose
tools that have specific guidelines in their papers on how to analyze the found
CNVs to identify DNVs. Finally, even though all tools in the CNV arm use read
depth to infer copy number, the 3 tools chosen have different and complementary
(in my opinion) methods to interpret that data.

For CNV detection we are using a combination of software tools, focusing on consensus calls, given the lack of agreement in the results of many tools that look for CNVs in exome sequencing data. XHMM is a fairly popular tool for CNV detection, and uses a combination of PCA and a Hidden Markov Model to detect rare CNVs based on a batched-comparison principle. A different approach is taken in CNVkit, which uses read depth information in targeted reads and the nonspecifically captured off-target reads to infer CNVs. Survey research papers have shown that the calls by these CNV methods can be complementary and, as in the case for the SNV arm, we might choose to drop or even add different tools as we start to apply them to our recently acquired data (e.g. Excavator2Tool, ExomeCopy, ExomeDepth)

Based on my notes, I tried but decided to later drop CNVnator due to its poor
implementation and lack of clarity of performance on WES. 

Lastly, the goal was to combine the calls across tools to hopefully reduce the
number of false positive calls, which seems to be an issue in the literature.

The two papers I was using to construct my pipeline were:

https://molecularcytogenetics.biomedcentral.com/articles/10.1186/s13039-017-0333-5

and

http://onlinelibrary.wiley.com/doi/10.1002/humu.22537/full

There is also this one:

https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2374-2

and this critiscism by Les:

https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-016-0336-6


## Current results
We have 99 WES sequences, representing 21 different nuclear families. They are split into 21 affected and 29 unaffected trios (some families have more than one unaffected sibling).

The gist of the analysis is finding DNVs in each trio structure. Then, within families, mark which DNVs are present in affected but not in the unaffected siblings. Finally, we compute statistics based on how many of those DNVs appeared across affected trios, and find biological plausibility for the DNVs that pop up.

I also tried PennCNV to look at CNV using the bead arrays, but I found no
results there.

For CNVs, we have called CNVS in all samples in both XHMM and CNVkit. We're now
assessing which variants are rare and de novo in all trios, and possibly test
other WES CNV tools if there is no consensus in the results.