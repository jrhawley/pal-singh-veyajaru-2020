# 2018-11-29

I've tried to use DEseq2 to find differential cloci between the MB6-treated and Ctrl samples.
At the end, I'm getting a strange, curved histogram for the distribution of p-values, instead of a mostly flat histogram.

I'm not completely clear on where to go from here within DEseq2, so I'm going to try using EdgeR as another analysis tool

# 2018-11-30

EdgeR has been installed, and I'm trying it out now.
I'm getting a much more "standard" set of results, it looks like, with 5.5% of sites being significantlt differentially accessible between the 2 conditions.
I'm going to ask Zhaleh if she has any thoughts on why these things are so different from each other.

# 2018-12-07

I'm going to try using DiffBind.
It's a Bioconductor package developed for differential TF binding detection from ChIP-seq data, but ATAC-seq should work analogously, and has been used by others in previous works.

Using DiffBind with the DEseq2 engine, the histogram of p-values is much more reasonable than when I generated it myself, both for the raw and filtered peaks.
There may be an error in my own implementation, so I think I should stick to DiffBind, instead of doing it myself.

The results between the raw and filtered cases are very different, surprisingly.
There are > 2000 significantly differentially accessible sites in the filtered peaks and ~ 1400 in the raw peaks.
I suspect that the decrease in sites is due to decreased power by performing for tests (more peaks means more tests, means more stringent filtering).

# 2018-12-10

I've switched the annotation so that group 1 is MB6, making all downstream results of the form "MB6 - Ctrl", which makes more intuitive sense, at a glance.
With the 2005 significantly differentially accessible loci between the MB6 and Ctrl samples, I'm going to annotate them with GENCODE v19 (latest version for hg19, which is what the FASTQs were aligned to).

# 2018-12-11

Using the edgeR-based analysis with DiffBind on the filtered peaks, we find 2003 significantly differentially accessible loci between the two conditions.
1572 of these have increased accessibility in the MB6 samples, and the remaining 431 of them have decreased accessibility.
Globally, the MB6 samples have less accessible chromatin than the Ctrl samples (see `../2018-11-28_global-accessbility/`).

I tried using GREAT to see if there was any quick ontology interpretation that I could make.
Using a q-value threshold of 0.05 and a log2FC threshold of 1.5, there are 193 significantly differentially accessible sites, and no terms were returned by GREAT.

# Conclusions

I haven't mapped any of the significantly differentially accessible sites to genes or proteins from Rashim's model, so the most conclusive this that I can say is that there are 193 significantly differentially accessible sites, 161 decreased, and 32 increased, in the MB6 samples compared to the Ctrls.
Additionally, there is a global reduction in accessibility.
I've sent this information along to Rashim, and hopefully he is able to interpret these results.
