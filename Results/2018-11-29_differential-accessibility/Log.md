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

