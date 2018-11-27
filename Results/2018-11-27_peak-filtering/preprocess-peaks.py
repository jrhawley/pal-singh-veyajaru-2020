"""
preprocess-peaks
==========

1. Remove chrM, and non-canonical chromosomes
2. Sort the peaks
3. Remove summits, only keeping unique peaks and retaining the largest -log10(q)
among a set of identical peaks
4. Count number of peaks per threshold and generate plots
"""

from __future__ import division, absolute_import, print_function
import os.path
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import gzip
import pybedtools as pbt

# ==============================================================================
# Constants
# ==============================================================================
LOG10Q_THRESHOLDS = np.linspace(1, 5, 5)
peak_metadata = pd.read_table("config.tsv", index_col=False)
peak_metadata["Filtered"] = [
    os.path.join("BedGraphs", c + "_Rep" + str(r) + ".filtered.bedGraph")
    for (c, r) in zip(peak_metadata["Condition"], peak_metadata["Replicate"])
]
peak_metadata["Sorted"] = [
    os.path.join("BedGraphs", c + "_Rep" + str(r) +
                 ".filtered.sorted.bedGraph")
    for (c, r) in zip(peak_metadata["Condition"], peak_metadata["Replicate"])
]
peak_metadata["Unique"] = [
    os.path.join("BedGraphs", c + "_Rep" + str(r) +
                 ".filtered.sorted.unique.bedGraph")
    for (c, r) in zip(peak_metadata["Condition"], peak_metadata["Replicate"])
]
CONDITIONS = np.unique(peak_metadata["Condition"])
REPLICATES = np.unique(peak_metadata["Replicate"])
# ==============================================================================
# Functions
# ==============================================================================


def filter_chr(narrowPeak, outfile):
    """
    Filter mitochondrial and non-canonical chromosomes from a narrowPeak file

    Parameters
    ----------
    narrowPeak : str
        Path to narrowPeak file
    outfile : str
        Path to output file
    """
    chrs = ["chr" + str(c) for c in range(1, 23)] + ["chrX", "chrY"]
    f_in = gzip.open(narrowPeak, "rt")
    f_out = open(outfile, "w")
    for line in f_in:
        splitline = line.rstrip().split("\t")
        chrom = splitline[0]
        start = splitline[1]
        end = splitline[2]
        logq = splitline[8]
        if chrom in chrs:
            f_out.write("\t".join([chrom, start, end, logq]) + "\n")
    f_in.close()
    f_out.close()


def remove_dups(bg, outfile):
    """
    Remove duplicate peaks due to --call-summits, and keep the largest -log10(q)

    Parameters
    ----------
    bg : str
        Path to input peak bedGraph file
    outfile : str
        Path to output file
    """
    f_in = open(bg, "r")
    f_out = open(outfile, "w")
    # read first line
    prev = f_in.readline().rstrip().split("\t")
    for line in f_in:
        splitline = line.rstrip().split("\t")
        # check if chr, start, and end are the same as the previous line
        # (can do this because input is sorted)
        if splitline[0:3] == prev[0:3]:
            if splitline[3] > prev[3]:
                # replace -log10(q) if this peak has a larger value
                prev[3] = splitline[3]
        else:
            # print prev to outfile if this is a new locus
            f_out.write("\t".join(prev) + "\n")
            # replace prev for next comparison
            prev = splitline
    f_in.close()
    f_out.close()


def file_len(file):
    """
    Count the number of lines in the file

    Parameters
    ----------
    file : str
        Path to file
    """
    return(sum(1 for line in open(file)))


def filter_q_thresh(bg, name):
    """
    Filter bedGraph created from narrowPeak by different q-value thresholds

    Parameters
    ----------
    bg : str
        Path to input peak bedGraph file
    name : str
        Basename for output files
    """
    peaks = pd.read_table(bg, index_col=False, header=None, names=[
                          "chr", "start", "end", "logq"])
    # check that "Filter" folder exists. if not, make it
    if not os.path.exists("Filter"):
        os.mkdir("Filter")
    for t in LOG10Q_THRESHOLDS:
        thresh_name = "logq_" + str(t)
        # check that "logq_t" folder exists. if not, make it
        if not os.path.exists(os.path.join("Filter", thresh_name)):
            os.mkdir(os.path.join("Filter", thresh_name))
        peaks[peaks["logq"] >= t].to_csv(
            os.path.join("Filter", thresh_name, name))

# ==============================================================================
# Main
# ==============================================================================
print("Filtering chromosomes")
for i in range(len(peak_metadata)):
    print("\t", peak_metadata["Condition"].iloc[i],
          peak_metadata["Replicate"].iloc[i])
    filter_chr(
        peak_metadata["narrowPeak"].iloc[i],
        peak_metadata["Filtered"].iloc[i]
    )

print("Sorting BedGraphs")
for i in range(len(peak_metadata)):
    print("\t", peak_metadata["Condition"].iloc[i],
          peak_metadata["Replicate"].iloc[i])
    bed = pbt.BedTool(peak_metadata["Filtered"].iloc[i])
    bed_sorted = bed.sort()
    bed_sorted.saveas(peak_metadata["Sorted"].iloc[i])

print("De-duplicating summits and peaks")
for i in range(len(peak_metadata)):
    print("\t", peak_metadata["Condition"].iloc[i],
          peak_metadata["Replicate"].iloc[i])
    remove_dups(
        peak_metadata["Sorted"].iloc[i],
        peak_metadata["Unique"].iloc[i]
    )

print("Filtering peaks at multiple thresholds")
# count number of peaks in each filtered file
counts = pd.concat([peak_metadata[["Condition", "Replicate"]]
                    for i in range(len(LOG10Q_THRESHOLDS))],
                   ignore_index=True)
counts["Threshold"] = pd.Series(
    np.repeat(LOG10Q_THRESHOLDS, len(peak_metadata)), index=counts.index)
counts["Peaks"] = 0

for i in range(len(peak_metadata)):
    cond = peak_metadata["Condition"].iloc[i]
    rep = peak_metadata["Replicate"].iloc[i]
    print("\t", cond, rep)
    base = cond + "_Rep" + str(rep) + ".bedGraph"
    # filter_q_thresh(peak_metadata["Unique"].iloc[i], base)
    for t in LOG10Q_THRESHOLDS:
        filt_file = os.path.join("Filter", "logq_" + str(t), base)
        counts.loc[(counts["Condition"] == cond) & (counts["Replicate"] == rep) & (
            counts["Threshold"] == t), "Peaks"] = file_len(filt_file)

# print counts to output file
counts.to_csv(os.path.join("Filter", "peak-counts.tsv"), sep="\t")

print("Plotting peak counts per threshold")
plt.figure(figsize=(12, 12))
g = sns.FacetGrid(counts, row="Condition", col="Replicate", margin_titles=True)
g.map(plt.bar, "Threshold", "Peaks", tick_label=LOG10Q_THRESHOLDS)
g.fig.subplots_adjust(right=.95)
plt.savefig(os.path.join("Filter", "peak-counts.png"))
plt.close()

print("Done")
