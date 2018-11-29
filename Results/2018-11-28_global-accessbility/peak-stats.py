"""
peak-stats.py
==========

Basic summary statistics on a set of narrowPeak files
"""

from __future__ import division, absolute_import, print_function
import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# ==============================================================================
# Constants
# ==============================================================================

# ==============================================================================
# Constants
# ==============================================================================
peak_metadata = pd.read_table("config.tsv", index_col=False)
CONDITIONS = np.unique(peak_metadata["Condition"])
REPLICATES = np.unique(peak_metadata["Replicate"])

# ==============================================================================
# Functions
# ==============================================================================


def file_len(file):
    """
    Count the number of lines in the file

    Parameters
    ----------
    file : str
        Path to file
    """
    return(sum(1 for line in open(file)))


# ==============================================================================
# Main
# ==============================================================================
# Counting peak sizes
# -------------------------------------
# aggregate data placeholder for peaks from each sample
print("Counting peaks")

# read first sample to initialize data frame
cond = peak_metadata.loc[0, "Condition"]
repl = peak_metadata.loc[0, "Replicate"]
peakfile = peak_metadata.loc[0, "narrowPeak"]
print("\t", cond, repl)
# count the number of peaks
peak_metadata.loc[0, "PeakCount"] = file_len(peakfile)

peaks = pd.read_table(
    peak_metadata.loc[0, "narrowPeak"],
    index_col=False,
    header=None,
    sep="\t",
    names=["chr", "start", "end"],
    dtype={"chr": str, "start": int, "end": int})
peaks["Condition"] = peak_metadata.loc[0, "Condition"]
peaks["Replicate"] = peak_metadata.loc[0, "Replicate"]

for i in range(1, len(peak_metadata)):
    cond = peak_metadata["Condition"].iloc[i]
    repl = peak_metadata["Replicate"].iloc[i]
    peakfile = peak_metadata["narrowPeak"].iloc[i]
    print("\t", cond, repl)
    # count the number of peaks
    peak_metadata.loc[i, "PeakCount"] = file_len(peakfile)
    df = pd.read_table(peakfile, index_col=False, header=None, sep="\t",
                       names=["chr", "start", "end"])
    df["Condition"] = cond
    df["Replicate"] = repl
    peaks = peaks.append(df, ignore_index=True)

# calculate peak sizes
peaks["width"] = peaks["end"] - peaks["start"]

# Counting peaks and bp
# -------------------------------------
# calculate total number of peaks and bp in peaks
peak_totals = peaks.groupby(["Condition", "Replicate"]).agg({
    "width": "sum", "start": "count"})
peak_totals = peak_totals.rename(columns={"width": "bp", "start": "peaks"})
# reset index for simpler plotting
peak_totals = peak_totals.reset_index()
melted_totals = pd.melt(peak_totals, value_vars=["peaks", "bp"],
                        id_vars=["Condition", "Replicate"],
                        var_name="Data", value_name="Count")

# ==============================================================================
# Plots
# ==============================================================================
# distribution of peak sizes
# -------------------------------------
print("Plotting peak size distribution")
plt.figure(figsize=(12, 12))
g = sns.catplot(data=peaks, kind="violin", x="Replicate",
                y="width", col="Condition", inner="quartile")
plt.savefig(os.path.join("Global", "peak-size-distribution.png"))
plt.close()

print("Plotting peak and bp counts")
fig = plt.figure(figsize=(12, 12))
sns.set(style="whitegrid")
g = sns.catplot(data=melted_totals, col="Data",
                x="Condition", y="Count", sharey=False,)
plt.savefig(os.path.join("Global", "peak-and-bp-counts.png"))
plt.close()
