import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pysam
import pickle
import itertools
import collections

def coverage_from_sam(file,targets=None):
    coverage = dict()
    if file.endswith('.sam'):
        samfile = pysam.AlignmentFile(file, "r")
    elif file.endswith('.bam'):
        samfile = pysam.AlignmentFile(file, "rb")
    else:
        raise ValueError("File must be .sam or .bam")
    for read in samfile:
        ref_name = read.reference_name
        if not ref_name:
            continue
        if targets is not None and ref_name not in targets:
            continue
        if ref_name not in coverage:
            coverage[ref_name] = np.zeros(samfile.get_reference_length(ref_name))
        for i in range(read.reference_start, read.reference_end):
            coverage[ref_name][i] += 1
    samfile.close()
    return coverage

def plot_coverage(coverage, start=None, end=None, figsize=(10, 6), title=None, save_path=None, log_scale=False):
    plt.figure(figsize=figsize)
    plt.plot(coverage, color='blue')
    if end is not None:
        plt.xlim(right=end)
    if start is not None:
        plt.xlim(left=start)
    if log_scale:
        plt.yscale('log')

    plt.xlabel('Position')
    plt.ylabel('Coverage')
    if title:
        plt.title(title)
    if save_path:
        plt.savefig(save_path)
    plt.show()

def iter_count(iterable):
    elem_count = itertools.count()
    collections.deque(zip(iterable, elem_count))
    return next(elem_count)