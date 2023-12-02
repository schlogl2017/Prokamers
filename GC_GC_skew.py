#! usr/bin/env python
################################################################################
# Script to calculate the GC content, GC skew and GC variance in Prokariotes   #
# Paulo Sérgio Schlögl                                                         #
# Version 01 - 05/10/2022                                                      #
# License                                                                      #
################################################################################
#!/usr/bin/env python
import os
import glob
import math
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from fasta_parser import parse_fasta


def gc_content(seq, percent = False):
    """
    Calculate G+C content, return percentage (as float between 0 and 100).
    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content.  The percentage is calculated against
    the full length.

    Inputs:
        sequence - a string representing a sequence (DNA/RNA/Protein)

    Outputs:
        gc - a float number or a percent value representing the gc content of the sequence.

    """
    seq = seq.upper()
    seq_len = len(seq)
    g_c = seq.count("G") + seq.count("C")
    gc = g_c / seq_len
    if percent:
        return gc * 100
    return gc


def get_gc_content_slide_window(sequence, window, step, percent =False):
    """
    Calculate the GC (G+C) content along a sequence. Returns a dictionary-like
    object mapping the sequence window to the gc value (floats).
    Returns 0 for windows without any G/C by handling zero division errors, and
    does NOT look at any ambiguous nucleotides.

    Inputs:
        sequence - a string representing a sequence (DNA/RNA/Protein)
        window_size - a integer representing the length of sequence to be analyzed
                   at each step in the full sequence.
        step - a integer representing the length of overlapping sequence allowed.

    Outputs:
        gc- a dictionary-like object mapping the window size sequence to it gc values
            as floating numbers
    """
    gc = defaultdict(float)
    for start, seq in get_sequence_chunks(sequence, window, step):
        gc[start] += gc_content(seq, percent = percent)
    return gc


def difference_gc(total_gc, gc_dict):
    """
    Calculates the difference between the mean GC content of window i, and the
    mean chromosomal GC content, as Di = GC i − GC
    """
    # iterate through all keys, gc calculated values for the
    # genome chunk
    # get the difference between the chromosomal mean and the chunks
    # add the difference to the appropriate chunk
    d_i = {
        chunk: gc_cnt - total_gc
        for chunk, gc_cnt in gc_dict.items()
    }
    return d_i


def get_chromosomal_gc_variation(difference_dic):
    """
    Calculates the chromosomal GC variation defined as the log-transformed average
    of the absolute value of the difference between the mean GC content of each
    non-overlapping sliding window i and mean chromosomal GC content.
    chromosomal_gc_variantion = log(1/N*sum(|Di|), where N is the maximum number of
    non-overlapping window size in bp in a sliding windows.
    """
    # get the number of chunks
    n = len(difference_dic)
    arr = np.fromiter(
        difference_dic.values(),
        dtype=np.float32,
        count=len(difference_dic),
    )
    var = np.log(np.sum(np.abs(arr)) / n)
    return var


def get_sequence_skew(sequence):
    """
    Calculates the difference between the total number of
    occurrences of G and the total number of occurrences of C in
    the first i elements of the sequence. 

    Inputs:
        sequence - string representing the sequence     
    
    Outputs:
    
        skew - an array-like object that represents the GC skew of the sequence.
    
    Ex: 
    > get_sequence_skew('ACAACGTAGCAGTAGCAGTAGT')
    [0, 0, -1, -1, -1, -2, -1, -1, -1, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 2, 2]
    
    """
    # make the sequence upper case
    sequence = sequence.upper()
    # start the array
    skew = [0]
    # iterates to the sequence elements and it indexes
    for idx, element in enumerate(sequence):
        # check if element[i] is a G
        # if so add 1
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        # if the element[i] is a C
        # add to the array -1
        elif sequence[idx] == 'C':
            skew.append(skew[idx] - 1)
        else:
            # if it is not G or C add 0
            skew.append(skew[idx])
    return skew


def get_min_skew(seq_skew):
    """
    Calculates a position in a sequence minimizing the skew.
    
    Inputs:
        seq_skew - a list representing the sequence skew.    
    
    Outputs:
    
        min_skew - an array-like object that represents the pistion 
                   where the GC skew is the minimized in the sequence.    
    
    Example:
    get_minimum_skew('ACAACGTAGCAGTAGCAGTAGT')
    [5]
    """
    min_skew = []
    skew = seq_skew
    # get the minimized skew values
    m_skew = min(skew)
    # iterates to the length of the sequence
    # to get the index positions
    for idx in range(len(sequence) + 1):
        # if the position i has the same value 
        # as the minimum appende to the array
        if skew[idx] == m_skew:
            min_skew.append(idx)
    return min_skew


def get_gc_chr_var(filename, window, step):
    gcVar = defaultdict(float)
    for Id, sequence in parse_fasta(filename):
        geno_id = "_".join(Id.split("_")[:2])
        gc = get_gc_content(sequence)
        gc_win = get_GC_by_slide_window(sequence,
                                        get_gc_content,
                                        window,
                                        step)
        gc_dif = difference_gc(gc, gc_win)
        gcVar[geno_id] = get_chromosomal_gc_variation(gc_dif)
    return gcVar


def histplot(data, 
             title, 
             geno_id, 
             xlabel,
             dpi=100,
             forma="png",
             path=None, 
             as_save=False):
    path_to_save = f"{path}/{geno_id}_{xlabel}"
    sns.histplot(data=diff_gc,kde=True)
    dif_vals = diff_gc.values()
    plt.xlim((min(dif_vals) - 0.1, max(dif_vals) + 0.05))
    plt.xlabel(f"{xlabel}")
    plt.title(f"{geno_id}")
    plt.tight_layout()
    if as_save:
        plt.savefig(path_to_save, dpi=dpi, format=forma)


def skew(sequence):
    """
    Fx ( G ) − Fx ( C ) / Fx ( G ) + Fx ( C )
    """
    sequence = sequence.upper()
    g = sequence.count('G')
    c = sequence.count('C')
    skew = (g - c) / (g + c)
    return skew




















