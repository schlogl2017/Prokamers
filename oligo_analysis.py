#!/usr/bin/env python
import os
import sys
import argparse
import time
import glob
import math
from collections import defaultdict
from itertools import product
import numpy as np
import pandas as pd
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
from fasta_parser import parse_fasta
import count_kmers


def base_stats(sequence, alphabet, as_count=False, as_dict=False):
    """Calculates de frequency or the number of bases in a sequence.
    
    Inputs:
    
        sequence - string representing the sequence
        alphabet - a alphabet (strings characters) that compound the string sequence
        as_count - boolean set as False
        as_dict - boolean set as False
    
    Output:
    
        counts - as default returns a numpy array as frequencies (floats) or
                 as a dictionary-like object
    
    Examples:
    
    > baseFreqs(seq, 'ACGT', asCounts = False, asDict = False)
    array([0.25, 0.25, 0.25, 0.25])

    as_count - True, returns a numpy array of counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = False)
    array([2, 2, 2, 2])

    as_dict - True and as_count as default (False) returns a dictionary as bases frequencies (float)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = False, asDict = True)
    {'A': 0.25, 'C': 0.25, 'G': 0.25, 'T': 0.25}

    as_count True and as_dict True, returns a dictionary as base counts (integer)
    > baseFreqs('ACGTACGT', 'ACGT', asCounts = True, asDict = True)
    {'A': 2, 'C': 2, 'G': 2, 'T': 2}
    """
    # make the sequence upper case
    seq = sequence.upper()
    # count all bases in sequence and collect as an array
    counts = np.array([seq.count(i) for i in alphabet])
    # if is onle the counts
    if as_count:
        freqs = counts
    # other wise as frequencies
    else:
        freqs = counts / sum(counts * 1.0)
    # or as a dictionary like object
    if as_dict:
        return dict(zip(alphabet, freqs)), len(seq)
    else:
        return freqs


def get_all_possible_kmers(alphabet, kmin, kmax):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet.

    Inputs:

        alphabet - a alphabet (strings characters) that compound the string sequence
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        kmers - list of all possible combinations of k-mers of length k with length
                between kmin and kmax.

    """
    kmers = [''.join(letters) for n in range(kmin, kmax + 1)
             for letters in product(alphabet, repeat=n)]
    return kmers


def get_strand_complement(sequence):
    """Returns the complement strand of the genome.
     
     Inputs:
        sequence - string representing the sequence   

    Outputs:
    
        sequence - string representing the complement of 
                   the string.    
    """
    # make the sequence upper case
    seq = sequence.upper()
    # table to change the complement characters
    change = str.maketrans('ACGT', 'TGCA')
    return seq.translate(change)


def get_reverse_complement(sequence):
    """
    Returns the reverse complement strand of the genome.

    Inputs:

        sequence - string representing the sequence.

    Outputs:

        reversed_complement_sequence - string representing the reversed
                                       sequence complement.
    """
    return get_strand_complement(sequence)[::-1]


def get_expected_higher_markov(kmer_list, kmer_counts):
    """
    Calculates the expected value for a list of kmers and their counts.
    
    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
    
    Output:
    
        expected - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
    
    The expected values are calculated as:
    'Expected = kmer[:-1] * kmer[1:] / kmer[1:-1]'
    """
    # initialize the container
    expected = defaultdict(int)
    # iterates through the list of kmers
    for kmer in kmer_list:
        # gets the suffix, prefix and middle kmers
        suf, pref, mid = kmer_counts[kmer[1:]], kmer_counts[kmer[:-1]], kmer_counts[kmer[1:-1]]
        # deal with sero division errors
        if mid == 0:
            expected[kmer] = expected.get(kmer, 0)
        else:
            # add the kmer and it expected values 
            expected[kmer] = expected.get(kmer, 0) + (pref * suf) / mid
    return expected


def get_variance(kmer_list, len_seq, kmer_expected):
    """
    Calculates the variance from a list of strings of length k.
    
    Inputs:
        kmer_list - list-like object representing all possible kmers of length k.
        kmer_expectd - a dictionary-like object mapping kmer to their calculated
                       expected values.
        len_seq - integer representing the length of the sequence where kmers were
                  counted.
    
    Outputs:
    
        variance - a dictionary-like object mapping kmer to their calculated
                   expectd variance.
                   
    Because the model for the count is the sum of N almost independent observations, 
    each with probability P(W), it can be well modeled as a binomial distribution, 
    with varianceThe variance is calculated as:
    E(C(W)) * (1 - E(C(W))/N)    
    """
    # gets the kmer length from the list 
    # of kmers
    k = len(kmer_list[0])
    # get all possible positions to count the kmers
    N = len_seq - k + 1
    # initialize the container
    variance = defaultdict(float)
    # iterates through the kmer list
    for kmer in kmer_list:
        # gets the expected value from the dictionary
        ex_val = kmer_expected[kmer]
        # deals with zero error division
        if ex_val == 0:
            # gets the zero value if expected is zero
            variance[kmer] = variance.get(kmer, 0.0)
        else:
            # calculates the variance
            var = ex_val * (1 - ex_val / N)
            # add the kmer and the variance values in the container
            variance[kmer] = variance.get(kmer, 0.0) + var
    return variance    


def get_standard_deviation(variance):
    """
    Calaculates the standard deviation from the kmers expected values.
    
    Inputs:
        variance - a dictionary-like object mapping kmer to their calculated
                   expectd variance.
    
    Outputs:
        
        std - a dictionary-like object mapping kmer to their calculated
                   expectd std.
    
    The variance is calculated as:
    sigma(W) = sqrt(Expected) * (1 - Expected/len(seq) -k + 1))
    """
    # initialize the container
    std = defaultdict(float)
    # iterates through the kmer list
    for kmer in variance:
        # deal with zero divisions errors
        if variance[kmer] == 0.0:
            std[kmer] = std.get(kmer, 0.0)
        else:
            # add the kmer and the calculated std in the
            # container
            sd = math.sqrt(variance[kmer])
            std[kmer] = std.get(kmer, 0.0) + sd
    return std   


def z_scores(kmer_exp, kmer_counts, std):
    """
    Calculates the z scores to under/over represented kmers from a sequence.
    The score is calculaated as:
    
    Z(W) = (C(W) – E(C(W))) / sigma(W), where 
    C(w) - observed values
    E(C(w)) - represents the expected value from a kmer
    sigma - represents the standard deviation
    
    Inputs:
        kmer_exp - dictionary-like object mapping kmer of length k to their 
                   calculated expected values.
        kmer_counts - dictionary-like object mapping kmer to their counts. The
                      kmer lengths must be between kmin (kmx-2) and kmax.
        std - a dictionary-like object mapping kmer to their calculated
                   expectd std.
    
    Outputs:
        z_scores - dictionary-like object mapping kmer to their z_scores.
    """
    # initialize the container
    z_scores = defaultdict(float)
    # iterates through the kmer keys
    for kmer in kmer_exp:
        # gets the kmer std value
        sd = std[kmer]
        # deals with zero error division
        if sd == 0.0:
            z_scores[kmer] = z_scores.get(kmer, 0.0)
        else:
            # calculates the z score and add 
            # the kmer and the z score values to the container
            z = (kmer_counts[kmer] - kmer_exp[kmer]) / sd
            z_scores[kmer] = z
    return z_scores 
    

# Calculate tetranucleotide values for each input sequence
def calc_oligo_usage(fasta_files, kmin, kmax):
    """
    Calculate the oligonucleotide and it corresponding Z-score 
    for each input sequence.
    
    Inputs:
        fasta_files - list of fasta files to be analyzed
        kmin - integer representing the minimum length of the
               oligonucleotides.
        kmax - integer representing the maximum length of the
               oligonucleotides.
    Outputs:
        orgs_oligos -  a dictionary of dictionaries, mapping the
        species/genus to the kmer counts (kmer as keys and counts as values
        in the nested dictinonary).
        Ex: 
        {'Thermodesulfatator': defaultdict(float,
             {'AAAA': -4.432324117689154,
              'AAAC': 3.9058350768342853,
              'AAAG': -0.7145235274286397,
              'AAAT': 4.004908092189418,
              'AACA': 9.363779538974663,
              'AACC': -8.309560927685952,
              'AACG': -10.038301152300557,
              'AACT': 8.322627567820732,
              'AAGA': 12.962109060036239,
              'AAGC': -11.595188214766853,
              'AAGG': -5.770556875526479,
              'AAGT': 3.673330047942202,
              'AATA': 12.445578038351664,
              'AATC': -11.555648757484228,
              'AATG': -8.120108540533565,
              'AATT': 3.504440216329279,...}}
    """
    # need a list of kmer of length k to analysis
    kmer_list = get_all_possible_kmers('ACGT', kmax, kmax)
    orgs_oligos = {}
    for filename in fasta_files:
        for Id, seq in parse_fasta(filename):
            N = len(seq)
            org = filename.split('/')[2]
            print(f"Calculating oligonucleotide z-scores for {org}")
            # count mers from kmin (kmax -2) to kmax
            kmer_counts = count_kmers.count_kmers(seq, kmin, kmax)
            # get the expected kmer values from the kmers of length kmax
            # Expected = kmer[:-1] * kmer[1:] / kmer[1:-1]
            kmer_exp = get_expected_higher_markov(kmer_list, kmer_counts)
            # calculate the variance
            # E(C(W)) * (1 - E(C(W))/N), N = len(seq)-k+1
            kmer_var = get_variance(kmer_list, N, kmer_exp)
            # calculate the standard deviation
            # sigma(W) = sqrt(Expected) * (1 - Expected/len(seq) -k + 1))
            kmer_sigma = get_standard_deviation(kmer_var)
            # calculate  the z-scores
            # Z(W) = (C(W) – E(C(W))) / sigma(W)
            kmer_zscr = z_scores(kmer_exp, kmer_counts, kmer_sigma)
            orgs_oligos[org] = kmer_zscr
    return orgs_oligos


def get_obs_exp_ratio(kmer_counts, kmer_expected, kmer_list):
    """ 
    Calculates the odds ratio of the kmers observed values from
    the expected values.
    """
    ratio = defaultdict(float)
    for kmer in kmer_list:
        o = kmer_counts[kmer]
        e = kmer_expected[kmer]
        if e == 0.0:
            ratio[kmer] = 0.0
        else:
            ratio[kmer] = o/e
    return ratio
    
    
# Write the set of tetranucleotide frequency Z scores to a plain text
# tab-separated table, in the output directory
def write_oligo_usage(dir_out, tab_filename, oligo_z):
    """ Writes the Z score for each tetranucleotide to a plain text
        tab-separated table with one row for each tetranucleotide, and one
        column for each input sequence.
        - filename is the location of the file to which the tetranucleotide
              Z-scores should be written
        - tetra_z is a dictionary containing the tetranucleotide Z-scores
              for each input sequence
    """
    try:
        fh = open(os.path.join(dir_out, tab_filename), 'w')
        print(f"Writing Z-scores to {fh.name}")
    except:
        print("Could not open %s for writing (exiting)")
        sys.exit(1)
    orgs = sorted(oligo_z.keys())
    kmers = sorted(list(oligo_z.values())[0].keys())
    # Write headers
    print(f"# oligo_usage.py {time.asctime()}")
    print("# Oligonucleotides frequency Z-scores")
    print('kmer\t' + '\t'.join(orgs), file=fh)
    for kmer in kmers:
        outstr = [kmer] + ["%.2f" % tetra_z[org][kmer] for org in orgs]
        print('\t'.join(outstr), file=fh)
    fh.close()    
    
    
def do_calc_oligo_usage_window(sequence, kmin, kmax, window, step):
    # 
    # get all possible k length kmers
    kmer_list = get_all_possible_kmers('ACGT', kmax, kmax)
    orgs_oligos = defaultdict(dict)
    index = []
    difs = []
    N = len(seq)
    # count all kmers of length kmin - kmax
    obs = count_kmers.count_kmers(seq, kmin, kmax)
    # get the expectancy with kmes of length kmax
    exp = get_expected_higher_markov(kmer_list, obs)
    # get the observed/expected ratio
    ob_exp_ratios = get_obs_exp_ratio(obs, exp, kmer_list)
    # chunk the sequence in chunks of length window and overlap step
    chunks = list(get_sequence_chunks(seq, window, step, True))
    for chunk in chunks:
        # get starts of the chunks
        index.append(chunk[1])
        # count mers from kmin (kmax -2) to kmax
        obs_chunk = count_kmers.count_kmers(chunk[0], kmin, kmax)
        # get the expected kmer values from the kmers of length kmax
        # Expected = kmer[:-1] * kmer[1:] / kmer[1:-1]
        exp_chunk = get_expected_higher_markov(kmer_list, obs_chunk)
        # get the difference observed/expected -  ratio
        dif_sum = 0
        for kmer in kmer_list:
            o = obs_chunk[kmer]
            e = exp_chunk[kmer]
            r = ob_exp_ratios[kmer]
            if e == 0.0:
                dif_sum += 0.0
            else:
                dif_sum += abs((o/e) - r)
        difs.append(dif_sum)
    # calculate the z-scores
    z_scores = []
    arr = np.asarray(difs)
    m = np.mean(arr)
    sd = np.std(arr)
    for i, d in enumerate(difs):
        if sd == 0.0:
            z_scores.append(0.0)
        else:
            z = (d - m) / sd
            z_scores.append(z)
    return index, z_scores    


def get_kmer_freq(seq_len, kmer_counts, kmax):
    freq = defaultdict(float)
    gouped_kmer = toolz.groupby(len, c)
    for i in range(2, kmax + 1):
        for kmer in gouped_kmer[i]:
            freq[k] = kmer_counts[kmer]/(seq_len - i + 1)
    return freq


def do_kmer_counts(sequence, alphabet, kmin, kmax, as_subset=None, as_rev_seq=False):
    # if to count the forward and reverse sequence
    if as_rev_seq:
        sequence = sequence + get_reverse_complement(sequence)
    # if you want to count only a subset of the sequence
    if as_subset != None:
        sequence = sequence[subset[0]:subset[1]]
    # count kmers
    kmer_counts = count_kmers.count_kmers(sequence, kmin, kmax)
    return kmer_counts
    
    
def do_kmer_counts_window(sequence, kmin, kmax, win_len, step):
    """
    Count the kmers of length kmin to kmax from substrings of
    length win_len with overlaps of length step.
    
    Inputs:
        sequence - a string representing a sequence
        kmin - integer representing the minimum length of the
               oligonucleotides.
        kmax - integer representing the maximum length of the
               oligonucleotides.
        win_len - a integer representing the lengths of subsequences.
        step - a integer representing the lengths of overlapping bases.
    
    Outputs:
        kmer_windows - a dictionary of dictionaries mapping the subsequences
                       index to the kmer counts in that subsequence.
    
    Ex:
    defaultdict(dict, {(0, 10): defaultdict(int,
                         {'AA': 1,'AC': 1,'AG': 2,'AT': 0,'CA': 1,'CC': 2,
                          'CG': 0,'CT': 0,'GA': 2,'GC': 0,'GG': 0,'GT': 0,
                          'TA': 0,'TC': 0,'TG': 0,'TT': 0,'AAA': 0,'AAC': 0,
                          'AAG': 0, 'AAT': 0, 'ACA': 0, 'ACC': 1, 'ACG': 0, 
                          'ACT': 0, 'AGA': 2, 'AGC': 0, 'AGG': 0, 'AGT': 0, 
                          'ATA': 0, 'ATC': 0,...}}
    """
    # create subsequences
    chunks = get_sequence_chunks(sequence, win_len, step=step, as_overlap = True)
    kmer_windows = defaultdict(dict)
    # for each subsequence count the kmers
    for chunk in chunks:
        seq = chunk[0]
        start = chunk[1]
        end = chunk[2]
        # put the keys as start-end and the kmer counts in the dictionary
        kmer_windows[start, end] = count_kmers.count_kmers(chunk[0], kmin, kmax)
    return kmer_windows    
    
    
def get_sequence_chunks(
        sequence,
        window,
        step=1,
        as_overlap = True
):
    """
    Function to sub-sequence in overlap windows of length window and overlaping as length step, 
    other wise non overlapping.

    Inputs:
        sequence - a string representing a DNA sequence.
        as_overlap - boolean that represents the overlap length.
        step - a integer representing the lengths of overlapping bases.
        window - a integer representing the lengths of subsequences.
               Default = 1

    Outputs:
        subseq - a string representing a slice of the sequence with or without
                 overlapping characters.
    
    Examples:
    # overlapping
    > list(get_sequence_chunks('ACCAGTGGATTGAGGAGATATAG', 10,  1))
    [('ACCAGTGGAT', 0, 10),
     ('CCAGTGGATT', 1, 11),
     ('CAGTGGATTG', 2, 12),
     ('AGTGGATTGA', 3, 13),
     ('GTGGATTGAG', 4, 14),
     ('TGGATTGAGG', 5, 15),
     ('GGATTGAGGA', 6, 16),
     ('GATTGAGGAG', 7, 17),
     ('ATTGAGGAGA', 8, 18),
     ('TTGAGGAGAT', 9, 19),
     ('TGAGGAGATA', 10, 20),
     ('GAGGAGATAT', 11, 21),
     ('AGGAGATATA', 12, 22),
     ('GGAGATATAG', 13, 23)]
    # non-overlapping
    > list(get_sequence_chunks('ACCAGTGGATTGAGGAGATATAG', 10, False))
    [('ACCAGTGGAT', 0, 10), ('TGAGGAGATA', 10, 20)]
    """
    sequence = sequence.upper()

    len_range = range(0, len(sequence) - window + 1, step)
    
    if not as_overlap:
        # overlap sequence length
        len_range = range(0, len(sequence) - window + 1, window)

    # iterates to the overlap region
    for i in len_range:
       # creates the substring
        subseq = sequence[i:i + window]
        start = i
        end = i + window
        yield subseq, start, end    
    
    
def do_get_expected_higher_markov(sequence, alphabet, kmer_list, kmax, as_subset=None, as_rev_seq=False):
    if as_rev_seq:
        sequence = sequence + get_reverse_complement(sequence)
    if as_subset != None:
        sequence = sequence[subset[0]:subset[1]]
    kmer_counts = count_kmers.count_kmers(sequence, kmax-2, kmax)
    expec = get_expected_higher_markov(kmer_list, kmer_counts)
    return {k:v for (k,v) in expec.items() if v > 0.0}    
    
    
    
def do_calculate_oligo_usage(filename, kmer_list, kmin, kmax, win_len, step, subset=False):
    for Id, sequence in fasta_parser.parse_fasta(filename):
        if subset:
            sequence = sequence[subset[0]:subset[1]]
        # returns three lists of kmes,start indexes of the seq chunks and z_scores
        return calculate_oligo_usage(sequence, kmer_list, kmin, kmax, win_len, step)    
        
        
        
        
     
        
        
        
        
        
        
            
