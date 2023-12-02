#! usr/bin/env python
#  Functions to work with kmer statistics
#  Autor: Paulo Sérgio Schlögl
#  24/06/2019
#  Last Modification: 27/09/2022
import re
import math
from itertools import product, islice, pairwise
import csv
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.stats import norm


def get_kmer_sub_strs(kmer):
    pref = kmer[:-1]
    mid = kmer[1:-1]
    suf = kmer[1:]
    return pref, mid, suf
    

def get_data_from_csv(filename):
    data = defaultdict(float)
    with open(filename, 'r') as fh:
        csv_data = csv.reader(fh)
        for row in csv_data:
            d1, d2 = row[0], float(row[1])
            data[d1] = data.get(d1, 0.0) + d2
    return data


def group_by_length(kmer_counts, k):
    selected =  defaultdict(int)
    for kmer, count in kmer_counts.items():
        if len(kmer) == k:
            selected[kmer] = selected.get(kmer, 0) + count
    return selected


def calculate_expected_kmers(filename, k):
    kmer_counts = defaultdict(int)
    with open(filename, "r") as csvfile:
        reader = csv.reader(csvfile)
        # escape the header of the csv file
        header = csvfile.readline()
        for line in reader:
            kmer = line[0]
            count = int(line[1])
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + count
    selected = group_by_length(kmer_counts, k)
    kmer_list = selected.keys()
    expected = get_expected_higher_markov(kmer_list, kmer_counts)
    dfc = pd.DataFrame(selected.items(), columns=['kmer', 'observed'])
    dfe = pd.DataFrame(expected.items(), columns=['kmer', 'expected'])
    df = dfc.merge(dfe, on = "kmer")
    return df  


def merge_jelly_dump_kmer(name, filenames, rg):
    grouped_files = []
    for filename in filenames:
        m = get_assembly_ids(filename, rg)
        if name == m:
            grouped_files.append(filename)
    dfs = [pd.read_table(filename, 
                         sep="\s+", 
                         names=["kmer", 
                                "observed"]) for filename in grouped_files]
    df = pd.concat(dfs)
    df.index = df["kmer"].str.len()
    df = df.sort_index(ascending=True).reset_index(drop=True)
    mask = df['kmer'].str.len() == 1
    length = df[mask]['observed'].sum()
    df['observed_freq'] = df["observed"] / length
    return df


def odds_ratio(kmer_list, kmer_counts, kmer_expected):
    """
    Need to test
    """
    odds = defaultdict(float)
    for kmer in kmer_list:
        odd = kmer_counts[kmer] / kmer_expected[kmer]
        odds[kmer] = odds.get(kmer, 0.0) + odd
    return odds


def odds_ratio_kmer(kmer, kmer_counts):
    kmer_list = list(kmer)
    base_mul = 1
    for base in kmer_list:
        base_mul *= kmer_counts[base]
    odds = kmer_counts[kmer] / base_mul
    return odds


def make_odds_ratio(kmer_counts, kmer_list):
    odds_ratio = defaultdict(float)
    for kmer in kmer_list:
        odd = odds_ratio_kmer(kmer, kmer_counts)
        odds_ratio[kmer] = odds_ratio.get(kmer, 0.0) + odd
    return odds_ratio


def kmer_mean(kmer_freqs, alphabet, k):
    """
    x^ = 1/N*sum(xi); N = 4**k, xi = freq kmer i
    """
    total = sum(kmer_freqs.values())
    len_alph = len(alphabet)
    mean = total/len_alph**k
    return mean


def kmer_slice(kmer, k):
    """ >>> list(kmer_slice('ATCGA', 2))
        [('A', 'T'), ('T', 'C'), ('C', 'G'), ('G', 'A')]
    """
    itr = iter(kmer)
    res = tuple(islice(itr, k))
    if len(res) == k:
        yield res
    for b in itr:
        res = res[1:] + (b,)
        yield res


def slicing(kmer, k):
    for km in kmer_slice(kmer, k):
        kmer_sliced = "".join(km)
        yield kmer_sliced


def get_kmer_frequency(kmer_data):
    kmer_freq = defaultdict(float, [(kmer, 0.0) for kmer in kmer_data])
    total_freq = 0
    for kmer in kmer_freq:
        kmer_freq[kmer] = kmer_data[kmer]
        total_freq += kmer_data[kmer]
    return {k: (cnt/total_freq) for k, cnt in kmer_freq.items()}


def get_freqs(kmer_counts, alphabet, k):
    freqs = defaultdict(float)
    for kmer in get_all_possible_kmers(alphabet, k, k):
        freqs[kmer] = freqs.get(kmer, 0.0) + kmer_counts[kmer]
    total = sum(freqs.values())
    return {kmer: count/total for kmer, count in freqs.items()}


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




def get_expected_kmer_by_zom(kmer_list, base_freqs, len_seq):
    zom = defaultdict(float)
    for kmer in kmer_list:
        ze = expected_kmer_by_zom(kmer, base_freqs, len_seq)
        zom[kmer] = zom.get(kmer, 0.0) + ze
    return zom


def second_order_markov(kmer, kmer_counts):
    p = kmer_counts[kmer[:3]]
    s = kmer_counts[kmer[1:4]]
    m1 = kmer_counts[kmer[2:5]]
    m2 = kmer_counts[kmer[3:]]
    d1 = kmer_counts[kmer[1:3]]
    d2 = kmer_counts[kmer[2:4]]
    d3 = kmer_counts[kmer[3:5]]
    exp = (p * s * m1 * m2) / (d1 * d2 * d3)
    return exp
    

def third_order_markov(kmer, kmer_counts):
    p = kmer_counts[kmer[:4]]
    s = kmer_counts[kmer[1:5]]
    s2 = kmer_counts[kmer[2:]]
    d1 = kmer_counts[kmer[1:4]]
    d2 = kmer_counts[kmer[2:5]]
    exp = (p * s * s2) / (d1 * d2)
    return exp   


def get_assembly_ids(string, regex):
    assm = regex.search(string).group()
    return assm
# rg =  re.compile(r"(GC[AF]_\d+\.\d)")


def get_expected_higher_markov(kmer_data, kmer_counts):
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
    expected = defaultdict(float)
    variance = defaultdict(float)
    for kmer in kmer_data:
        expected[kmer] = expected.get(kmer, 0.0)
        variance[kmer] = variance.get(kmer, 0.0)
        pref, mid, suf = get_kmer_sub_strs(kmer)
        p = kmer_counts[pref]
        m = kmer_counts[mid]
        s = kmer_counts[suf]
        if m == 0:
            expected[kmer] = 0.0
            variance[kmer] = 0.0
        else:
            exp = (p * s) / m
            expected[kmer] = round(exp, 5)
            variance[kmer] = round((exp * (m - p) * (m - s)) / (m ** 2), 5)
    return expected, variance


def get_variance_new(kmer_list, kmer_counts, len_seq, kmer_expected):
    """
    Calculates the variance from a list of strings of length k.

    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_expect - a dictionary-like object mapping kmer to their calculated
                       expected values.

    Outputs:

        variance - a dictionary-like object mapping kmer to their calculated
                   expect variance.
    kmer_var <- kmer_exp * (m - p) * (m - s) / (m^2)
    Because the model for the count is the sum of N almost independent observations,
    each with probability P(W), it can be well modeled as a binomial distribution,
    with varianceThe variance is calculated as:
    E(C(W)) * (1 - E(C(W))/N)
    """
    k = len(kmer_list[0])
    N = len_seq - k + 1
    variance = defaultdict(float)
    for kmer in kmer_list:
        suf = kmer_counts[kmer[1:]]
        pref = kmer_counts[kmer[:-1]]
        mid = kmer_counts[kmer[1:-1]]
        ex_val = kmer_expected[kmer]
        if ex_val == 0:
            variance[kmer] = variance.get(kmer, 0.0)
        else:
            #var = exp * (m - p) * (m - s) / (m^2)
            var = (ex_val * (mid - pref) * (mid - suf)) / (mid**2)
            variance[kmer] = variance.get(kmer, 0.0) + var
    return variance


def get_variance(kmer_list, len_seq, kmer_expected):
    """
    Calculates the variance from a list of strings of length k.

    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        kmer_expect - a dictionary-like object mapping kmer to their calculated
                       expected values.

    Outputs:

        variance - a dictionary-like object mapping kmer to their calculated
                   expect variance.

    Because the model for the count is the sum of N almost independent observations,
    each with probability P(W), it can be well modeled as a binomial distribution,
    with varianceThe variance is calculated as:
    E(C(W)) * (1 - E(C(W))/N)
    """
    k = len(kmer_list[0])
    N = len_seq - k + 1
    variance = defaultdict(float)
    for kmer in kmer_list:
        ex_val = kmer_expected[kmer]
        if ex_val == 0:
            variance[kmer] = variance.get(kmer, 0.0)
        else:
            var = ex_val * (1 - ex_val / N)
            variance[kmer] = variance.get(kmer, 0.0) + var
    return variance


def get_standard_deviation(variance):
    """
    Calculates the standard deviation from the kmers expected values.

    Inputs:
        variance - a dictionary-like object mapping kmer to their calculated
                   expected variance.

    Outputs:

        std - a dictionary-like object mapping kmer to their calculated
                   expected std.

    The variance is calculated as:
    sigma(W) = sqrt(variance))
    """
    # initialize the container
    std = defaultdict(float)
    # iterates through the kmer keys
    for kmer in variance:
        # deals with zero error division
        if variance[kmer] == 0.0:
            std[kmer] = std.get(kmer, 0.0)
        else:
            # calculates the standard deviation and add
            # the kmer and the standard deviation values to the container
            sd = math.sqrt(variance[kmer])
            std[kmer] = std.get(kmer, 0.0) + sd
    return std


def get_z_scores(kmer_data, expected_kmers, std):
    """Calculates the z_score of all palindromes.
    Input:
    palindrome_lst = list of palindromes (str)
    counts = dictionary of kmer counts
    expected = dictionary of kmer expected values
    length_sequence = length of sequence (int)
    Output:
    z_score dictionary where key are palindromes and values are the calculated z_score (float)
    The z_scores are calculated as:
        Z(W) = (C(W)) - E(C(W)) / sigma(W)
    And sigma as:
        sigma(W) = sqrt(E(C(W))) * (1 - E(C(W)/N))
    """
    z_score = defaultdict(float)
    for kmer in kmer_data:
        if expected_kmers[kmer] == 0.0:
            z_score[kmer] = 0.0
        else:
            sigma = std[kmer]
            obs = kmer_data[kmer]
            exp = expected_kmers[kmer]
            z = round((obs - exp) / sigma, 5)
            z_score[kmer] = z_score.get(kmer, 0.0) + z
    return z_score


def get_p_values(z_scores_kmers):
    """
    Calculates the p value for all kmers.
    The calculation is done as:
    over represented: P(z > t) = erfc(t/sqrt(2))/2
    under represented: P(z > t) = erfc(-t/sqrt(2))/2
    t: threshold

    Inputs:
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.

    Outputs:
        p_vals - dictionary-like object mapping kmer to their p values.
    """
    # initialize the container
    p_vals = defaultdict(float)
    # iterates through the kmer keys
    for kmer in z_scores_kmers:
        # calculates the p values to under represented
        # kmers (negative z scores)
        # add the kmer and p values to the container
        if z_scores_kmers[kmer] < 0.0:
            under = math.erfc(-z_scores_kmers[kmer] / math.sqrt(2)) / 2
            p_vals[kmer] = p_vals.get(kmer, 0.0) + under
        else:
            # add the kmer and p values to the container to over represented
            # and all other kmers
            others = math.erfc(z_scores_kmers[kmer] / math.sqrt(2)) / 2
            p_vals[kmer] = p_vals.get(kmer, 0.0) + others
    return p_vals


def get_scipy_p_values(z_scores_kmers):
    """
    Calculates the p value for all kmers.
    The calculation is done as:
    over represented: P(z > t) = erfc(t/sqrt(2))/2
    under represented: P(z > t) = erfc(-t/sqrt(2))/2
    t: thresholder
    
    Inputs:
        z_scores_kmers - dictionary-like object mapping kmer to their z_scores.
        
    Outputs:
        p_vals - dictionary-like object mapping kmer to their p values.
    """
    p_vals = defaultdict(float)
    for kmer in z_scores_kmers:
        p = norm.sf(abs(z_scores_kmers[kmer])) * 2
        p_vals[kmer] = p_vals.get(kmer, 0.0) + p
    return p_vals


def get_e_values(p_vals):
    """
    Calculates the variance from a list of strings of length k.

    Inputs:
        kmer_list - list-like object representing all kmer counted in a sequence.
                    The kmers have a length k.
        p_vals - a dictionary-like object mapping kmer to their calculated
                       p values.

    Outputs:

        e_values - a dictionary-like object mapping kmer to their calculated
                   e values.
    """
    # number of tested hypothesis
    hyp_num = len(p_vals)
    # initialize the container
    e_values = defaultdict(float)
    # iterates through the kmer list
    for kmer in p_vals:
        # gets the p values from the input container
        hyp = hyp_num * p_vals[kmer]
        # calculates the e values and add the kmer
        # and the e values to the container
        e_values[kmer] = e_values.get(kmer, 0.0) + hyp
    return e_values


def significance_indexes(kmer_list, e_values):
    """ 
    The significance index is simply a negative logarithm
    conversion of the E-value (in base 10).
    
    """
    sig_idx = defaultdict(float)
    for kmer in kmer_list:
        # calculates the significance indexes
        sig = -math.log10(e_values[kmer])
        sig_idx[kmer] = sig_idx.get(kmer, 0.0) + sig
    return e_values


def get_kmer_stats(kmer_data,
                   observed_freq,
                   expected, 
                   expected_freq,
                   z_scores, 
                   e_vals, 
                   p_vals):
    data = []
    for kmer in kmer_data:
        data.append((kmer,
                     kmer_data[kmer],
                     observed_freq[kmer],
                     expected[kmer],
                     expected_freq[kmer],
                     z_scores[kmer],
                     e_vals[kmer],
                     p_vals[kmer]))
    df = pd.DataFrame(data,
                      columns=['kmer',
                               'observed',
                               'observed_freq',
                               'expected',
                               'expected_freq',
                               'z_score',
                               'e_value',
                               'p_value']).sort_values(by='z_score').reset_index(drop=True)
    df['Rank'] = df['z_score'].rank(method='max')
    return df


def get_inter_kmer_distance(seq, kmer):
    lk = len(kmer)
    starts = []
    for i, _ in enumerate(seq):
        subseq = seq[i:i+lk]
        if subseq == kmer:
            starts.append(i)
    dist = [end - start for start, end in pairwise(starts)]
    return dist, starts


def kmer_modelling(kmer_list, 
                   kmer_count, 
                   names, 
                   seq_len_dic, 
                   data=None):
    dict_data = defaultdict(dict)
    df = pd.DataFrame()
    for name in names:
        dict_data[name] = dict_data.get(name, {})
        # calculate the kmer statistics
        counts = kmer_count[name]
        exp_kmer = get_expected_higher_markov(kmer_list, counts)
        kmer_freq = get_kmer_frequency(kmer_list, counts)
        kmer_var = get_variance(kmer_list, 
                                seq_len_dic[name], 
                                exp_kmer)
        kmer_std = get_standard_deviation(kmer_var)
        kmer_zsc = get_z_scores(kmer_list, 
                                counts, 
                                exp_kmer, 
                                seq_len_dic[name])
        kmer_pval = get_p_values(kmer_zsc)
        kmer_eval = get_e_values(kmer_list, kmer_pval)
        kmer_sig_idx = significance_indexes(kmer_list, kmer_eval)
        kmer_counts = {
                       k:cnt for k, cnt in counts.items() if k in  kmer_list                      
                      }
        # select the apropriate statistics to return
        if data == "z_scores":
            dict_data[name] = get_scores_evalues(kmer_list, kmer_zsc, kmer_eval)
            dfz = pd.DataFrame(dict_data).T.rename(columns={0: "z_score", 1: "e_values"})
            df =  dfz.T.reset_index().rename(columns={'index':'kmer'})
        if data == "freq":
            dict_data[name] = kmer_freq
            dff = pd.DataFrame(dict_data).reset_index()
            df = dff.rename(columns={'index':'kmer'})
        if data == "expected":
            dict_data[name] = exp_kmer
            dfe = pd.DataFrame(dict_data).reset_index()
            df = dfe.rename(columns={'index':'kmer'})
        if data == "sig_idx":
            dict_data[name] = kmer_sig_idx
            dfsg = pd.DataFrame(dict_data).reset_index()
            df = dfe.rename(columns={'index':'kmer'})
        elif data == None:
            dict_data[name] = kmer_counts
            dfc = pd.DataFrame(dict_data).reset_index()
            df = dfc.rename(columns={'index':'kmer'})
    return df


def get_scores_evalues(kmer_list, zscores_dict, evalues_dict):
    z_e = defaultdict(dict, [(k,{}) for k in kmer_list])
    for kmer in kmer_list:
        for (k,v), (k2,v2) in zip(zscores_dict.items(), evalues_dict.items()):
            if kmer == k and kmer == k2:
                z_e[kmer] = (v, v2)
    return z_e  


def get_data_from_csv(filename):
    data = defaultdict(float)
    with open(filename, 'r') as fh:
        csv_data = csv.reader(fh)
        for row in csv_data:
            d1, d2 = row[0], float(row[1])
            data[d1] = data.get(d1, 0.0) + d2
    return data


def get_seq_length(filenames):
    lengths = defaultdict(int)
    for filename in filenames:
        name = filename.split("/")[-1].strip("fasta.sizes")
        lengths[name] = lengths.get(name, 0)
        with open(filename, "r") as fh:
            for line in fh:
                line_spl = line.strip().split("\t")
                lengths[name] = int(line_spl[1])
    return lengths


def seq_lengths(filenames):
    lengths = defaultdict(int)
    for filename in filenames:
        with open(filename, 'r') as fh:
            for line in fh:
                name, length = line.split("\t")
                lengths[name] = lengths.get(name, 0) + int(length)
    return lengths


def get_seq_len_from_csv(filename):
    lengths = defaultdict(float)
    with open(filename_len, "r") as fh:
        header = fh.readline()
        for line in fh:
            line_str = line.strip()
            name, pb = line_str.split(",")[0], float(line_str.split(",")[-1])
            lengths[name] = lengths.get(name, 0.0) + pb
    return lengths


def expect_palidromes(gc_seq, gc_pal, len_seq, len_pal):
    """
    Species-specific Typing of DNA Based on Palindrome Frequency Patterns
    E STELLE Lamprea-Burgunder, P HILIPP Ludin, and P ASCAL Mäser
    """
    one = (gc_seq/2)**(gc_pal*lenpal)
    two = ((1 - gc_seq)/2) ** (1-gc_pal) *  len_pal
    exp_pal = one * two * len_seq
    return exp_pal


def R_pal(pal_count, exp_pal):
    return pal_count/exp_pal


def ouv(kmer_count, kmer_expected):
    ouv = 0
    for kmer in kmer_expected:
        ouv += kmer_count[kmer] - kmer_expected[kmer]
    return ouv


def get_names(filename):
    filename_spl = filename.split("/")
    superkindom, group, namer = filename_spl[2], filename_spl[3], filename_spl[-1]
    # for bac/archaea/viruses
    if len(namer) == 31:
        name = namer.replace("_km_1_10_chr.csv", "")
    # for mito/plastids
    else:
        name = namer.replace("_km_1_6_chr.csv", "")
    return superkindom, group, name


def get_data_from_counts(kmer_counts, k, probs=False):
    data = defaultdict(float)
    for kmer in get_all_possible_kmers(iupac_dna, k, k):
        data[kmer] = data.get(kmer, 0.0) + kmer_counts[kmer]
    if probs:
        return {k:v/sum(data.values()) for k, v in data.items()}
    return data


def kmer_count_bases(kmer):
    kbases = defaultdict(int)
    for base in kmer:
        kbases[base] = kbases.get(base, 0) + 1
    return kbases


def dinuc_relative_abundance(dinuc_data, kmer_freqs):
    """
    rel_ab <= 0.78 and rel_ab >= 1.23 as suitable benchmarks for 
    assessing whether a dinucleotide is signiﬁcantly over-or 
    under-represented in a DNA sequence.
    (Karlin and Cardon, 1994; Karlin and Ladunga, 1994)
    """
    kmer_rel_ab = defaultdict(float)
    for kmer in dinuc_data.keys():
        b1 = kmer_freqs[kmer[0]]
        b2 = kmer_freqs[kmer[0]]
        rel_ab = dinuc_data[kmer] / (b1 * b2)
        kmer_rel_ab[kmer] = kmer_rel_ab.get(kmer, 0.0) + rel_ab
    return kmer_rel_ab


def select_pal_sum_counts(kmer_counts):
    """
    Select and add up the frequencies of the palindromic kmers,
    found in a DNA sequence.
    """
    # select the palindrome tuples
    # [('TATA', 'ATAT'),..('CTAG', 'GATC')]
    palindromes = []
    for km in kmer_counts:
        kmr = get_reverse_complement(km)
        if km == kmr:
            palindromes.append((km, kmr[::-1]))
    # now it selects the palindromes that are
    # lexically classified
    # {'AATT', 'ACGT', 'AGCT', 'ATAT', ...,'CTAG'}
    # as ATAT cames first than TATA..
    set_pal = set()
    for tup in palindromes:
        if tup[0] < tup[1]:
            set_pal.add(tup[0])
        elif tup[1] < tup[0]:
            set_pal.add(tup[1])
    # Now this part add up the frequencies of the 
    # palindrome pairs Ex, TATA + ATAT freqs
    sum_pal = defaultdict(float)
    for km in set_pal:
        comp = get_strand_complement(km)
        sum_freq = kmer_counts[km] + kmer_counts[comp]
        sum_pal[km] = sum_pal.get(km, 0.0) + sum_freq
    # Returns a new dictionary with counts of the
    # palindromes and the non-palindromes kmers
    # select only the lexically palindromes
    pal = set([item for sublist in palindromes for item in sublist])
    #print(pal)
    # returns a new dictionary with palindromes ad up and
    # all other no palindromic kmers
    new_dict = {k: v for k, v in kmer_counts.items() if k not in pal}
    new_dict.update(sum_pal)
    return new_dict


def get_wildcard_kmers(kmer):
    N = ["A", "C", "G", "T"]
    wildcards = []
    if len(kmer) == 3:
        for n in N:
            knm = kmer[0] + n + kmer[2]
            wildcards.append(knm)
    elif len(kmer) == 4:
        for n in N:
            # XYZW
            wildcards.append(kmer[0] + n + kmer[2]) # XNZ
            wildcards.append(kmer[1] + n + kmer[3]) # YNW
            wildcards.append(kmer[:2] + n + kmer[3]) # XYNW
            wildcards.append(kmer[0] + n + kmer[2:]) # XNZW            
            for n2 in N:
                wildcards.append(kmer[0] + n + n2 + kmer[3]) # XN1N2W
    return wildcards


def dinuc_relative_abundance(dinuc_data, kmer_freqs):
    """
    rel_ab <= 0.78 and rel_ab >= 1.23 as suitable benchmarks for 
    assessing whether a dinucleotide is signiﬁcantly over-or 
    under-represented in a DNA sequence.
    (Karlin and Cardon, 1994; Karlin and Ladunga, 1994)
    """
    kmer_rel_ab = defaultdict(float)
    for kmer in dinuc_data.keys():
        b1 = kmer_freqs[kmer[0]]
        b2 = kmer_freqs[kmer[0]]
        rel_ab = dinuc_data[kmer] / (b1 * b2)
        kmer_rel_ab[kmer] = kmer_rel_ab.get(kmer, 0.0) + rel_ab
    return kmer_rel_ab


def diff_dinuc_relative_abundance(spc1_dinuc, spc2_dinuc):
    # two species dinuc_relative_abundance as input 
    vec1 = np.array(list(spc1_dinuc.values()))
    vec2 = np.array(list(spc2_dinuc.values()))
    dif = 1/16*sum(abs(vec1 - vec2))
    return dif


def trinuc_markov_bias(trimer, kmer_counts):
    tri = kmer_counts[trimer]
    mid = kmer_counts[trimer[1]]
    di1 = kmer_counts[trimer[:2]]
    di2 = kmer_counts[trimer[1:]]
    return (tri * mid) / (di1 * di2)


def trinuc_relative_abundance(trinuc_data, kmer_counts):
    """
    rel_ab <= 0.78 and rel_ab >= 1.23 as suitable benchmarks for 
    assessing whether a dinucleotide is signiﬁcantly over-or 
    under-represented in a DNA sequence.
    (Karlin and Cardon, 1994; Karlin and Ladunga, 1994)
    """
    trimer_rel_ab = defaultdict(float)
    for kmer in trinuc_data.keys():
        rel_ab = trinuc_markov_bias(kmer, kmer_counts)
        trimer_rel_ab[kmer] = trimer_rel_ab.get(kmer, 0.0) + rel_ab
    return trimer_rel_ab


def tetra_markov_bias(tetramer, kmer_counts):
    tet = kmer_counts[tetramer]
    mid = kmer_counts[tetramer[1:3]]
    tri1 = kmer_counts[tetramer[:3]]
    tri2 = kmer_counts[tetramer[1:]]
    return (tet * mid) / (tri1 * tri2)


def tetra_relative_abundance(tetranuc_data, kmer_counts):
    """
    rel_ab <= 0.78 and rel_ab >= 1.23 as suitable benchmarks for 
    assessing whether a dinucleotide is signiﬁcantly over-or 
    under-represented in a DNA sequence.
    (Karlin and Cardon, 1994; Karlin and Ladunga, 1994)
    """
    tetramer_rel_ab = defaultdict(float)
    for kmer in tetranuc_data.keys():
        rel_ab = tetra_markov_bias(kmer, kmer_counts)
        tetramer_rel_ab[kmer] = tetramer_rel_ab.get(kmer, 0.0) + rel_ab
    return tetramer_rel_ab


def tmz(zscores, kmer_counts):
    """
    
    """
    tmz = defaultdict(float)
    l = kmer_counts["A"] + kmer_counts["C"] + kmer_counts["G"] + kmer_counts["T"]
    for km, zsc in  zscores.items():
        tz = zsc / math.sqrt(l)
        tmz[km] = tmz.get(km, 0.0) + tz
    return tmz


def select_pal_sum_counts(kmer_counts):
    """
    Select and add up the frequencies of the palindromic kmers,
    found in a DNA sequence.
    """
    # select the palindrome tuples
    # [('TATA', 'ATAT'),..('CTAG', 'GATC')]
    palindromes = []
    for km in kmer_counts:
        kmr = get_reverse_complement(km)
        if km == kmr:
            palindromes.append((km, kmr[::-1]))
    # now it selects the palindromes that are
    # lexically classified
    # {'AATT', 'ACGT', 'AGCT', 'ATAT', ...,'CTAG'}
    # as ATAT cames first than TATA..
    set_pal = set()
    for tup in palindromes:
        if tup[0] < tup[1]:
            set_pal.add(tup[0])
        elif tup[1] < tup[0]:
            set_pal.add(tup[1])
    # Now this part add up the frequencies of the 
    # palindrome pairs Ex, TATA + ATAT freqs
    sum_pal = defaultdict(float)
    for km in set_pal:
        comp = get_strand_complement(km)
        sum_freq = kmer_counts[km] + kmer_counts[comp]
        sum_pal[km] = sum_pal.get(km, 0.0) + sum_freq
    # Returns a new dictionary with counts of the
    # palindromes and the non-palindromes kmers
    # select only the lexically palindromes
    pal = set([item for sublist in palindromes for item in sublist])
    #print(pal)
    # returns a new dictionary with palindromes ad up and
    # all other no palindromic kmers
    new_dict = {k: v for k, v in kmer_counts.items() if k not in pal}
    new_dict.update(sum_pal)
    return new_dict









   
