#!/usr/bin/env python
# -*- coding: utf-8 -*-
# expected z-scores p-values e-values
# usage: kmer_model.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR] [-t TYPE_SEQ]
#                      [-ki KMIN] [-ka KMAX]
# kmer_model.py: error: the following arguments are required: -di/--dir_in
# python kmer_model.py -di Results -sd kmer_counts -do Results -t chr -ki 4 -ka 7 (kmers len 6)
import os
import sys
import time
import argparse
import math
import glob
import csv
from termcolor import colored
from itertools import product
from collections import defaultdict
import pandas as pd
from system_utils import make_me_a_folder
from alphabet import iupac_dna


def get_data_from_csv(filename):
    data = defaultdict(float)
    with open(filename, 'r') as fh:
        csv_data = csv.reader(fh)
        for row in csv_data:
            d1, d2 = row[0], float(row[1])
            data[d1] = data.get(d1, 0.0) + d2
    return data


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
    expected = defaultdict(int)
    for kmer in kmer_list:
        suf, pref, mid = kmer_counts[kmer[1:]], kmer_counts[kmer[:-1]], kmer_counts[kmer[1:-1]]
        if mid == 0:
            expected[kmer] = expected.get(kmer,
                                          0)
        else:
            expected[kmer] = expected.get(kmer,
                                          0) + int((pref * suf) / mid)
    return expected


def get_kmer_frequency(kmer_list, kmer_counts):
    kmer_freq = defaultdict(float, [(kmer, 0.0) for kmer in kmer_list])
    total_freq = 0
    for kmer in kmer_list:
        kmer_freq[kmer] = kmer_counts[kmer]
        total_freq += kmer_counts[kmer]
    return {k: (cnt/total_freq) for k, cnt in kmer_freq.items()}


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
    sigma(W) = sqrt(Expected) * (1 - Expected/len(seq) -k + 1))
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


def get_z_scores(kmer_list, kmer_counts, expected_kmers, len_seq):
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
    for kmer in kmer_list:
        if expected_kmers[kmer] == 0.0:
            z_score[kmer] = 0.0
        else:
            sigma = math.sqrt(expected_kmers[kmer]) * (1 - expected_kmers[kmer] / (2 * len_seq))
            z = (kmer_counts[kmer] - expected_kmers[kmer]) / sigma
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


def get_e_values(kmer_list, p_vals):
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
    hyp_num = len(kmer_list)
    # initialize the container
    e_values = defaultdict(float)
    # iterates through the kmer list
    for kmer in kmer_list:
        # gets the p values from the input container
        p = hyp_num * p_vals[kmer]
        # calculates the e values and add the kmer
        # and the e values to the container
        e_values[kmer] = e_values.get(kmer, 0.0) + p
    return e_values


def get_kmer_stats(kmer_list, 
                   kmer_count, 
                   expected, 
                   frequency,
                   z_scores, 
                   e_vals, 
                   p_vals):
    data = []
    for kmer in kmer_list:
        data.append((kmer,
                     kmer_count[kmer],
                     expected[kmer],
                     frequency[kmer],
                     z_scores[kmer],
                     e_vals[kmer],
                     p_vals[kmer]))
    df = pd.DataFrame(data,
                      columns=['kmer',
                               'observed',
                               'expected',
                               'frequency',
                               'z_score',
                               'e_value',
                               'p_value']).sort_values(by='z_score').reset_index(drop=True)
    return df


def parse_arguments():
    """Parse the command line arguments to kmer_model script.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    The resulting results are csv files from each genome contained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(
        description="""A script to model all kmers of length kmin<=k<=kmax
        from bacterial genomes/plasmids.""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory root. In my case the name is conjugated with a subdir. Ex Results/Bacteria')

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files. Ex. Results/Bacteria')

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Name for a subdirectory, ex., kmer_counts.')

    parser.add_argument('-t',
                        '--type_seq',
                        type=str,
                        dest='type_seq',
                        help='String representing the type of the sequence, ex. chr or plm')

    parser.add_argument('-ki',
                        '--kmin',
                        type=int,
                        dest='kmin',
                        help='Integer representing the minimum length to kmers')

    parser.add_argument('-ka',
                        '--kmax',
                        type=int,
                        dest='kmax',
                        help='Integer representing the maximum length to kmers')

    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time.time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f'\nThe working directory: {cwd}\n',
                  'green',
                  attrs=['bold']))
    # passing the arguments to the script
    args = parse_arguments()
    # name of the input directory, ex. Data/Genomes_splitted
    dir_in = args.dir_in
    # name of the sub directory to save the final result
    # Chromosomes/Plasmids
    sub_dir = args.sub_dir
    # name of the root directory to save the final result
    dir_out = args.dir_out
    # minimum kmer length
    kmin = args.kmin
    # maximum kmer length
    kmax = args.kmax
    # type of sequence (chromosomes/plasmids)
    type_seq = args.type_seq
    # alphabet
    alphabet = iupac_dna
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)
    # get all sequence lengths
    data = dir_in.split("/")[1]
    len_seq = get_data_from_csv(f'{dir_in}/lengths/{data}_lengths.csv')
    # get the csv files path Ex Results/Bacteria/kmer_counts/GCF_.../
    # Results/Archaea/kmer_counts/*/*_kmer_2_8_chr.csv
    filenames = glob.glob(f'{dir_in}/{sub_dir}/*/*_kmer_2_8_{type_seq}.csv')
    # initialize the file counter
    cnt_files = 0
    # input the file paths and print it to show where the script is doing
    # 'Results/Mito/kmer_counts/NC_041285.1_Chelodina_burrungandjii_mito/NC_041285.1_Chelodina_burrungandjii_mito_kmer_2_8_mito.csv'
    for filename in filenames:
        # NC_041285.1_Chelodina_burrungandjii_mito
        filename_splited = filename.split('/')
        name = filename_splited[3]
        print(colored(f"Start working with {name}\n", attrs=['bold']))
        # getting the kmers counts from csvs files
        kmer_count = get_data_from_csv(filename)
        # get the k-mer list for analysis, k = 6
        kmer_list = get_all_possible_kmers(alphabet, kmin, kmax)
        # calculating the expected number for all k-mers
        expected = get_expected_higher_markov(kmer_list, kmer_count)
        # gettiing the kmer frequncies
        frequencies = get_kmer_frequency(kmer_list, kmer_count)
        # getting the z-scores
        zscrs = get_z_scores(kmer_list, kmer_count, expected, len_seq[name])
        # get the p-values from k-mers
        pvals = get_p_values(zscrs)
        # get the k-mers e-values
        evals = get_e_values(kmer_list, pvals)
        # saving the final results as a csv file
        df = get_kmer_stats(kmer_list, 
                            kmer_count, 
                            expected, 
                            frequencies, 
                            zscrs, 
                            evals, 
                            pvals)
        df['rank'] = df['z_score'].rank(method='max')
        csv_name = f'{name}_k{kmax}_{type_seq}.model.csv'
        # Results/kmer_model/Mito/NC_041285.1_Chelodina_burrungandjii_mito
        sub_dir_out = dir_in.split("/")[1]
        full_path = os.path.join('Results', f'{sub_dir_out}','kmer_model', name)
        complete_path = f'{full_path}/{csv_name}'
        if os.path.isfile(complete_path):
            pass
        elif not os.path.exists(full_path):
            os.makedirs(full_path)
        df.to_csv(f'{full_path}/{csv_name}', index=False)
        print(f'Number of kmer (kmin-{kmin}/kmax-{kmax}) from {name}: {len(expected)}\n')
        # add to the count file
        cnt_files += 1
    # the final time
    end = time.time()
    # print some info
    print(colored(f"Total number of files: {cnt_files}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
