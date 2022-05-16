#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage: kmer_counting.py [-h] -di Path [-sd path] [-do path] -ki KMIN -ka KMAX
#                         [-t SEQ_TYPE]
# kmer_counting.py: error: the following arguments are required: -di/--dir_in, -ki/--kmin, -ka/--kmax
import argparse
import glob
import os
import sys
import time
from termcolor import colored
from collections import defaultdict
import pandas as pd
import fasta_parser
import count_kmers


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description='kmer count')
    parser.add_argument('-di',
                        '--dir_in',
                        metavar='Path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory to the input files.(Data/name)')

    parser.add_argument('-sd',
                        '--sub_dir',
                        metavar='path',
                        type=str,
                        required=False,
                        dest='sub_dir',
                        help='Subdirectory of the input files.(Chrmosomes)')

    parser.add_argument('-do',
                        '--dir_out',
                        metavar='path',
                        type=str,
                        required=False,
                        dest='dir_out',
                        help='Directory to the ouput files.(Results)')

    parser.add_argument('-ki',
                        '--kmin',
                        type=int,
                        required=True,
                        dest='kmin',
                        help='Integer representing the minimum length of the kmer substring.')

    parser.add_argument('-ka',
                        '--kmax',
                        type=int,
                        required=True,
                        dest='kmax',
                        help='Integer representing the maximum length of the kmer substring.')

    parser.add_argument('-e',
                        '--extension',
                        type=str,
                        required=False,
                        dest='extension',
                        help='File type.(fna.gz/fna)')

    parser.add_argument('-t',
                        '--seq_type',
                        type=str,
                        required=False,
                        dest='seq_type',
                        help='Sequence type.(Chromosomes/plasmids)')

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
    dir_in = args.dir_in
    sub_dir = args.sub_dir
    kmin = args.kmin
    kmax = args.kmax
    extension = args.extension
    seq_type = args.seq_type
    dir_out = args.dir_out
    filenames = []
    if sub_dir != None:
        # if there are a subdirectory (Chrmosome or Plasmid)
        filenames += glob.glob(f"{dir_in}/*/{sub_dir}/*.{extension}")
    else:
        filenames += glob.glob(f"{dir_in}/*/*.{extension}")
    seq_lengths = defaultdict(int)
    for filename in filenames:
        for Id, sequence in fasta_parser.parse_fasta(filename):
            name = filename.split('/')[-1][:-4]
            print(colored(f"Working with genome {name}\n{filename}\n",
            'green'))
            seq_len = len(sequence)
            seq_lengths[(name, Id)] = seq_len
            print(colored(f"counting kmer of length {kmin}-{kmax}", attrs=['bold']))
            kmer_count = count_kmers.count_kmers(sequence, kmin, kmax)
            kmer_count_clean = {k: cnt for k, cnt in kmer_count.items() if set(k).issubset(set("ACGT"))}
            df = pd.DataFrame(kmer_count_clean.items(), columns=['kmer', 'observed'])
            full_path = os.path.join(dir_out, 'kmer_counts', name)
            file_name = f'{name}_kmer_{kmin}_{kmax}_{seq_type}.csv'
            if not os.path.exists(full_path):
                os.makedirs(full_path)
            df.to_csv(f'{full_path}/{file_name}', header=True, index=False)
    df_l = pd.DataFrame(seq_lengths.items(), columns=['name', 'length'])
    path = os.path.join(dir_out, 'lengths')
    outfile = f'lenghts_{seq_type}.csv'
    if not os.path.exists(path):
        os.makedirs(path)
    df_l.to_csv(f'{path}/{outfile}', header=False, index=False)
    # the final time
    end = time.time()
    # print some info
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
