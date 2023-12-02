#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage: kmer_counting.py [-h] -di Path [-sd path] [-ssd path] [-do path] -ki
#                        KMIN -ka KMAX [-e EXTENSION] [-t SEQ_TYPE] [-r REV]
#
#kmer count
#options:
#  -h, --help            show this help message and exit
#  -di Path, --dir_in Path
#                        Directory to the input files.(Genomes/Bacteria)
#  -sd path, --sub_dir path
#                        Subdirectory of the input
#                        files.(Acidobacteria/Acidobacteriia)
#  -ssd path, --ssub_dir path
#                        Subdirectory of the input files.(CHR/CDS)
#  -do path, --dir_out path
#                        Directory to the ouput files.(Results)
#  -ki KMIN, --kmin KMIN
#                        Integer representing the minimum length of the kmer
#                        substring.
#  -ka KMAX, --kmax KMAX
#                        Integer representing the maximum length of the kmer
#                        substring.
#  -e EXTENSION, --extension EXTENSION
#                        File type.(fna.gz/fna)
#  -t SEQ_TYPE, --seq_type SEQ_TYPE
#                        Sequence type.(Chr/plas)
#  -r REV, --rev_seq REV
#                        If reversed sequence counting needed.
#python kmer_counting.py -di Genomes/Bacteria/Acidobacteria/Chr -do Results/kmer_counts -ki 1 -ka 10 -e fna -t chr

import argparse
import glob
import os
import sys
import time
from pathlib import Path
from termcolor import colored
from collections import defaultdict
import pandas as pd
import fasta_parser
from kmer_counter.count_kmers import mer_count
from sequence_utils import get_reverse_complement


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description="kmer count")
    parser.add_argument("-di",
                        "--dir_in",
                        metavar="Path",
                        type=str,
                        required=True,
                        dest="dir_in",
                        help="Directory to the input files.(Genomes/Bacteria)")

    parser.add_argument("-sd",
                        "--sub_dir",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="sub_dir",
                        help="Subdirectory of the input files.(Acidobacteria/Acidobacteriia)")

    parser.add_argument("-ssd",
                        "--ssub_dir",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="ssub_dir",
                        help="Subdirectory of the input files.(CHR/CDS)")

    parser.add_argument("-do",
                        "--dir_out",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="dir_out",
                        help="Directory to the ouput files.(Results)")

    parser.add_argument("-ki",
                        "--kmin",
                        type=int,
                        required=True,
                        dest="kmin",
                        help="Integer representing the minimum length of the kmer substring.")

    parser.add_argument("-ka",
                        "--kmax",
                        type=int,
                        required=True,
                        dest="kmax",
                        help="Integer representing the maximum length of the kmer substring.")

    parser.add_argument("-e",
                        "--extension",
                        type=str,
                        required=False,
                        dest="extension",
                        help="File type.(fna.gz/fna)")

    parser.add_argument("-t",
                        "--seq_type",
                        type=str,
                        required=False,
                        dest="seq_type",
                        help="Sequence type.(chr/plas)")
                 
    parser.add_argument("-r",
                        "--rev_seq",
                        action='store_true',
                        required=False,
                        dest="rev",
                        help="If reversed sequence counting needed.")
    
    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time.time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f"\nThe working directory: {cwd}\n",
                  "red",
                  attrs=["bold"]))
    # passing the arguments to the script
    args = parse_arguments()
    # Genomes/Bacteria
    dir_in = args.dir_in
    # Acidobacteria/Acidobacteriia
    sub_dir = args.sub_dir
    # Chr/CDS/Plas
    ssub_dir = args.ssub_dir
    # minimum kmer length
    kmin = args.kmin
    # maximum kmer length
    kmax = args.kmax
    # file extention .fna.gz etc
    extension = args.extension
    # sequence type chromosome/plasmid/cds
    seq_type = args.seq_type
    # count only reversed sequence type chromosome/plasmid/cds
    rev = args.rev    
    # folder to save the results
    dir_out = args.dir_out
    supergroup = dir_in.split("/")[1] # Bacteria
    group = sub_dir.split("/")[0] # Acidobacteria
    subgroup = sub_dir.split("/")[1] # Acidobacteriia
    full_path = os.path.join(dir_out, "kmer_counts", supergroup, group, subgroup, ssub_dir)
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    # Genomes_new/Bacteria Pseudomonadota/Acidithiobacillia CHR/ chr fna
    filenames = glob.glob(f"{dir_in}/{sub_dir}/{ssub_dir}/*_{seq_type}.{extension}")
    #print(filenames)
    for filename in filenames:
        my_path = Path(filename)
        name = my_path.stem
        name = name.replace(f"_{extension}", "")
        myfilename = f"{name}_kmer_{kmin}_{kmax}_{seq_type}.csv"
        if Path(f"{full_path}/{myfilename}").is_file() and os.stat(f"{full_path}/{myfilename}").st_size != 0:
            print(colored(f"{myfilename} already exists! Done counting!!!!!!!\n", 
                          "red", attrs=["bold"]))
            pass
        for Id, sequence in fasta_parser.parse_fasta(filename):
            print(colored(f"Working with genome {name}\n",
            "green"))
            print(colored(f"\nCounting kmer of length {kmin}-{kmax}", "green", attrs=["bold"]))
            df = pd.DataFrame()
            if rev:
                seq = sequence + get_reverse_complement(sequence)
                kmer_count_rev = mer_count(seq, kmin, kmax)
                dfr = pd.DataFrame(kmer_count_rev.items(), columns=["kmer", "observed"])
                df = pd.concat([df, dfr])
            else:
                kmer_count = mer_count(sequence, kmin, kmax)
                dff = pd.DataFrame(kmer_count.items(), columns=["kmer", "observed"])
                df = dff
            df.to_csv(f"{full_path}/{myfilename}", header=True, index=False)
    # the final time
    end = time.time()
    # print some info
    print(colored(f"Total time for the script finishes: {round(end - start, 2)}.",
                  "red",
                  attrs=["bold"]))
    print(colored("Done!",
                  "green",
                  attrs=["bold"]))


if __name__ == "__main__":
    sys.exit(main())
