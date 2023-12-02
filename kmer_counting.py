#!/usr/bin/env python
# -*- coding: utf-8 -*-
# usage: kmer_counting.py [-h] -di Path [-sd path] [-do path] -ki KMIN -ka KMAX
#                         [-t SEQ_TYPE]
# kmer_counting.py: error: the following arguments are required: -di/--dir_in, -ki/--kmin, -ka/--kmax
#python scr/kmer_counting.py -di Genomes/Bacteria/Acidobacteria/Chr -do Results/kmer_counts -ki 1 -ka 10 -e fna -t chr

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
sys.path.append("/media/paulosschlogl/Paulo/pauloscchlogl/Genome_kmers/kmer_counter")
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
                        help="Sequence type.(Chr/plas)")
                        
    parser.add_argument("-r",
                        "--rev_seq",
                        type=bool,
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
    # sequence type chromosome/plasmid/cds
    rev = args.rev    
    # folder to save the results
    dir_out = args.dir_out
    supergroup = dir_in.split("/")[1] # Bacteria
    group = sub_dir.split("/")[0] # Acidobacteria
    subgroup = sub_dir.split("/")[1] # Acidobacteriia
    full_path = os.path.join(dir_out, "kmer_counts", supergroup, group, subgroup, ssub_dir)
    if not os.path.exists(full_path):
        os.makedirs(full_path)
    # Genomes/Bacteria Acidobacteria/Acidobacteriia CHR/CDS/PLAS # _chr.fna.gz
    # _cds.fna _pls.fna.gz
    filenames = glob.glob(f"{dir_in}/{sub_dir}/{ssub_dir}/*_{extension}")
    #print(filenames)
    for filename in filenames:
        my_path = Path(filename)
        name = my_path.stem
        name = name.replace(f"_{extension}", "")
        myfilename = f"{name}_kmer_{kmin}_{kmax}_{seq_type}.csv"
        if Path(f"{full_path}/{myfilename}").is_file():
            print(colored(f"{myfilename} already exists! Done counting!!!!!!!\n", 
                          "red", attrs=["bold"]))
            pass
        for Id, sequence in fasta_parser.parse_fasta(filename):
            print(colored(f"Working with genome {name}\n{filename}\n",
            "green"))
            print(colored(f"\nCounting kmer of length {kmin}-{kmax}", "green", attrs=["bold"]))
            if rev:
                kmer_count_rev = mer_count(get_reverse_complement(sequence), kmin, kmax)
                dfr = pd.DataFrame(kmer_count_rev_clean.items(), columns=["kmer", "observed_rev"])
                df = pd.merge(dff, dfr, on="kmer")
            else:
                kmer_count = mer_count(sequence, kmin, kmax)
                dff = pd.DataFrame(kmer_count_clean.items(), columns=["kmer", "observed"])
                df.to_csv(f"{full_path}/{myfilename}", header=False, index=False)
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
