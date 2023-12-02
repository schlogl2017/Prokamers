#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Script to calculate the expected number of kmers by a High Markov Order
# Author: Paulo Sérgio Schlögl
# Data: 15/05/2023
# Last modification: 29/08/2023
# copyrights 

#usage: expected_kmer_by_hom.py [-h] -di Path [-sd path] [-p PATTERN]
#                               [-do path] -k K_LEN [-e EXTENSION]
#Calculates the expected kmer values by a higher order markov model.
#options:
#  -h, --help            show this help message and exit
#  -di Path, --dir_in Path
#                        Directory to the input files.(Jelly/Bacteria/Pseudomon
#                        adota/Acidithiobacillia)
#  -sd path, --sub_dir path
#                        Subdirectory of the input files.(CHR/PLAS)
#  -p PATTERN, --pattern PATTERN
#                        Pattern representing the strands of the input
#                        sequences.(fwr/rev)
#  -do path, --dir_out path
#                        Directory to the ouput files.(Results)
#  -k K_LEN, --k_len K_LEN
#                        Integer representing the length of the kmer substring.
#  -e EXTENSION, --extension EXTENSION
#                        File type.(csv/tsv)


import argparse
import glob
import os
import sys
import csv
import re
import time
from pathlib import Path
from termcolor import colored
from collections import defaultdict
import pandas as pd
from kmer_model_funcs import get_expected_higher_markov, calculate_expected_kmers_and_variance, group_by_length, get_assembly_ids
from alphabet import iupac_dna


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. 
    These options can be accessed using args.option. 
    For example args.Path stores the Path provided.
    """
    parser = argparse.ArgumentParser(description="""
                                     Calculates the expected kmer values
                                     by a higher order markov model.

                                     """)
    parser.add_argument("-di",
                        "--dir_in",
                        metavar="Path",
                        type=str,
                        required=True,
                        dest="dir_in",
                        help="Directory to the input files.(Jelly/Bacteria/Pseudomonadota/Acidithiobacillia)")

    parser.add_argument("-sd",
                        "--sub_dir",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="sub_dir",
                        help="Subdirectory of the input files.(CHR/PLAS)")

    parser.add_argument("-p",
                        "--pattern",
                        type=str,
                        required=False,
                        dest="pattern",
                        help="Pattern representing the strands of the input sequences.(fwr/rev)")

    parser.add_argument("-do",
                        "--dir_out",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="dir_out",
                        help="Directory to the ouput files.(Results)")

    parser.add_argument("-k",
                        "--k_len",
                        type=int,
                        required=True,
                        dest="k_len",
                        help="Integer representing the length of the kmer substring.")

    parser.add_argument("-e",
                        "--extension",
                        type=str,
                        required=False,
                        dest="extension",
                        help="File type.(csv/tsv)")
    
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
    # CHR/PLAS
    sub_dir = args.sub_dir
    # Chr/CDS/Plas
    pattern = args.pattern
    # kmer length
    k_len = args.k_len
    # file extention .fna.gz etc
    extension = args.extension
    # regex to get the assembly number
    rg_geno =  re.compile(r"(GC[AF]_\d+\.\d)")
    rg_chr_name = re.compile(r"\w+\d+\.\d+_")
    # alphabet
    alphabet = iupac_dna
    # folder to save the results
    #  Archaea Unclas_Archaeon Unclas_Archaeon
    dir_out = args.dir_out
    supergr = dir_in.split("/")[1] # Bacteria/Archaea
    gr = dir_in.split("/")[2] # Pseudomonadota/Unclas_Archaeon
    subgr = dir_in.split("/")[3] # Acidithiobacillia
    # Jelly_counts/Archaea/Unclas_Archaeon/Unclas_Archaeon/*/CHR/CP010426.1_fwr.csv
    filenames = glob.glob(f"{dir_in}/*/{sub_dir}/*_{pattern}.{extension}")
    #print(filenames)
    #num_files = len(filenames)
    cnt = 0
    for filename in filenames:
        #my_path = Path(filename)
        name = re.search(rg_chr_name, filename).group(0).replace("_", "")
        geno_id = get_assembly_ids(filename, rg_geno)
        myfilename = f"{name}_kmer_{k_len}_expec_{pattern}.csv"
        full_path = os.path.join(dir_out, # HMO_exp
                                 supergr, # Archaea
                                 gr, # Unclas_Archaeon
                                 subgr, # Unclas_Archaeon
                                 geno_id, #GCA_000830275.1
                                 sub_dir # CHR/PLAS
                                 )
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        if Path(f"{full_path}/{myfilename}").is_file() and os.stat(f"{full_path}/{myfilename}").st_size != 0:
            print(colored(f"{myfilename} already exists! Done calculating!!!!!!!\n", 
                          "red", attrs=["bold"]))
            pass
        print(f"Calculating the expected kmer (k = {k_len}) values from {name}\n")
        df = calculate_expected_kmers_and_variance(filename, k_len, alphabet, True)
        df.to_csv(f"{full_path}/{myfilename}", header=False, index=False)
        cnt += 1
    # the final time
    end = time.time()
    # print some info
#    print(colored(f"Total number of files: {num_files}", 
#                  "green", 
#                  attrs=["bold"]))
    print(colored(f"Number of files manipulated: {cnt}", 
                  "green", 
                  attrs=["bold"]))
    print(colored(f"Total time for the script finishes: {round(end - start, 2)}.",
                  "green",
                  attrs=["bold"]))
    print(colored("Done!",
                  "red",
                  attrs=["bold"]))


if __name__ == "__main__":
    sys.exit(main())


