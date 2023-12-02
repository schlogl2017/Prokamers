#!/usr/bin/env python
# -*- coding: utf-8 -*-
# expected z-scores p-values e-values
#usage: kmer_model.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR] [-ssd SSUB_DIR]
#                     [-t TYPE_SEQ] [-ki KMIN] [-ka KMAX]
#kmer_model.py: error: the following arguments are required: -di/--dir_in
# python scr/kmer_model.py -di Results/kmer_counts -sd kmer_model -ssd Archaea/Thermoplasmatota -do Results -t chr -ki 3 -ka 3
# python scr/kmer_model.py -di Results/kmer_counts -do Results -sd Archaea -ssd Euryarchaeota -t # chr -ki 7 -ka 7

import re
import os
import sys
import time
import argparse
import glob
from termcolor import colored
import pandas as pd
from system_utils import make_me_a_folder
from alphabet import iupac_dna
from kmer_model_funcs import *


def parse_arguments():
    """Parse the command line arguments to kmer_model script.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    The resulting results are csv files from each genome contained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(
        description="""A script to model all kmers of length kmin <= k <= kmax
        from bacterial genomes/plasmids.
        Usage:python scr/kmer_model.py -di Results/kmer_counts -sd kmer_model -ssd Archaea/Thermoplasmatota -do Results -t chr -ki 3 -ka 3
        """,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-di",
                        "--dir_in",
                        metavar="path",
                        type=str,
                        required=True,
                        dest="dir_in",
                        help="Directory root. Ex Results/kmer_counts")

    parser.add_argument("-do",
                        "--dir_out",
                        type=str,
                        dest="dir_out",
                        help="directory name for resulting files. Ex. Results")

    parser.add_argument("-sd",
                        "--sub_dir",
                        type=str,
                        dest="sub_dir",
                        help="Name for a subdirectory, ex., Archaea.")

    parser.add_argument("-ssd",
                        "--ssub_dir",
                        type=str,
                        dest="ssub_dir",
                        help="Name for a subdirectory, ex., Asgard.")

    parser.add_argument("-t",
                        "--type_seq",
                        type=str,
                        dest="type_seq",
                        help="String representing the type of the sequence, ex. chr/plm/cds")

    parser.add_argument("-ka",
                        "--kmax",
                        type=int,
                        dest="kmax",
                        help="Integer representing the maximum length to kmers")

    parser.add_argument("-ki",
                        "--kmin",
                        type=int,
                        dest="kmin",
                        help="Integer representing the maximum length to kmers")

    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time.process_time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f"\nThe working directory: {cwd}\n",
                  "green",
                  attrs=["bold"]))
    # passing the arguments to the script
    args = parse_arguments()
    # name of the input directory
    dir_in = args.dir_in # Results/kmer_counts
    # name of the sub directory to save the final result
    sub_dir = args.sub_dir # Archaea
    # Subsubdir where the files lives
    ssub_dir = args.ssub_dir # Asgard
    # name of the root directory to save the final result
    dir_out = args.dir_out # Results
    # maximum kmer length
    kmax = args.kmax
    # maximum kmer length
    kmin = args.kmin    
    # type of sequence (chr/pls/cds)
    type_seq = args.type_seq
    # alphabet
    alphabet = iupac_dna
    rg =  re.compile(r"(GC[AF]_\d+\.\d)")
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored("The directory to save the files already exists!",
                      "red",
                      attrs=["bold"]))
        pass
    else:
        make_me_a_folder(dir_out)
    # -di Results/kmer_counts -do Results -sd Bacteria -ssd Acidithiobacillia -t chr -ka 7
    filenames = glob.glob(f"{dir_in}/{sub_dir}/{ssub_dir}/CHR/*_{type_seq}.csv")
    print(f"{dir_in}/{sub_dir}/{ssub_dir}/CHR/*_{type_seq}.csv")
    cnt_files = 0
    for filename in filenames:
        superkindom, group, name = get_names(filename)
        name = rg.search(name).group()
        print(colored(f"Start working with {name}\n", attrs=["bold"]))
        # read the kmer count data in
        kmer_counts = get_data_from_csv(filename)
        # select the kmer data according with the kmer length_sequence
        kmer_data = group_by_length(kmer_counts, kmax)
        # calculate the expected kmer values and it variance
        kmer_exp, kmer_var = get_expected_higher_markov(kmer_data, kmer_counts)
        # calculate the standard deviation of the expected kmer values
        kmer_std = get_standard_deviation(kmer_var)
        # Calculates the z-scores
        kmer_zscores = get_z_scores(kmer_data, kmer_exp, kmer_std)
        # Calculates the p-values
        kmer_pvalues = get_scipy_p_values(kmer_zscores)
        # Calculates the evalues
        kmer_evalues = get_e_values(kmer_pvalues)
        # Calculates the observed kmer frequency
        kmer_freq = get_kmer_frequency(kmer_data)
        # Calculates the expected kmer frequency
        kmer_exp_freq = get_kmer_frequency(kmer_exp)
        # make a dataframe
        kmer_df = get_kmer_stats(kmer_data, kmer_freq, kmer_exp, kmer_exp_freq, kmer_zscores, kmer_pvalues, kmer_evalues)
        csv_name = f"{name}_k{kmax}_{type_seq}.model.csv"
        # Results/kmer_model
        full_path = os.path.join(f"{dir_out}", "kmer_model", superkindom, group)
        # Results/kmer_counts/Archaea/Asgard
        complete_path = f"{full_path}"
        if os.path.isfile(complete_path):
            pass
        elif not os.path.exists(full_path):
            os.makedirs(full_path)
        kmer_df.to_csv(f"{full_path}/{csv_name}", index=False)
        print(f"Number of kmer (kmin-{kmin}/kmax-{kmax}) from {name}: {len(kmer_data)}\n")
        cnt_files += 1
    # the final time
    end = time.process_time()
    # print some info
    print(colored(f"Total number of files: {cnt_files}.\n",
                  attrs=["bold"]))
    print(colored(f"Total time for the script: {round(end - start, 2)}.",
                  "red",
                  attrs=["bold"]))
    print(colored("Done!",
                  "green",
                  attrs=["bold"]))


if __name__ == "__main__":
    sys.exit(main())

