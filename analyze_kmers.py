#! usr/bin/env python


import os
import sys
import time
import argparse
import glob
from termcolor import colored
import pandas as pd
from system_utils import make_me_a_folder
from alphabet import iupac_dna
from kmer_model_func import get_all_possible_kmers, get_expected_higher_markov, get_kmer_frequency, get_variance, get_standard_deviation, get_z_scores, get_p_values, get_e_values, get_kmer_stats, get_data_from_csv, seq_lengths, get_names, get_scipy_p_values


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

    parser.add_argument("-ki",
                        "--kmin",
                        type=int,
                        dest="kmin",
                        help="Integer representing the minimum length to kmers")

    parser.add_argument("-ka",
                        "--kmax",
                        type=int,
                        dest="kmax",
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
    print(dir_in)
    # name of the sub directory to save the final result
    sub_dir = args.sub_dir # Archaea
    print(sub_dir)
    # Subsubdir where the files lives
    ssub_dir = args.ssub_dir # Asgard
    print(ssub_dir)
    # name of the root directory to save the final result
    dir_out = args.dir_out # Results
    # minimum kmer length
    kmin = args.kmin
    # maximum kmer length
    kmax = args.kmax
    # type of sequence (chr/pls/cds)
    type_seq = args.type_seq
    # alphabet
    alphabet = iupac_dna
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored("The directory to save the files already exists!",
                      "red",
                      attrs=["bold"]))
        pass
    else:
        make_me_a_folder(dir_out)
    # get all sequence lengths
    # Results/GC_slidewindow/Mito/protists/NC_008288.1.fasta.sizes
    #sizes = glob.glob(f"{dir_out}/GC_slidewindow/{ssub_dir}/*.sizes")
    #viruses Others
    file_names = glob.glob(f"Results/Lengths/{sub_dir}/{ssub_dir}/{ssub_dir}_chr_lengths.tsv")
    print(file_names)
    len_seqs = seq_lengths(file_names)
    # Results/GC_slidewindow/Archaea/Asgard/GCA_008000775.1.sizes
    #len_seq = get_seq_length(sizes)
    #print(len_seq)
    # for each separated genomes
    #dir_in = Results/kmer_counts / ssub_dir = Archaea/Asgard _chr.csv
    filenames = glob.glob(f"{dir_in}/{sub_dir}/{ssub_dir}/*_{type_seq}.csv")
    cnt_files = 0
    # input the file paths and print it to show where the script is doing
    # Results/kmer_counts/Archaea/Asgard/GCA_008000775.1_kmer_1_10_cds.csv
    for filename in filenames:
        superkindom, group, name = get_names(filename)
        #to_join = spl_name[0] #, spl_name[1] # [GCA , 008000775.1]
        #name = "_".join(to_join) # GCA_008000775.1
        print(colored(f"Start working with {name}\n", attrs=["bold"]))
        # getting the kmers counts from csvs files
        kmer_count = get_data_from_csv(filename)
        # get the k-mer list for analysis, k = 6
        kmer_list = get_all_possible_kmers(alphabet, kmin, kmax)
        # calculating the expected number for all k-mers
        # the minimum k would be k = 3
        expected = get_expected_higher_markov(kmer_list, kmer_count)
        # gettiing the kmer frequncies
        frequencies = get_kmer_frequency(kmer_list, kmer_count)
        # getting the z-scores
        zscrs = get_z_scores(kmer_list, kmer_count, expected, len_seqs[name])
        # get the p-values from k-mers
        pvals = get_scipy_p_values(zscrs)
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
        df["rank"] = df["z_score"].rank(method="max")
        csv_name = f"{name}_k{kmax}_{type_seq}.model.csv"
        # Results/kmer_model
        full_path = os.path.join("Results", "kmer_model", superkindom, group)
        # Results/kmer_counts/Archaea/Asgard
        complete_path = f"{full_path}"
        if os.path.isfile(complete_path):
            pass
        elif not os.path.exists(full_path):
            os.makedirs(full_path)
        df.to_csv(f"{full_path}/{csv_name}", index=False)
        print(f"Number of kmer (kmin-{kmin}/kmax-{kmax}) from {name}: {len(expected)}\n")
        cnt_files += 1
    # the final time
    end = time.process_time()
    # print some info
    print(colored(f"Total number of files: {cnt_files}\n.",
                  attrs=["bold"]))
    print(colored(f"Total time for the script: {round(end - start, 2)}.",
                  "red",
                  attrs=["bold"]))
    print(colored("Done!",
                  "green",
                  attrs=["bold"]))


if __name__ == "__main__":
    sys.exit(main())

