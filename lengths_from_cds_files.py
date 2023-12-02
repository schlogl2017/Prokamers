#! usr/bin/env python

import argparse
import glob
import os
import sys
import time
from termcolor import colored
from collections import defaultdict
import pandas as pd
import fasta_parser


def get_cds_lens(filename):
    if "Bacteria" or "Archaea" or "Viruses" in filename:
        name = "_".join(filename.split("/")[-1].split("_")[:2])
    else:
        name = filename.split("/")[-1].split("_")[0]
    for Id, sequence in fasta_parser.parse_fasta(filename):
        gene = Id.split(" ")[0]
        seq_len = len(sequence)
        print(f"{name}\t{gene}\t{seq_len}")


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. For example opt.alp stores the alphabet provided.
    """
    parser = argparse.ArgumentParser(description="Calculate the length of a multi fasta file.")
    # Genomes/Bacteria/
    parser.add_argument("-di",
                        "--dir_in",
                        metavar="Path",
                        type=str,
                        required=True,
                        dest="dir_in",
                        help="Directory to the input files.")
    # Acidobacteria/Acidobacteriia
    parser.add_argument("-sd",
                        "--sub_dir",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="sub_dir",
                        help="Subdirectory of the input files.(Chrmosomes)")
    # Resuls/Lengths
    parser.add_argument("-do",
                        "--dir_out",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="dir_out",
                        help="Directory to the ouput files.(Results)")
    # *_cds_from_genomic.fna.gz / _CDS.fasta
    parser.add_argument("-e",
                        "--extension",
                        type=str,
                        required=False,
                        dest="extension",
                        help="File type.(_cds_from_genomic.fna.gz/_CDS.fasta.gz)")
    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time.time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f"\nThe working directory: {cwd}\n",
                  "grey",
                  attrs=["bold"]))
    # passing the arguments to the script
    args = parse_arguments()
    dir_in = args.dir_in
    sub_dir = args.sub_dir
    extension = args.extension
    dir_out = args.dir_out
    supergroup = dir_in.split("/")[1] # Genomes/Bacteria/
    group = sub_dir.split("/")[0] # Acidobacteria/Acidobacteriia
    sub_group = sub_dir.split("/")[1] # Acidobacteria/Acidobacteriia
    # Genomes/Mito/Animals/Amphibians/CDS/KX686108.1_CDS.fasta.gz
    # Genomes/Bacteria/Acidobacteria/Acidobacteriia/CDS/GCA_000014005.1_ASM1400v1_cds_from_genomic.fna.gz
    filenames = glob.glob(f"{dir_in}/{sub_dir}/CDS/*_{extension}")

    cds_lengths = defaultdict(dict)
    for filename in filenames:
        if "Bacteria" or "Archaea" or "Viruses" in filename:
            name = "_".join(filename.split("/")[-1].split("_")[:2])
        else:
            name = filename.split("/")[-1].split("_")[0]
            cds_lengths[name] = cds_lengths.get(name, {})
        for Id, sequence in fasta_parser.parse_fasta(filename):
            gene = Id.split(" ")[0]
            cds_lengths[name][gene] = len(sequence)
    #df = pd.DataFrame.from_dict(cds_lengths)
    df = pd.DataFrame(cds_lengths)
    print(df.T)
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
