#!/usr/bin/env python

#usage: GC_slide_window.py [-h] -di Path [-s SEQTYPE] [-do DIR_OUT] -w WINDOW

#Calculates the over/under represented kmers of length k from genomes using higher order markov models through
#R'mes software.

#options:
#  -h, --help            show this help message and exit
#  -di Path, --dir_in Path
#                        Directory to the input files.Ex. (Genomes_new/Archaea/C_Heimdallarchaeota)
#  -s SEQTYPE, --seqtype SEQTYPE
#                        Subdirectory of the input files.(PLAS/CHR)
#  -do DIR_OUT, --dir_out DIR_OUT
#                        Directory to the ouput files.(Results)
#  -w WINDOW, --window WINDOW
#                        Integer representing the length of the kmer substring.

import argparse
import glob
import os
import sys
import csv
import re
import time
from pathlib import Path
from termcolor import colored
from fasta_parser import get_ID_number, get_files


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. 
    These options can be accessed using args.option. 
    For example args.Path stores the Path provided.
    """
    parser = argparse.ArgumentParser(description="""
                                     Calculates the GC contente in slide windows
                                     of length window from genomes.
                                     """)
                                     
    parser.add_argument("-di",
                        "--dir_in",
                        metavar="Path",
                        type=str,
                        required=True,
                        dest="dir_in",
                        help="""Directory to the input files.Ex. 
                        (Genomes_new/Archaea/C_Heimdallarchaeota)""")

    parser.add_argument("-sd",
                        "--sub_dir",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="sub_dir",
                        help="Subdirectory of the input files.(CHR/PLAS)")

    parser.add_argument("-do",
                        "--dir_out",
                        type=str,
                        required=False,
                        dest="dir_out",
                        help="Directory to the ouput files.(Results)")

    parser.add_argument("-w",
                        "--window",
                        type=int,
                        required=True,
                        dest="window",
                        help="Integer representing the length of the kmer substring.")

    parser.add_argument("-e",
                        "--extension",
                        type=str,
                        required=False,
                        dest="extension",
                        help="File type.(csv/tsv)")

    parser.add_argument("-p",
                        "--pattern",
                        type=str,
                        required=False,
                        dest="pattern",
                        help="Pattern representing the strands of the input sequences.(fwr/rev)")
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
    # kmer length
    window = args.window
    # Chr/CDS/Plas
    pattern = args.pattern
    # file extention .fna.gz etc
    extension = args.extension
    # regex to get the assembly number
    rg =  re.compile(r"(GC[AF]_\d+\.\d)")
    # folder to save the results
    dir_out = args.dir_out
    # Genomes_new/Archaea/Asgard_group/C_Heimdallarchaeota
    superg = dir_in.split("/")[1] # Archaea
    gr = dir_in.split("/")[2] # Asgard_group
    subgr = dir_in.split("/")[3] # C_Heimdallarchaeota
    # Genomes_new/Archaea/Asgard_group/C_Heimdallarchaeota/GCA_020351745.1/CHR/CP070665.1.fna
    print(f"{dir_in}/*/{sub_dir}/*.{extension}")
    filenames = glob.glob(f"{dir_in}/*/{sub_dir}/*.{extension}")
    #print(filenames)
    cnt = 0
    for filename in filenames:
        # GCA_020351745.1
        infolder = os.path.dirname(filename)
        geno_id = get_ID_number(filename, rg)
        my_path = Path(filename)
        name = my_path.stem
        full_path = os.path.join(dir_out, 
                                 superg, 
                                 gr, 
                                 subgr,
                                 geno_id,
                                 sub_dir)
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        print(f"Calculating the GC content from {geno_id} file {name} in a window of {window}\n")
        os.system(f"samtools faidx {filename}")
        os.system(f"cut -f 1,2 {filename}.fai > {infolder}/{name}.sizes")
        os.system(f"bedtools makewindows -g {infolder}/{name}.sizes -w {window} > {infolder}/{name}.{window}.bps.bed")
        os.system(f"bedtools nuc -fi {filename} -bed {infolder}/{name}.{window}.bps.bed > {full_path}/{name}.{window}_gc_nucs.txt")
        cnt += 1
    # the final time
    end = time.time()
    # print some info
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




