#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 19:47:38 2023

@author: paulo
"""
import os
import argparse
import glob
import sys
import re
import time
import math
from pathlib import Path
from termcolor import colored
from collections import defaultdict
import pylab
from matplotlib import pyplot as plt
from matplotlib import cm
#from kmer_model_funcs import get_assembly_ids
from data_to_plot_R import get_counts_files, get_lengths_files_csv, oligo_count_dict


def get_len(path_root):
    en_dict = defaultdict(int)
    for filename in get_lengths_files_csv(path_root, 
                                          pat="_lengths.csv"):
        with open(filename) as csvfile:
            filereader = csv.reader(csvfile)
            for row in filereader:
                ID = row[0]
                length = int(row[1])
                len_dict[ID] = length
    return len_dict


def probabilities(kmer_counts, seq_len, k):
    probs = defaultdict(float)
    N = (seq_len - k + 1)
    for km, cnt in kmer_counts.items():
        probs[km] = cnt / N
    return probs


def chaos_game_representation(probabilities, k):
    size = int(math.sqrt(4**k))
    chaos = [[0] * size for _ in range(size)]
    max_x = size
    max_y = size
    pos_x = 1
    pos_y = 1
    for km, value in probabilities.items():
        for char in key:
            if c == "T":
                pos_x += max_x / 2
            elif c == "C":
                pos_y += max_y / 2
            elif c == "G":
                 pos_x += max_x / 2
                 pos_y += max_y / 2
            max_x = max_x / 2
            max_y /= 2
        chaos[pos_y-1][pos_x-1] = value
        max_x = size
        max_y = size
        pos_x = 1
        pos_y = 1
    return chaos


def parse_arguments():
    """Parse the command line arguments to the the script.
    Sets up the argparse command-line parser and calls it. 
    These options can be accessed using args.option. 
    For example args.Path stores the Path provided.
    """
    parser = argparse.ArgumentParser(description="""
                                     Create the data to plot kmer
                                     distribution on prokaryotic genomes.

                                     """)
    parser.add_argument("-d1",
                        "--dir1",
                        metavar="Path",
                        type=str,
                        required=True,
                        dest="dir1",
                        help="Directory to the count input file. Ex Genomes_new")

    parser.add_argument("-d2",
                        "--dir2",
                        metavar="Path",
                        type=str,
                        required=True,
                        dest="dir2",
                        help="Directory to the input lenght file. Ex Counts")

    parser.add_argument("-p",
                        "--pattern",
                        type=str,
                        required=False,
                        dest="pattern",
                        help="Pattern of the input file.")

    parser.add_argument("-do",
                        "--dir_out",
                        metavar="Path",
                        type=str,
                        required=False,
                        dest="dir_out",
                        help="Directory to the ouput files. Ex Oligo_plots/Data")

    parser.add_argument("-k",
                        "--k_len",
                        type=int,
                        required=True,
                        dest="k_len",
                        help="Integer representing the length of the kmer substring.")
    
    parser.add_argument("-t",
                        "--seqtype",
                        type=str,
                        required=False,
                        dest="seqtype",
                        help="String representing th type of sequence. Ex PLAS/CHR")

    parser.add_argument("-n",
                        "--names",
                        type=str,
                        required=True,
                        dest="names",
                        help="String representing th type of sequence. Ex Archaea and Bacteria")
                        
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
    # Lengths directory files
    dir1 = args.dir1
    # Acidobacteria/Acidobacteriia
    dir2 = args.dir2
    # pattern of the file
    pattern = args.pattern
    # kmer length
    k_len = args.k_len
    # sequence type CHR or PLAS
    seqtype = args.seqtype
    # names Bacteria Archaea
    names = args.names
    # regex to get the assembly number
    rg =  re.compile(r"(GC[AF]_\d+\.\d)")
    # length files directory
    dir_count = f"{dir1}/{a}"
    dir_len = f"{dir1}/{b}"
    
    if pattern != None:
        data = oligo_count_dict(dir_count, k_len, rg, seqtype=seqtype, pat=pattern)
    else:
        data = oligo_count_dict(dir_count, k_len, rg, seqtype=seqtype, pat="_kmer.counts")
    lengths = get_len(dir_len)
    prob = probabilities(data, seq_len, k_len)
    chaos = chaos_game_representation(prob, k_len)
    pylab.title(f"Chaos game representation for {k_len}-mers")
    pylab.imshow(chaos, interpolation="nearest", cmap=cm.gray_r)
    pylab.savefig(os.path.splitext(file)[0]+"chaos3.png")
    pylab.show()
    end = time.time()
    # print some info
    print(colored(f"Total time for the script finishes: {round(end - start, 2)}.",
                  "green",
                  attrs=["bold"]))
    print(colored("Done!",
                  "red",
                  attrs=["bold"]))


if __name__ == "__main__":
    sys.exit(main())










    PATH = os.getcwd()
    filenames = glob.glob(f"{dir_in}/{sub_dir}/{ssub_dir}/*.{extension}")
    rg =  re.compile(r"(GC[AF]_\d+\.\d)")
    for file in filelist:
        f = open(file)
        s1 = f.read()
        data = "".join(s1.split("\n")[1:])
        f3 = oligo_count_dict(path_root, k, rg, seqtype="CHR", pat="_kmer.counts")
        f4 = oligo_count_dict(path_root, k, rg, seqtype="CHR", pat="_kmer.counts")

        f3_prob = probabilities(f3, 3)
        f4_prob = probabilities(f4, 4)

        chaos_k3 = chaos_game_representation(f3_prob, 3)
        pylab.title("Chaos game representation for 3-mers")
        pylab.imshow(chaos_k3, interpolation="nearest", cmap=cm.gray_r)
        pylab.savefig(os.path.splitext(file)[0]+"chaos3.png")
        pylab.show()

        chaos_k4 = chaos_game_representation(f4_prob, 4)
        pylab.title("Chaos game representation for 4-mers")
        pylab.imshow(chaos_k4, interpolation="nearest", cmap=cm.gray_r)
        pylab.savefig(os.path.splitext(file)[0]+"chaos4.png")
        pylab.show()
