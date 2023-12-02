#!/usr/bin/env python
# usage: gc_content_gc_skew.py [-h] -di path [-do DIR_OUT] [-sd SUB_DIR]
#                              [-w WINDOW] [-s STEP] [-t SEQ_TYPE]
# gc_content_gc_skew.py: error: the following arguments are required: -di/--dir_in
import os
import sys
import glob
from time import time
import argparse
from termcolor import colored
import pandas as pd
import fasta_parser
from system_utils import make_me_a_folder
from gc_gc_skew import *


def parse_arguments():
    """Parse the command line arguments to the gc_content_gc_skew script.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    The resulting results are csv files from each genome contained in the genus directory
    with a list of all kmers generated from these genomes (chromosomes and plasmids).
    """
    parser = argparse.ArgumentParser(
        description="""A script to calculate the GC content (G+C) and the GC skew (G-C/G+C)
         from bacterial genomes/plasmids.""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory root. In my case the name is conjugated with a subdir')

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Name for a subdirectory, ex., Chromosomes.')

    parser.add_argument('-w',
                        '--window',
                        type=int,
                        dest='window',
                        help='Integer representing the size of the sub sequence')

    parser.add_argument('-s',
                        '--step',
                        type=int,
                        dest='step',
                        help='Integer representing the size of the overlap sequence')

    parser.add_argument('-t',
                        '--seq_type',
                        type=str,
                        dest='seq_type',
                        help='Type of sequence, ex, chromosome or plasmid')
    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time()
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
    # min and max kmer length
    window = args.window
    step = args.step
    # type of sequence (chr/plsm
    seq_type = args.seq_type
    # get the fasta files
    # glob.glob(f'Data/Genomes_splitted/*/Chromosomes/*_chr.fna.gz')
    filenames = glob.glob(f'{dir_in}/*/{sub_dir}/*_chr.fna.gz')
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)
    # get the files by names and create a dataframe

    for filename in filenames:
        name = filename.split('/')[2]
        for Id, sequence in fasta_parser.parse_fasta(filename):
            print(f'Calculating the GC content from sequence {Id}')
            # returns a float number
            gc = gc_content(sequence)
            # returns a dictionary
            gc_slide_win = get_gc_content_slide_window(sequence, window, step)
            df_gc_slide_window = pd.DataFrame(gc_slide_win.items(), columns=['window', 'gc(%)'])
            # returns a dictionary
            gc_diff = difference_gc(gc, gc_slide_win)
            df_gc_window = pd.DataFrame(gc_diff.items(), columns=['window', 'diff'])
            df_gc_slw_dif_slw = pd.merge(df_gc_slide_window, df_gc_window, how='left', on='window')
            # returns a float number
            gc_var = get_chromosomal_gc_variation(gc_diff)
            var_series = pd.Series(gc_var, index=[name])
            gc_series = pd.Series(gc, index=[name])
            df_gc_var_gc_total = pd.DataFrame([var_series, gc_series],
                                              columns=[name]).T.rename(columns={0: 'gc_var',
                                                                                1: 'gc_total'}).reset_index().rename(
                columns={'index': 'name'})
            print(f'Writing the csv file from genus: {name}')
            full_path = os.path.join(dir_out, 'gc_data', name)
            csv1 = f'{Id}_gc_slw_gc_diff{seq_type}.csv'
            csv2 = f'{Id}_gc_var_gc_{seq_type}.csv'
            if os.path.exists(full_path):
                pass
            else:
                os.makedirs(full_path)
            df_gc_slw_dif_slw.to_csv(f'{full_path}/{csv1}', index=False)
            df_gc_var_gc_total.to_csv(f'{full_path}/{csv2}', index=False)
    end = time()
    # print some info
    print(colored(f"Total number of genus/species analyzed: {len(filenames)}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
