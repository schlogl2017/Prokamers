#!usr/bin/env python
# python splitting_fasta_files.py -p Data/bacteria/Acidisarcina -d Results/test/ -o plasmids/chromosome
# -*- coding: utf-8 -*-
import time
from termcolor import colored
import os
import sys
import glob
import argparse
from fasta_parser import split_sequences_from_fasta_file
from fasta_parser import write_fasta_file, punctuation_strip
from system_utils import make_me_a_folder, get_fasta_files


def parse_arguments():
    """Parse the command line arguments to the the split fasta script.
    Sets up the argparse command-line parser and calls it. These options can be accessed
    using opt.option. The resulting results of split fasta script are two
    python dictionaries (plasmids, chromosome) mapping genome/plasmids id with all
    sequences.
    """
    parser = argparse.ArgumentParser(description='A script to split compressed bacterial fasta files.')
    parser.add_argument('-p',
                        '--path',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='path',
                        help='Path to the files')
    parser.add_argument('-d',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')
    parser.add_argument('-e',
                        '--extension',
                        type=str,
                        dest='extension',
                        help='Extension of the files.Ex, .fna.gz')
    parser.add_argument('-o',
                        '--outfile',
                        type=str,
                        dest='outfile',
                        help='Name for output files.')
    return parser.parse_args()


def main():
    start = time.time()
    cwd = os.getcwd()

    print(colored(f'\nThe working directory: {cwd}\n', attrs=['bold']))

    opt = parse_arguments()
    dir_name = opt.path
    outfile = opt.outfile.split('/')
    dir_out = opt.dir_out
    extension = opt.extension
    filenames = filenames = glob.glob(f'{dir_name}/*.{extension}')
    print(filenames)
    
    if os.path.exists(dir_out):
        pass
    else:
        make_me_a_folder(dir_out)

    cnt_files = 0
    for filename in filenames:
        # name of the taxon directory, ie. Acidisarcina
        genome_id = filename.split('/')[-1][:15]
        print(colored(f"Results for file: {genome_id}", attrs=['bold']))
        plasmids, chromosome = split_sequences_from_fasta_file(filename)
        # checking the data obtained
        le_pl, le_ch = len(list(plasmids.values())), len(list(chromosome.values()))
        print(colored(f"Number of plasmids: {le_pl}", attrs=['bold']))
        print(colored(f"Number of chromosomes: {le_ch}", attrs=['bold']))
        if le_pl != 0:
            plasm_names = plasmids.keys()
            full_path_plasmids = os.path.join(dir_out, outfile[0])
            # checking if there are a path to save the data
            if not os.path.exists(full_path_plasmids):
                os.makedirs(full_path_plasmids)
            # saving the data
            for name in plasm_names:
                write_fasta_file(plasmids, f'{full_path_plasmids}' + '/' + name + '.fna')
        else:
            pass
        if le_ch != 0:
            chromosome_names = chromosome.keys()
            full_path_chromosome = os.path.join(dir_out, outfile[1])
            # checking if there are a path to save the data
            if not os.path.exists(full_path_chromosome):
                os.makedirs(full_path_chromosome)
            # saving the data
            for name in chromosome_names:
                write_fasta_file(chromosome, f'{full_path_chromosome}' + '/' + name + '.fna')
            cnt_files += 1

    end = time.time()

    print(colored(f"Total number of files: {cnt_files}", attrs=['bold']))
    print(colored(f'Total time for the script: {end - start}', attrs=['bold']))
    print(colored('Done', attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
