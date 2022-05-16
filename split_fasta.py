#!usr/bin/env python
#usage: split_fasta.py [-h] -p path [-d DIR_OUT] [-o OUTFILE]
# python splitting_fasta_files.py -p Data/bacteria/Acidisarcina -d Results/test/ -o plasmids/chromosome
# -*- coding: utf-8 -*-
import time
import os
import sys
import glob
import argparse
from termcolor import colored
from fasta_parser import split_sequences_from_fasta_file, get_name
from fasta_parser import write_fasta_file, punctuation_strip
from system_utils import make_me_a_folder, get_fasta_files


def get_data_to_fasta(dic, full_path, extension):
    name = ''.join([punctuation_strip(k) for k in dic])
    if os.path.exists(full_path):
        pass
    else:
        make_me_a_folder(full_path)
    write_fasta_file(dic, f'{full_path}/{name}.{extension}')


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
    # glob.glob('Bacteria/*.fna.gz')
    extension = opt.extension
    filenames = glob.glob(f'{dir_name}/*.{extension}')
    outfile = opt.outfile.split('/')
    dir_out = opt.dir_out
    if os.path.exists(dir_out):
        pass
    else:
        make_me_a_folder(dir_out)

    cnt_files = 0
    for filename in filenames:
        # Bacteria/GCA_000210675.1_ASM21067v1_genomic.fna.gz
        # ex GCA_000210675.1
        genome_id = filename.split('/')[-1][:15]
        print(colored(f"Splitting file: {genome_id}", attrs=['bold']))
        plasmids, chromosome = split_sequences_from_fasta_file(filename)
        full_path_plasmids = os.path.join(dir_out, genome_id, outfile[0])
        # saving the data
        if not plasmids:
            pass
        else: 
            get_data_to_fasta(plasmids, full_path_plasmids, 'fna')

        full_path_chromosome = os.path.join(dir_out, genome_id, outfile[1])
        # saving the data
        if not chromosome:
            pass
        else:
            get_data_to_fasta(chromosome, full_path_chromosome, 'fna')

        cnt_files += 1

    end = time.time()

    print(colored(f"Total number of splited files is {cnt_files}", attrs=['bold']))
    print(colored(f'The list of filenames is {len(filenames)}', attrs=['bold']))
    print(colored(f'Total time for the script: {end - start}', attrs=['bold']))
    print(colored('Done', attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
