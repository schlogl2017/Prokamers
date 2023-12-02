#! usr/bin/env python
# Script to read a cds file and change the name of the file and
# the header
import os
import sys
import gzip
import glob
import time
from termcolor import colored
import argparse
from Bio import SeqIO
from fasta_parser import parse_fasta


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
                        help="Directory to the input files.")

    parser.add_argument("-sd",
                        "--sub_dir",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="sub_dir",
                        help="Subdirectory of the input files.")

    parser.add_argument("-do",
                        "--dir_out",
                        metavar="path",
                        type=str,
                        required=False,
                        dest="dir_out",
                        help="Directory to the ouput files.")

    parser.add_argument("-e",
                        "--extension",
                        type=str,
                        required=False,
                        dest="extension",
                        help="File type.(fna.gz/fna)")
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
    # Genomes/Archaea
    dir_in = args.dir_in
    # Asgard/Cand_Heimdallarchaeota
    sub_dir = args.sub_dir
    # cds_from_genomic.fna.gz/CDS.fasta(Mito/Plastids)
    extension = args.extension
    # Genomes
    dir_out = args.dir_out
    # Genomes/ ""Archaea""
    supergroup = dir_in.split("/")[1]
    # "Asgard" /Cand_Heimdallarchaeota
    group = sub_dir.split("/")[0]
    # Asgard/ "Cand_Heimdallarchaeota"
    subgroup = sub_dir.split("/")[1]
    #Genomes/Archaea/Asgard/Cand_Heimdallarchaeota/CDS/GCA_003144275.1_ASM314427v1_cds_from_genomic.fna.gz
    filenames = glob.glob(f"{dir_in}/{sub_dir}/CDS/*_{extension}")
    for filename in filenames:
        print(filename)
        name = "_".join(filename.split("/")[-1].split("_")[:2])
        # Genomes/Archaea/Asgard/Cand_Heimdallarchaeota
        full_path = os.path.join(dir_out, supergroup, group, subgroup, "CDS")
        file_name = f"{name}_cds.fna"
        if not os.path.exists(full_path):
            os.makedirs(full_path)
        with open(f"{full_path}/{file_name}", "w") as fout:
            for rec in SeqIO.parse(gzip.open(filename, "rt"), "fasta"):
            #for Id, sequence in parse_fasta(filename ):
                #print(sequence[:30])
                #Id = Id.split(" ")[0]
                print(colored(f"Working with genome {name}\n{filename}\n", "grey"))
                fout.write(f">{rec.id}\n")
                fout.write(f"{rec.seq}\n")
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
