#!/usr/bin/python

import urllib.request
import os
import sys
import time


def write_fasta_file(sequence, name, out_file, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    sequence_dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    out_file : str
        Path to write the sequences to.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    with open(out_file, 'w') as fout:
        fout.write(f'>{name}\n')
        for i in range(0, len(sequence), wrap):
            fout.write(f'{sequence[i:i + wrap]}\n')
                

if len(sys.argv) != 3:
    print("USAGE: fetch_genome.py <genome_id_list> <out_dir>")
    sys.exit(1)

url_template = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=%s&rettype=fasta&retmode=text"

if os.path.exists(sys.argv[2]):
    pass
else:
    os.mkdir(sys.argv[2])

c = 0

for Id in open(sys.argv[1]):
    Id = Id.strip()
    if Id == "":
        continue

    sys.stdout.write(f"\nFetching genome {Id}\n")
    sys.stdout.flush()
    gbk_out_file = os.path.join(sys.argv[2], Id + ".fa")
    if os.path.exists(gbk_out_file):
        print("Already fetched")
        pass
    else:
        with open(gbk_out_file, "w") as fh:
            genome = urllib.request.urlopen(url_template % Id).read().decode("utf-8") 
            fh.write(genome)
    c += 1
    print("Done\n")
    time.sleep(1.0)
print(f'The number of downloaded files was {c}')


