#!/usr/bin/env python


import sys
import glob
from fasta_parser import parse_fasta, get_fasta_headers


data = sys.argv[1]
ext = sys.argv[2]

if len(sys.argv) != 3:
    print('USAGE: python get_header.py < path to fasta files > < file extention >')
    sys.exit()
    
print(f'{data}/*.{ext}')
filenames = glob.glob(f'{data}/*.{ext}')
print(filenames)

for filename in filenames:
    print(get_fasta_headers(filename))


