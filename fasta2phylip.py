#! /usr/bin/env python
#
# v1
# v1.1 2022-10-21 python3 update

import sys
import argparse
from Bio import AlignIO

'''fasta2phylip.py example.aln

convert fasta alignments to relaxed phylip alignments
allows more than 10 character in names

or use -f for other formats, including:
'clustal' 'maf' 'nexus'
'''

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', help="fasta format file")
	parser.add_argument('-f', '--out-format', default="phylip-relaxed", help="output format, default: phylip-relaxed")
	parser.add_argument('-F', '--in-format', default="fasta", help="input format, default: fasta")
	parser.add_argument('--force-format', action="store_true", help="force 10-character phylip and no spaces")
	parser.add_argument('-r', '--reverse', action="store_true", help="use defaults for phylip to fasta")
	args = parser.parse_args(argv)

	if args.reverse:
		args.in_format = "phylip-relaxed"
		args.out_format = "fasta"
		if args.input_file.rsplit(".",1)[-1]=="fasta" or args.input_file.rsplit(".",1)[-1]=="aln":
			sys.stderr.write("\nWARNING input file does not appear to be phylip\n\n")

	with open(args.input_file, 'r') as input_handle:
		alignments = AlignIO.parse(input_handle, args.in_format)
		if args.force_format:
			for alignment in alignments:
				print( " {} {}".format( len(alignment), alignment.get_alignment_length() ), file = wayout)
				for seqrec in alignment:
					outname = str(seqrec.id) + "         "
					print( "{} {}".format(outname[0:10], str(seqrec.seq) ), file = wayout)
		else:
			AlignIO.write(alignments, wayout, args.out_format)

if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
