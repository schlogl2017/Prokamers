#!/usr/bin/env python
#
# last updated to work with sizecutter v3 2014-03-13
# updated to skip errors due to empty files 2015-03-03
# python3 update 2023-04-10


'''get_multifile_size_info.py  last modified 2023-04-10

    collect size summary info from multiple fasta files

get_multifile_size_info.py file1.fasta file2.fasta file3.fasta > summary.tab
'''

import sys
import time
# requires additional script sizecutter.py
import sizecutter


def main(argv, wayout):
	input_list = argv[1:]
	
	print( "{}\t{}\t{}\t{}\t{}\t{}\t{}".format("Filename", "Sequences", "Total", "Mean", "Median", "Longest", "n50") , file=wayout )

	for input_file in input_list:
		sizecutter_output = sizecutter.main(["-n","-q",input_file], sys.stderr)
		if sizecutter_output:
			seqcount, seqmass, seqmean, seqmedian, seqn50, longest_seq = sizecutter_output
			print( "{}\t{}\t{}\t{:.1f}\t{}\t{}\t{}".format(input_file, seqcount, seqmass, seqmean, seqmedian, longest_seq, seqn50 ) , file=wayout )
			print( "Size information from {} written to {}:  {}".format(input_file, wayout, time.asctime()) , file=sys.stderr )
		else:
			print( "No sequences returned from {}:  {}".format(input_file, time.asctime()) , file=sys.stderr )

if __name__ == "__main__":
	main(sys.argv, sys.stdout)
