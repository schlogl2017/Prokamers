#! /usr/bin/env python
#
# counttetranucleotides.py v1
# program to read in fasta file and output tetranucleotides frequencies

import sys
import itertools
import argparse
import time
import string
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

usage='''
counttetranucleotides.py -n trimmed_seqs.fasta
or
counttetranucleotides.py -n trimmed_seqs.fasta > codon_table.txt
'''

def make_codon_dict(codonlist):
	# set values of all of them to zero, in case some codon is never used
	codondict = {c:0 for c in codonlist}
	return codondict

def reverse_complement(sequence,revcomptable):
	rs = sequence[::-1].translate(revcomptable)
	return rs

def make_rev_comps(quadcounter,kmerlist):
	# for fast string translation
	revcomptable = string.maketrans("ACGT","TGCA")
	# dictionary to sum all kmers and their reverse complements, as well as palendromes
	revcompdict = {}
	usedkmers = []
	for k in kmerlist:
		rc = reverse_complement(k, revcomptable)
		# for palendromic kmers
		if rc == k:
			revcompdict[k] = quadcounter.get(k,0)
			usedkmers.append(k)
		# as long as both are not counted yet, then count them together
		else:
			if k not in usedkmers and rc not in usedkmers:
				combikey = "%s+%s" % (k,rc)
				revcompdict[combikey] = quadcounter.get(k,0) + quadcounter.get(rc,0)
				usedkmers.extend([k,rc])
	return revcompdict

def print_last_dict(revcompdict,tetracount, wayout):
	for x in sorted(revcompdict.keys()):
		print >> wayout, "%s %d %.4f" % (x, revcompdict[x], revcompdict[x]*100/float(tetracount) )

def main(argv, wayout):
	#print argv
	if not len(argv):
		argv.append('-h')

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=usage)
	parser.add_argument('-n','--nucleotides', help='nucleotide fasta file')
	parser.add_argument('-q','--quiet', action='store_true', help='reduce most stderr output')
	args = parser.parse_args(argv)
	
	tetracount = 0
	seqcount = 0
	miscounts = 0

	# generate all combinations of 4 letters from ACTG, the order is normally random so they must be sorted
	allquads = [''.join(x) for x in sorted(set(itertools.product("ACTG",repeat=4)))]
	quadcounter = make_codon_dict(allquads)

	print >> sys.stderr, "Counting tetranucleotides", time.asctime()
	for seqrec in SeqIO.parse(args.nucleotides,'fasta'):
		seqcount += 1
		# this will add up all of the tetranucleotides to the counter and either finish or hit an error at any N's
		for x in range(0,len(seqrec.seq)-3):
			quad = str(seqrec.seq[x:x+4])
			try:
				quadcounter[quad]+=1
				tetracount += 1
			# in case some sequence is not in the dictionary, such as with N or X
			except KeyError:
				if not args.quiet:
					print >> sys.stderr, "Error of %s in %s" % (quad, seqrec.id)
				miscounts += 1
		# do not count those which make an error

	print >> sys.stderr, "%d sequences with %d 4-mers" % (seqcount,tetracount), time.asctime()
	#print >> sys.stderr, "Counted %d sequences with Ns" % (nbases)
	print >> sys.stderr, "Counted %d incomplete 4-mers" % (miscounts)

	print >> sys.stderr, "Combining reverse complements", time.asctime()
	revcompdict = make_rev_comps(quadcounter, allquads)

	print_last_dict(revcompdict,tetracount, wayout)


if __name__ == "__main__":
	main(sys.argv[1:], sys.stdout)
