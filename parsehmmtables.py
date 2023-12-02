#!/usr/bin/env python
#
# v1.0 parse series of hmmsearch tabular files

'''
parsehmmtables.py v1.0 2016-01-13

    usage:
parsehmmtables.py -i *.meNOG.tab -f sequences.fasta > meNOG.sorted_hits.tab

    -f sequences.fasta can be either the transcriptome or translated proteins

    set evalue with -e, 1e-15 gives about 1 hit per sequence on average
    in actual trees, close orthologs range from 1e-150 to 1e-300

    generate hmmsearch tabular output files for all NOGs:
for FILE in ~/db/meNOG_hmm/meNOG.*hmm; do hmmsearch --cpu 4 -E 1e-5 --domtblout hmmsearch/`basename $FILE`.meNOG.tab --noali $FILE transcript_peptides.fasta >> meNOG.log; done
'''

# hmmsearch tblout output appears like
#                                                                        --- full sequence ---- --- best 1 domain ---- --- domain number estimation ----
# target name           accession  query name                 accession    E-value  score  bias   E-value  score  bias   exp reg clu  ov env dom rep inc description of target
#   ------------------- ----------       -------------------- ---------- --------- ------ ----- --------- ------ -----   --- --- --- --- --- --- --- --- ---------------------
# m.13957  -          meNOG.ENOG410V42G.meta_raw -            1.5e-95  322.1   0.0   1.7e-95  321.8   0.0   1.0   1   0   0   1   1   1   1 m.13957
#
# NOG IDs can also have meNOG.ENOG410V44Q.clustalo_raw format
#
# note that tabular output for hmmsearch is not tab-delimited
# and must be parsed using regular expressions
#
# max evalue for hmmsearch is around 300, maybe 350

import sys
import os
import argparse
import time
import re
import glob
from collections import defaultdict

def fasta_to_names(fastafile):
	'''make dict where each sequence ID in a fasta file is a key and an empty list is the value'''
	namedict = {}
	for line in open(fastafile,'r'):
		if line[0]==">":
			seqid = line[1:].split(" ")[0]
			namedict[seqid] = []
	return namedict

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', nargs="*", help="hmmsearch tabular output files or directory")
	parser.add_argument('-a','--annotations', help="annotations tsv file")
	parser.add_argument('-t','--temp', help="directory for temporary files ['temp/']", default="temp/")
	parser.add_argument('-d','--domains', type=int, help="max allowed domains [20]", default=20)
	parser.add_argument('-e','--evalue', type=float, help="evalue cutoff [1e-100]", default=1e-100)
	parser.add_argument('-f','--fasta', help="fasta reference file, whatever contains the base names")
	parser.add_argument('-F','--full-sequence', action="store_true", help="use full sequence evalue instead of best domain")
	parser.add_argument('-H','--histogram', action="store_true", help="make histogram of evalues across all files")
	parser.add_argument('-s','--strict', action="store_true", help="exclude even single domain hits that are below evalue, otherwise keep weak hits if it is the only hit")
	parser.add_argument('-v','--verbose', action="store_true", help="verbose output")
	args = parser.parse_args(argv)

	EVALUERE = "[\d\.]+e-\d+" # this will not catch decimal values like 0.0025 or 0
	NOGRE = "meNOG\.(.+)\.\w+_raw"

	histodict = defaultdict(int)
	nogdict = defaultdict(int)

	if args.fasta: # make the dictionary with an empty list for each possible transcript
		print >> sys.stderr, "Reading reference fasta", time.asctime()
		namedict = fasta_to_names(args.fasta)
		print >> sys.stderr, "Reference fasta contained {} sequences".format(len(namedict) ), time.asctime()
	else: # if no fasta file is present, only make new lists each time a new transcript is found
		namedict = defaultdict(list)

	if os.path.isdir(args.input[0]):
		globstring = "{}*".format(args.input[0])
		inputfiles = glob.iglob(globstring)
	elif os.path.isfile(args.input[0]):
		inputfiles = args.input
	else:
		print >> sys.stderr, "Unknown input files, exiting", time.asctime()
		sys.exit()

	filecounter = 0
	searchcount = 0
	print >> sys.stderr, "Reading input files", time.asctime()
	for htf in inputfiles:
		filecounter += 1
		if args.verbose and filecounter%1000==0:
			print >> sys.stderr, filecounter, time.asctime()
		for line in open(htf,'r'):
			line = line.rstrip()
			if line and not line[0]=="#": # remove empty and comment lines
				searchcount += 1
				targetname = line.split(" ",1)[0] # always first term
				splitname = targetname.split("|")[0] # this sometimes causes problems related to TransDecoder naming
				nogquery = re.search(NOGRE, line).group(1)
				nogdict[nogquery] += 1
				evalues = re.findall(EVALUERE, line)
				if not evalues:
					bestevfloat = 0.0
				else:
					if args.full_sequence:
						bestevfloat = float(evalues[0])
					else:
						bestevfloat = float(evalues[-1])
				if args.strict and bestevfloat > args.evalue:
					continue
				try:
					namedict[splitname].append( (nogquery, bestevfloat) )
				except KeyError: # for cases where the | split removes critical key finding information
					namedict[splitname].append( (nogquery, bestevfloat) )
				# add all values to histogram regardless of evalue cutoff
				if args.histogram:
					try:
						bestevalue = int(evalues[-1].split("e-")[1])
					except IndexError: # empty list should only occur for evalue 0
						bestevalue = 350
					histodict[bestevalue] += 1
	print >> sys.stderr, "Parsed {} input files".format(filecounter), time.asctime()
	print >> sys.stderr, "Parsed {} hmm hits".format(searchcount), time.asctime()
	print >> sys.stderr, "From {} NOGs".format(len(nogdict) ), time.asctime()

	print >> sys.stderr, "Filtering domains", time.asctime()
	proteincount = 0
	nomatchcount = 0
	for sn in sorted(namedict.keys(), key=lambda x: int(x.split(".")[1]) ): ### TODO make sortable if ValueError if prots have non systematic names
		proteincount += 1
		hitlist = []
		if namedict[sn]: # if list is not empty
			for nepair in sorted(namedict[sn], key=lambda x: x[1]):
				if hitlist:
					if nepair[1] <= args.evalue and len(hitlist) < args.domains:
						hitlist.append("{}_{}".format(*nepair) )
				else:
					hitlist.append("{}_{}".format(*nepair) )
			print >> sys.stdout, "{}\t{}".format(sn, "|".join(hitlist) )
		else: # if no args.fasta, then this will never print, thus ignoring no matches
			nomatchcount += 1
			print >> sys.stdout, "{}\tNo_matches".format(sn)
	print >> sys.stderr, "Counted {} peptides".format(proteincount), time.asctime()
	print >> sys.stderr, "Counted {} sequences with no matches".format(nomatchcount), time.asctime()

	if args.histogram:
		hitcount = sum(histodict.values())
		print >> sys.stderr, "Generating histogram of evalues", time.asctime()
		for k in sorted(histodict.keys()):
			print >> sys.stdout, k, histodict[k], hitcount
			hitcount -= histodict[k]

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
