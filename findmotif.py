#! /usr/bin/env python3
#
# v13.2 add gzip 2023-10-07
# v13 move to python3 2022-12-01
# v10 added comma to separate multiple motifs in same protein 2016-02-23
# v9 added primer mode and split motif search 2015-03-11
# v8 added trim and split options, allows for ambiguous DNA in fasta file 2015-03-05
# v7 added count of motifs, to track of multiple hits in one scaffold 2014-04-29
# v6 fixed bug for exclusion of reversed complements 2014-02-28
# v5 added ability to search reverse complement motifs 2013-10-28
# v4 changed nucl search to re, to conform to protein searching 2013/08/29
# v3 changed protein re search to allow coloration of matches 2013/08/28
# v2 redone to use standard python libraries and integrate in pipelines 20-mar-2013
#
# exclude or keep sequences that contain a motif

"""
FINDMOTIF.PY v13.2 last modified 2023-10-07

  print out fasta seqs if they contain (or not) a particular motif

  ### FOR NUCLEOTIDE MOTIFS ###
    use ambiguous DNA sequences for nucleotide searches:
    ATCGU - S[GC]W[AT]K[TG]M[AC]R[AG]Y[CT] - BDHV - N

    ambiguous nucleotides can be in motif or input file
    # for example:

findmotif.py -m GCSTACSYSNN -t n est_seqs.fasta

    to search the reverse complement as well for nucleotides, add -r
    -t n is not needed for unambiguous DNA with no reverse complement
    however for most nucleotide searches it is advisable to use -r and -t n

    use ',' to indicate a gap of unknown size, such as between two motifs:
-m TTGTCTTCAANGA,TAAAAGTCGTAACAA

    for actual primers, where one is a reverse complement (the reverse primer)
    use -p to search for the amplified region, for example:
    # V1-V2 bacterial 16S primers, 27F / 338R
findmotif.py -m AGAGTTTGATCCTGGCTCAG,GCTGCCTCCCGTAGGAGT -p -t n -r seqs.fa

  # V3-V4 prok 16S primers, Pro341F / Pro805R , from Takahashi 2014
    CCTACGGGNBGCASCAG,GACTACNVGGGTATCTAATCC
  # V3-V5 bacterial 16S primers , 338F / 907rBAC , from Zehr 2003
    ACWCCTACGGGWGGCAGCA,CCCGTCAATTCCTTTGAGTTT
  # V4 prok 16S primers, 515F-­Y / 806RB , from Pichler 2018
    GTGYCAGCMGCCGCGGTAA,GGACTACNVGGGTWTCTAAT
  # V4–V5 bacterial 16S, 515F-Y / 926R
    GTGYCAGCMGCCGCGGTAA,CCGYCAATTYMTTTRAGTTT
  # cpn60 chaperonin primers, H729 / H730
    GAIIIIGCIGGIGAYGGIACIACIAC,YKIYKITCICCRAAICCIGGIGCYTT

  # fungal 18S V1-V3 nu-SSU-0062-5 / nu-SSU-0531-3
    CCATGCATGTCTAAGTWTAA,CAATTGTTCCTCGTTAAG
  # fungal 18S V7-V9 nu-SSU-1333-5 / nu-SSU-1647-3
    CGATAACGAACGAGACCT,ANCCATTCAATCGGTANT
  # fungal 18S ITS3-ITS4 , from White 1990 
    GCATCGATGAAGAACGCAGC,TCCTCCGCTTATTGATAGC
  # fungal ITS (ITS5, ITS4)
    TCCTCCGCTTATTGATATGC,GGAAGTAAAAGTCGTAACAAGG

  ### FOR PROTEIN MOTIFS ###
    use normal protein code, and X for any amino acid:
    ACDEFGHIKLMNPQRSTVWY - X
    # for example, to find this conserved NRPS motif:

findmotif.py -m KIRGXRIEL prot_seqs.fasta

  ### GENERAL USAGE ###
    to trim around a motif, or to take pieces at each motif
    use -d and -s, noting that -s must be used with -d
    # for example:

findmotif.py -m CCTCATT -r -t n est_seqs.fasta -d 30 -s

    stdin can be used with - or can be gzipped
    results go to stdout
"""

import sys
import os.path
import re
import argparse
import time
import gzip
from Bio import SeqIO

def wrap_sequence(sequence, wrap):
	# space out the string with line breaks every n characters, each '\n' counts as one character
	spacedseq = "".join(sequence[i:i+wrap] + "\n" for i in range(0,len(sequence), wrap))
	return spacedseq

def make_motif_regex(motif, transtable):
	nucmotif = ''
	# cycle through dictionary and replace all ambiguous bases with regular expression characters
	for i in motif:
		try:
			nucmotif += transtable[i]
		# for cases with unusual characters, such as '+' or '*'
		except KeyError:
			nucmotif += i
	return nucmotif

def colorize_sequence(spacedseq, regexhits, wrap, offset=0, coloroff='\x1b[0m', redcolor='\x1b[31m'):
	seqtolist = list(spacedseq)
	for reseq in sorted(regexhits, key=lambda x: x.end(), reverse=True):
	# for each hit, insert the red and off escape per the positions from the regs
		offposition = reseq.regs[0][1]
		onposition = reseq.regs[0][0]
		# correct for trimming by offsetting by the position of the trimstart
		if offset:
			offposition -= offset
			onposition -= offset
		# insert at position defined from regex plus the wrap buffer
		# wrap buffer is effectively the number of lines to that position
		seqtolist.insert(offposition + offposition//wrap, coloroff)
		seqtolist.insert(onposition + onposition//wrap, redcolor)
	# rejoin the list as a string, and write
	joinedseq = ''.join(seqtolist)
	return joinedseq

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	# this handles the system arguments
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default = '-', help="fasta format file, can be .gz, can be stdin as -")
	parser.add_argument('-m', '--motif', help="motif, use X for any residue, use ',' to separate motifs in primer mode")
	parser.add_argument('-r', '--reversed', action='store_true', help="search the reverse complement as well for nucleotides")
	parser.add_argument('-t', '--type', choices="np", default="p", help="type - n (nucleotide) or p (protein - default)")
	parser.add_argument('-a', '--above', type=int, metavar='N', default=1, help="only keep sequences with more than N")
	parser.add_argument('-b', '--below', type=int, metavar='N', default=1000000, help="only keep sequences with fewer than N")
	parser.add_argument('-c', '--colorize', action='store_true', help="output the matches with color")
	parser.add_argument('-d', '--trim', type=int, help="trim matches at plus/minus N characters")
	parser.add_argument('-p', '--primer', action='store_true', help="matches to sequence between two motifs when the second one is reverse complemented, both primer seqs should be 5' to 3'")
	parser.add_argument('-f', '--fastq', nargs='?', const='fastq', default='fasta', help="flag to use fastq format sequences")
	parser.add_argument('-s', '--split', action='store_true', help="when using trim, split each found motif and include index as description")
	parser.add_argument('-v', '--verbose', action='store_true', help="verbose output")
	parser.add_argument('-w', '--wrap', default=60, type=int, help="characters per line, default 60")
	parser.add_argument('-x', '--exclude', action='store_true', help="instead exclude those with the motif")
	args = parser.parse_args(argv)

	seqcount = 0
	writecount = 0
	motifcount = 0
	formattype = args.fastq

	# strings for color escape signals
	redcolor = '\x1b[31m'
	greencolor = '\x1b[32m'
	boldon = '\x1b[1m'
	boldoff = '\x1b[22m'
	coloroff = '\x1b[0m'
	# because the stdin option allows multiple motifs via pipes,
	# if a motif is within 7-8 characters from the end of the line
	# the entire line may become colored due to disrupting the escape signals

	if args.split:
		if args.trim is None: # force default of 100bp if unspecified, but only for split mode
			args.trim = 100
			sys.stderr.write( "# split sequences at each motif, buffer -d not specified, using 100bp\n" )
		else:
			sys.stderr.write( "# split sequences at each motif, using {}bp\n".format(args.trim) )

	# generate the regular expression for the search
	if args.type == 'p':
		sys.stderr.write( "# Starting protein search for motif {}  {}\n".format( args.motif, time.asctime() ) )
		# set up regular expression for the protein motif
		protmotif = args.motif.replace("X","\w").replace("x","\w").replace(",","\w+")
		searchre = re.compile( protmotif )
	#only p and n are allowed, therefore type is otherwise 'n'
	else: # args.type == 'n'
		sys.stderr.write( "# Starting nucleotide search for motif {}  {}\n".format( args.motif, time.asctime() ) )
		# dictionary of ambiguous dna terms
		ambigdna = {'A':'A','T':'T','C':'C','G':'G','U':'U', 'S':'[GCS]', 'W':'[ATW]', 'K':'[TGK]', 'M':'[ACM]', 'R':'[AGR]', 'Y':'[CTY]', 'B':'[TCGB]', 'D':'[ATGD]', 'H':'[ATCH]', 'V':'[ACGV]', 'N':'[ATGCNSWKMRYBDHV]', ",":"[ATGCNSWKMRYBDHV]+"}
		# dictionary of dna terms for reverse complemented sequence,
		# hence should look opposite from above
		# A-TU C-G M-K Y-R V-B H-D are all switched, S W N are the same
		revambigdna = {'A':'T','T':'A','C':'G','G':'C','U':'A', 'W':'[GCS]', 'S':'[ATW]', 'M':'[TGK]', 'K':'[ACM]', 'Y':'[AGR]', 'R':'[CTY]', 'V':'[TCGB]', 'H':'[ATGD]', 'D':'[ATCH]', 'B':'[ACGV]', 'N':'[ATGCNSWKMRYBDHV]', ",":"[ATGCNSWKMRYBDHV]+"}

		# if using primer mode, reverse the second portion of the motif
		if args.primer:
			primersplit = args.motif.split(",")
			# primer mode requires only two motifs, any more or less are flagged as errors
			if len(primersplit) != 2:
				raise IndexError("Primer mode -p must have exactly two motifs separated by comma")
			# generate the regex motif as before, combining forward and reverse with any number of bases
			# the second primer is assumed to be reverse, so must be reverse complemented
			# a comma is added to the first string to allow for any number of bases in between the primers
			forwardregex = make_motif_regex(primersplit[0]+",", ambigdna) + make_motif_regex(primersplit[1][::-1], revambigdna)
			# same for the reverse search, where the forward primer must be reverse complemented
			# and the order must be switched
			reverseregex = make_motif_regex(primersplit[1]+",", ambigdna) + make_motif_regex(primersplit[0][::-1], revambigdna)
		# for normal mode
		else:
			forwardregex = make_motif_regex(args.motif, ambigdna)
			# read letters from motif backwards
			reverseregex = make_motif_regex(args.motif[::-1], revambigdna)
		# set up regular expression for the nucleotide motif
		searchre = re.compile(forwardregex)
		if args.reversed:
			# only generate reversed complement search if using reversed mode
			searchrc = re.compile(reverseregex)

	sys.stderr.write( "# Reading sequences from {}  {}\n".format( args.input_file.name, time.asctime() ) )
	if args.input_file.name.rsplit('.',1)[-1]=="gz":
		input_handler = gzip.open(args.input_file.name,'rt')
	else:
		input_handler = args.input_file
	for seq_record in SeqIO.parse(input_handler, formattype):
		seqcount += 1
		writeout = False
		if args.exclude:
			# if checking reversed complement, neither must contain the search to write the seq
			if args.reversed:
				if not searchre.search(str(seq_record.seq)) and searchrc.search(str(seq_record.seq)):
					writeout = True
			# otherwise, if only the motif is not found, then write the seq
			else:
				if not searchre.search(str(seq_record.seq)):
					writeout = True
		else:
			# normally returns iterator, so must make as list
			reresults = list(searchre.finditer(str(seq_record.seq)))
			# if searching the reverse complement, extend that list onto the first list
			if args.reversed and args.type=='n':
				reresults.extend(list(searchrc.finditer(str(seq_record.seq))))
			# count the total motif hits in this sequence
			motifcount += len(reresults)
			# if no results, ignore this sequence
			if len(reresults) >= args.above and len(reresults) < args.below:
				seqstring = str(seq_record.seq)
				# if splitting, behavior is different
				if args.split:
					# generate a sequence for each motif, rather than whole sequences
					for n, reseq in enumerate(reresults):
						# trim points are by each match, rather than global
						trimstart = reseq.start() - args.trim
						trimend   = reseq.end() + args.trim
						# if the trim position extends past the negative boundary, set start to 0
						if trimstart < 0:
							trimstart = 0
						splitstring = seqstring[trimstart:trimend]
						spacedseq = wrap_sequence(splitstring, args.wrap)
						# number each motif
						wayout.write(">{}_{}  {}-{}\n".format(seq_record.id, n, reseq.start(), reseq.end() ) )
						# detects if connected to a terminal or piped, or can set to color for piping
						if sys.stdout.isatty() or args.colorize:
							# colorize_sequence expects a list, so must make a list out of reseq
							# even though it is always only one item
							outseq = colorize_sequence(spacedseq, [reseq], args.wrap, trimstart)
							wayout.write(outseq)
						else:
							wayout.write(spacedseq)
						writecount += 1
				else:
					# find the first start position and the last end position as the boundaries for trimming
					# this interferes with the coloration because the index becomes wrong by trimming
					if args.trim:
						# start position is the minimum of all start positions across reresults
						# likewise end position is the max of all end positions across reresults
						trimstart = min([x.start() for x in reresults]) - args.trim
						trimend   = max([x.end() for x in reresults]) + args.trim
						# if the trim position extends past the negative boundary, set start to 0
						if trimstart < 0:
							trimstart = 0
						# then trim the sequence
						seqstring = seqstring[trimstart:trimend]
					else:
						trimstart = 0
					spacedseq = wrap_sequence(seqstring, args.wrap)
					# addition of line breaks must be done prior to adding the escape characters for color,
					# which otherwise would change where the breaks go and disrupt the output
					wayout.write( ">{}\n".format(seq_record.id) )
					# detects if connected to a terminal or piped, or can set to color for piping
					if sys.stdout.isatty() or args.colorize:
						outseq = colorize_sequence(spacedseq, reresults, args.wrap, trimstart)
						wayout.write(outseq)
					else:
						wayout.write(spacedseq)
					writecount += 1
		if writeout:
			writecount += 1
			wayout.write(seq_record.format(formattype))

	# if piping, then do not write out stderr lines
	if sys.stdout.isatty() or args.verbose:
		# in case motif is at the end, this will turn off coloring for all output
		sys.stderr.write(coloroff)
		sys.stderr.write( "# Counted {} sequences:  {}\n".format(seqcount, time.asctime() ) )
		sys.stderr.write( "# Wrote {} sequences\n".format(writecount) )
		sys.stderr.write( "# Found {} motifs\n".format(motifcount) )

			
if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
