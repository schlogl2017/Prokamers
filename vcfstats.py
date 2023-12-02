#!/usr/bin/env python
#
# vcfstats.py v1.0 2015-07-28
#
# reference for vcf from:
# http://www.1000genomes.org/wiki/analysis/variant%20call%20format/vcf-variant-call-format-version-41

'''
vcfstats.py v1.2 2022-01-02

vcfstats.py -i snps.vcf -g genomic_contigs.fasta -t transcripts.gtf -E exclude_list

    TO MAKE VCF FILE WITH FREEBAYES:
freebayes -f wga.fasta -u wga_bowtie2.sorted.bam | vcffilter -f "QUAL > 20" > wga_bowtie2.bam.sorted.vcf

    TO REQUIRE AT LEAST 10 OF 30 READS CONTAIN THE SNP:
freebayes -f SCIL_WGA_130802.fasta -i -X -u scil_wga_pe_370_600.sorted.bam | vcffilter -f "DP > 30 & AO > 10" > scil_wga_pe_370_600.dp30.ao10.vcf
    OR TO FILTER BY QUALITY SCORE (MAY PREDICT TOO MANY IN LOW COV REGIONS)
vcffilter -f "QUAL > 15" > scil_wga_pe_370_600.q15.vcf

    TO PREPARE SORTED BAM FILE AND INDEX:
bowtie2 -q -p 4 --no-discordant --no-sq --no-unal -x wga_build -1 R1.fastq.gz -2 R2.fastq.gz | samtools view -bT wga.fasta - | samtools sort - wga_bowtie2.sorted
samtools index wga_bowtie2.sorted.bam

    to add CDS evaluation, add -T transdecoder.gff3
    and write STDOUT to file with:
-T protein_models.gff3 > vcf.mutations.tab
    use -B to include BLOSUM matrix for mutation scoring
    FIND BLOSUM MATRIX AT:
http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt

    STDOUT in tabular format of 14 columns, with optional 15th column:
contig  transcript  position  refbase  altbase  strand  position-in-transcript
position-in-codon (0,1,2) ref-codon  alt-codon (rev-comp'd for - strand)
ref-AA  alt-AA  position-in-protein  Syn-or-Non
    optionally, BLOSUM score of the change
'''

#
import sys
import argparse
import time
import re
import gzip
from itertools import groupby,chain
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
#

BLOSUM62 = """#  Matrix made by matblas from blosum62.iij
#  * column uses minimum score
#  BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
#  Blocks Database = /data/blocks_5.0/blocks.dat
#  Cluster Percentage: >= 62
#  Entropy =   0.6979, Expected =  -0.5209
   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 
R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 
N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 
D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 
C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 
Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 
E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 
H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 
I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 
L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 
K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 
M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 
F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 
P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 
S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 
T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 
W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 
Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 
V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 
B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 
Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 
X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 
* -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 
"""

def ranges_to_sets(rangelist):
	'''convert list of tuples to set of base positions'''
	basepositionset = set([])
	for bounds in rangelist:
		basepositionset.update(range(bounds[0], bounds[1]+1 ) )
	return basepositionset

def combine_intervals(rangelist):
	'''convert list of tuples to non redundant intervals'''
	# sl = [(1,10), (1,6), (15,20), (1,10), (19,29), (30,35), (6,13), (40,48), (15,21), (42,50)]
	nrintervallist = []
	srtrangelist = sorted(rangelist) # sort list now, to only do this once
	interval = srtrangelist[0] # need to start with first time, which will be the same for the first bounds
	for bounds in srtrangelist:
		# since it is sorted bounds[0] should always be >= interval[0]
		if bounds[0] > interval[1]+1: # if the next interval starts past the end of the first + 1
			nrintervallist.append(interval) # add to the nr list, and continue with the next
			interval = bounds
		else:
			if bounds[1] > interval[1]: # bounds[1] <= interval[1] means do not extend
				interval = (interval[0], bounds[1]) # otherwise extend the interval
	else: # append last interval
		nrintervallist.append(interval)
	# should return [(3, 13), (15, 35), (40, 50)]
	return nrintervallist

def numbers_to_ranges(numberset):
	'''makes groups of value-positions, where values in a row end up as the same group
	n = [11,12,13,27,28,29]
	list(numbers_to_ranges(n))
	[(11, 13), (27, 29)]'''
	# enumerate() makes tuples of (0,11), (1,12), (2,13), (3,27), (4,28), (5,29)
	# list l is then 
	# [(0, 11), (1, 12), (2, 13)]
	# [(3, 27), (4, 28), (5, 29)]
	# so x[1]-x[0] gives 11,11,11,24,24,24
	# [(11, 13), (27, 29)]
	for i,l in groupby(enumerate(sorted(numberset)), lambda x: x[1]-x[0] ):
		l = list(l)
		print("{}".format(l) )
		# returns generator of first and last positions in the group
		yield l[0][1], l[-1][1]

def get_occupied_length(rangedict): # added in v1.3 for speed and memory improvements
	'''from a dictionary of list, where list items are tuples of intervals, return the sum of non redundant intervals'''
	# rangedict must be dictionary of lists, where list items are boundaries of features
	sys.stderr.write( "# Determining occupied bases  {}\n".format( time.asctime() ) )
	occlen = 0
	for rl in rangedict.values(): # rl in rangedict.values() should be as [(1,5), (25,29)]
		for nrintvl in combine_intervals(rl):
			occlen += nrintvl[1]+1-nrintvl[0] # calculate adjusted interval length
	return occlen

def sort_intron_exon(transcriptgtf, excludecontigs):
	'''read in GTF file and return dictionaries for introns and exons where keys are scaffolds and values are lists of intervals'''
	exonbyscaffold = defaultdict(list)
	genesbyscaffold = defaultdict(list)
	genesum, exonsum, intronsum = 0,0,0
	commentlines = 0
	sys.stderr.write( "# Reading exons from {}  {}\n".format(transcriptgtf, time.asctime() ) )
	for line in open(transcriptgtf,'r'):
		line = line.rstrip()
		if line: # ignore empty lines
			if line[0]=="#": # count comment lines, just in case
				commentlines += 1
			else:
				lsplits = line.split("\t")
				scaffold = lsplits[0]
				if excludecontigs and excludecontigs.get(scaffold, False):
					continue # skip anything that hits to excludable scaffolds
				feature = lsplits[2]
				boundaries = (int(lsplits[3]), int(lsplits[4]) ) # tuple of ints
				if feature=="transcript" or feature=="gene":
					genesbyscaffold[scaffold].append(boundaries)
				elif feature=="exon" or feature=="CDS": # allow CDS if exons are not present, such as for Augustus
					exonbyscaffold[scaffold].append(boundaries)
	if commentlines:
		sys.stderr.write( "# Counted {} comment lines  {}\n".format(commentlines, time.asctime() ) )
	exonindex = defaultdict(list)
	intronindex = defaultdict(list)
	sys.stderr.write( "# Determining exon and intron positions  {}\n".format( time.asctime() ) )
	for scaf, boundarylist in exonbyscaffold.items():
		baseset = ranges_to_sets(boundarylist) # get set of exon positions
		exonsum += len(baseset) # do counting now
		exonindex[scaf] = combine_intervals(boundarylist)
		geneset = ranges_to_sets(genesbyscaffold[scaf])
		genesum += len(geneset)
		intronset = geneset.difference(baseset) # intron is gene positions without exon positions
		intronsum += len(intronset)
		intronindex[scaf] = list( numbers_to_ranges(intronset) ) # merge introns into intervals
	sys.stderr.write( "# Counted bases for genes:{}, exons:{}, introns:{}  {}\n".format(genesum, exonsum, intronsum, time.asctime() ) )
	if intronsum != genesum - exonsum:
		sys.stderr.write( "# WARNING: intron sum {} does not equal {}\n".format(intronsum, genesum - exonsum) )
	return intronindex, exonindex, genesum, exonsum, intronsum

def find_in_interval(position, intervallist):
	for interval in sorted(intervallist):
		if interval[1] >= position >= interval[0]:
			return True
	else:
		return False

def get_exon_classes(exonclassfile, excludecontigs):
	'''read tabular exon info and return two dicts and two ints'''
	retintranges = defaultdict(list)
	skipexranges = defaultdict(list)
	retintsum = 0
	skipsum = 0
	excounter = 0
	sys.stderr.write( "# Reading special exon type positions from {}  {}\n".format( exonclassfile, time.asctime() ) )
	for line in open(exonclassfile,'r'):
		line = line.strip()
		if line and not line[0]=="#":
			excounter += 1
			# contig_715	gene.18005	5446	5507	IR
			lsplits = line.split('\t')
			scaffold = lsplits[0]
			exbounds = tuple( sorted( [ int(lsplits[2]), int(lsplits[3]) ] ) ) # order the two items
			featurelength = exbounds[1]+1-exbounds[0]
			exclass = lsplits[4]
			if exclass == "IR":
				retintranges[scaffold].append(exbounds)
				retintsum += featurelength
			else: # assume either skipped exon or cassette exon
				skipexranges[scaffold].append(exbounds)
				skipsum += featurelength
	sys.stderr.write( "# Counted {} special type exons  {}\n".format( excounter, time.asctime() ) )
	return retintranges, skipexranges, retintsum, skipsum

def count_contig_sum(genomiccontigs, exclusiondict=None, makedict=False):
	'''return total length of all contigs, except excludable ones'''
	sys.stderr.write( "# Parsing {}  {}\n".format(genomiccontigs, time.asctime() ) )
	contigsum = 0
	contigcount = 0
	if makedict:
		contigdict = {}
	for seqrec in SeqIO.parse(genomiccontigs,'fasta'):
		contigcount += 1
		if exclusiondict and exclusiondict.get(seqrec.id, False):
			continue
		contigsum += len(seqrec.seq)
		if makedict:
			contigdict[seqrec.id] = seqrec
	sys.stderr.write( "# Counted {} sequences  {}\n".format(contigcount, time.asctime() ) )
	sys.stderr.write( "# {} total bases\n".format(contigsum) )
	if makedict:
		return contigsum, contigdict
	return contigsum

def strand_from_exon(firstexon):
	'''get strand order from the first exon in a transcript'''
	if firstexon[0] < firstexon[1]: # then + strand
		return False
	else: # otherwise -, meaning reverse sort
		return True

def get_codon_stats(position, cdsbd, exonorder, reversestrand, scaffoldrec, refbase, altbase, blosummatrix, codonspan=2, debug=False ):
	'''takes transcript coding information and mutation information, then determines parameters for the SNP and returns a list (of strings) for output, and an integer of the codon position'''
	stringtopos = get_cds_str(position, cdsbd, exonorder, reversestrand, scaffoldrec) # string to snp, also offset by codonspan, which is 2
	posintrans = len(stringtopos) - codonspan # thus position in transcript is length of string minus codonspan
	snpposition = len(stringtopos)%3 # codon position
	codonoffset = -3 - snpposition
	aaposition = len(stringtopos)//3 # should be remainderless division by 3
	if debug:
		if reversestrand:
			print( ">{}_{}_rc{}_rc{}".format(position, len(stringtopos)%3, refbase, altbase) , file=sys.stdout )
		else:
			print( ">{}_{}_{}_{}".format(position, len(stringtopos)%3, refbase, altbase) , file=sys.stdout )
		print( "{}".format(stringtopos) , file=sys.stdout )
	refcodon = stringtopos[codonoffset:][:3]
	### TODO allow multiple snps in a list for double mutations, rather than one line for each
	altcodon = make_alt_codon(refcodon, altbase, snpposition, reversestrand)
	refaa = str(refcodon.translate() )
	altaa = str(altcodon.translate() )
	strandkey = "-" if reversestrand else "+"
	# output should be:
	# position on contig, base, base, strand, position in transcript, position in codon, codon, codon, AA, AA, position in protein
	outlist = [str(position), refbase, altbase, strandkey, str(posintrans), str(snpposition), str(refcodon), str(altcodon), refaa, altaa, str(aaposition)]
	outlist.append( "Syn" if refaa==altaa else "Non")
	if blosummatrix:
		outlist.append( str(blosummatrix[refaa][altaa] ) )
	return outlist, snpposition

def get_cds_str(position, matchedexon, exonorder, reversestrand, scaffoldrec, codonspan=2):
	'''given the position of the SNP, the coding exon that it matches
    the list of exons in order, the strand, and the sequence of the scaffold,
    generate the part of the transcript that is coding up to the SNP'''
	cdsstring = ""
	for exon in exonorder: ### TODO check if exons overlap, in case of two proteins called
		if exon==matchedexon: # indicates last exon up to the SNP
			if reversestrand:
				if position < 3: # for cases where index is 0,1 or 2, and position is negative, which returns an empty string
					cdsstring += scaffoldrec.seq[position-1:matchedexon[0]].reverse_complement() # omit codonspan
					cdsstring += "NN" # add the missing bases to allow for possible codon
				else:
					cdsstring += scaffoldrec.seq[position-1-codonspan:matchedexon[0]].reverse_complement()
			else:
				cdsstring += scaffoldrec.seq[matchedexon[0]-1:position+codonspan]
			return cdsstring
		else: # take all other exons as normal
			if reversestrand:
				cdsstring += scaffoldrec.seq[exon[1]-1:exon[0]].reverse_complement()
			else:
				cdsstring += scaffoldrec.seq[exon[0]-1:exon[1]]

def make_alt_codon(refcodon, altbase, snpposition, reversestrand):
	'''generate the alternate codon'''
	altstring = ""
	for i in range(3):
		if i==snpposition:
			if reversestrand: # must replace the reverse complement of the alt base
				altstring += str(Seq(altbase).reverse_complement())
			else:
				altstring += altbase
		else:
			altstring += str(refcodon)[i]
	return Seq(altstring)

def map_blosum(blosumfile):
	'''convert BLOSUM file to dictionary of dictionaries for each amino acid change'''
	blosummatrix = {}
	if os.path.isfile(blosumfile):
		sys.stderr.write( "# Reading BLOSUM matrix from {}  {}\n".format(blosumfile, time.asctime() ) )
		blosumiter = open(blosumfile,'r').readlines()
	else: # including if blosumfile is None
		sys.stderr.write( "# No file, using BLOSUM62 matrix  {}\n".format( time.asctime() ) )
		blosumiter = BLOSUM62.split("\n")
	for line in blosumiter:
		if line[0]=="#": # skip comment lines
			continue
		if line[0]==" ":
			aaorder = [line[i:i+3].strip() for i in xrange(1,len(line[1:-1]),3)]
		else:
			aa = line[0]
			mutvals = [int(line[i:i+3].strip()) for i in xrange(1,len(line[1:-1]),3)]
			blosummatrix[aa] = dict(zip(aaorder,mutvals))
	return blosummatrix

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', help="vcf file", required=True)
	parser.add_argument('-g','--genome', help="fasta format genomic contigs")
	parser.add_argument('-t','--transcripts', help="gtf/gff file of transcripts or gene models")
	parser.add_argument('-E','--exclude', help="file of list of contigs to exclude")
	parser.add_argument('-G','--ignore-gene', action="store_true", help="skip lines where type is gene")
	parser.add_argument('-H','--histogram', help="name for optional coverage histogram")
	parser.add_argument('-a','--AO', default=0, type=int, help="minimum AO value - occurrence of SNP")
	parser.add_argument('-d','--DP', default=0, type=int, help="minimum DP value - total reads mapped at that locus")
	parser.add_argument('-q','--qual', type=float, default=0.0, help="QUAL value must be greater than [0.0]")
	parser.add_argument('-r','--ref-count', default=1, type=int, help="minimum number of events to count alleles [1]")
	parser.add_argument('-T','--transdecoder', help="transdecoder genome gff file")
	parser.add_argument('-A','--augustus', help="AUGUSTUS gff file")
	parser.add_argument('-B','--blosum', help="BLOSUM matrix file, default is BLOSUM62")
	parser.add_argument('--exon-classes', help="optional tabular file of exon classes from splicevariantstats.py")
	args = parser.parse_args(argv)

	commentlines = 0

	if args.exclude:
		sys.stderr.write( "# Reading exclusion list {}  {}\n".format(args.exclude, time.asctime() ) )
		exclusiondict = {}
		for term in open(args.exclude,'r'):
			term = term.rstrip()
			if term[0] == ">":
				term = term[1:]
			exclusiondict[term] = True
		sys.stderr.write( "# Found {} contigs to exclude  {}\n".format( len(exclusiondict) , time.asctime() ) )
	else:
		exclusiondict = None

	sys.stderr.write( "# Parsing transcriptome {}  {}\n".format(args.transcripts, time.asctime() ) )
	intronindex, exonindex, genesum, exonsum, intronsum = sort_intron_exon(args.transcripts, exclusiondict)

	if args.exon_classes:
		irindex, skippedindex, irsum, sksum = get_exon_classes(args.exon_classes, exclusiondict)
	else:
		irindex = defaultdict(list) # empty dicts to search always returns False
		skippedindex = defaultdict(list)
		irsum, sksum = 0, 0

	if args.transdecoder or args.augustus:
		if args.transdecoder:
			proteingff = args.transdecoder
		else: # assume args.augustus, so only take augustus if transdecoder is not there
			proteingff = args.augustus
		sys.stderr.write( "# Parsing proteins from {}  {}\n".format(proteingff, time.asctime() ) )
		featurecounts = defaultdict(int)
		genebyRange = defaultdict(dict)
		CDSbyTranscript = defaultdict(list)
		for line in open(proteingff, 'r'):
			line = line.rstrip()
			if line: # ignore empty lines
				if line[0]=="#": # count comment lines, just in case
					commentlines += 1
				else:
					lsplits = line.split("\t")
					scaffold = lsplits[0]
					if args.exclude and exclusiondict.get(scaffold, False):
						continue # skip anything that hits to excludable scaffolds
					feature = lsplits[2]
					featurecounts[feature] += 1
					attributes = lsplits[8]
					if feature=="gene" and not args.augustus and not args.ignore_gene:
						transid = re.search('ID=([\w.-]+)\|?', attributes).group(1)
						genebounds = (int(lsplits[3]), int(lsplits[4]))
						genebyRange[scaffold][genebounds] = transid
					elif feature=="mRNA" and args.ignore_gene: # for TransDecoder
						transid = re.search('ID=([\w.-]+)\|?', attributes).group(1)
						genebounds = (int(lsplits[3]), int(lsplits[4]))
						genebyRange[scaffold][genebounds] = transid
					elif feature=="transcript" and args.augustus:
						transid= re.search('ID=([\w.-]+);', attributes).group(1)
						genebounds = (int(lsplits[3]), int(lsplits[4]))
						genebyRange[scaffold][genebounds] = transid
					elif feature=="CDS":
						if args.augustus:
							transid= re.search('ID=([\w.-]+).cds;', attributes).group(1)
						else:
							transid = re.search('ID=[cC][dD][sS].([\w.-]+)\|', attributes).group(1)
						if lsplits[6]=="+":
							cdsbounds = (int(lsplits[3]), int(lsplits[4]))
						else:
							cdsbounds = (int(lsplits[4]), int(lsplits[3]))
						CDSbyTranscript[transid].append(cdsbounds)
		if commentlines:
			sys.stderr.write( "# Counted {} comment lines  {}\n".format( commentlines, time.asctime() ) )
		for k in sorted(featurecounts.keys()):
			print( "{}\t{}".format(k, featurecounts[k]) , file=sys.stderr )
		codonpositioncounts = defaultdict(int)
		altbasecounts = defaultdict(int)

	# get overall genome stats
	if args.genome:
		if args.transdecoder or args.augustus:
			contigsum, genomedict = count_contig_sum(args.genome, exclusiondict, makedict=True)
			cdssetbyscaffold = defaultdict(set)
		else:
			contigsum = count_contig_sum(args.genome, exclusiondict)
			genomedict = None
	else:
		contigsum, genomedict = 0, None

	# get BLOSUM matrix from file
	blosummatrix = map_blosum(args.blosum) if args.blosum else None

	# DEBUG
	debug = False
	# read in vcf file and evaluate SNPs
	snpcounts = defaultdict(int)	
	typecounts = defaultdict(int)
	snptotal = 0
	snpsinexons = defaultdict(int)
	snpsinretint = defaultdict(int)
	snpsinskexs = defaultdict(int)
	snpsinintrons = defaultdict(int)
	snpsintergenic = defaultdict(int)
	excludesnps, lowqualsnps = 0,0
	exoncovdict, introncovdict, igenecovdict = defaultdict(int),defaultdict(int),defaultdict(int)

	if args.input.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write( "# Parsing {} as gzipped  {}\n".format(args.input, time.asctime() ) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write( "# Parsing {}  {}\n".format(args.input, time.asctime() ) )
	if args.AO or args.DP:
		sys.stderr.write( "# Finding SNPs counted at least {} times with total coverage at least {}\n".format(args.AO,args.DP) )
	for line in opentype(args.input, 'rt'):
		line = line.rstrip()
		if line: # ignore empty lines
			if line[0]=="#": # count comment lines, just in case
				commentlines += 1
			else:
				lsplits = line.split("\t")
				#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	unknown
				# contig1234	54321	.	A	G,T	150	.

				# INFO fields look like this:
				# AB=0;ABP=0;AC=0;AF=0;AN=2;AO=2;CIGAR=1X;DP=5;DPB=5;DPRA=0;
				# EPP=7.35324;EPPR=3.73412;GTI=0;LEN=1;MEANALT=1;MQM=8;MQMR=15.3333;NS=1;
				# NUMALT=1;ODDS=1.37774;PAIRED=0;PAIREDR=0;PAO=0;PQA=0;PQR=0;PRO=0;QA=13;QR=56;
				# RO=3;RPL=0;RPP=7.35324;RPPR=9.52472;RPR=2;RUN=1;SAF=2;SAP=7.35324;SAR=0;SRF=2;
				# SRP=3.73412;SRR=1;TYPE=snp
				try:
					info = dict([(field.split("=")) for field in lsplits[7].split(";")]) # make fast dict of INFO
				except ValueError: # dictionary update sequence element #0 has length 1; 2 is required
				# INFO for samtools mpileup
				# INDEL;IDV=8;IMF=0.571429;DP=14;VDB=0.000671399;SGB=-0.662043;MQ0F=0.0714286;
				# AF1=1;AC1=2;DP4=1,0,9,0;MQ=21;FQ=-54.5252;PV4=1,0.0214827,1,1
					info = dict([(field.split("=")) for field in lsplits[7].split(";")[1:]])
					typecounts[lsplits[7].split(";")[0]] += 1
				featuretype = info.get("TYPE",None)
				if featuretype:
					typecounts[info["TYPE"]] += 1
					if info["TYPE"]!="snp": # double SNPs TYPE=snp,snp are excluded here
						continue # currently ignore anything that is not a SNP
				### TODO potentially do filtering steps here instead of vcffilter
				scaffold = lsplits[0]
				snptotal += 1
				if args.exclude and exclusiondict.get(scaffold, False):
					excludesnps += 1
					continue
				position = int(lsplits[1])
				refbase = lsplits[3]
				### TODO potentially deal with double matches
				### TODO one of the two bases is the real one, and the refbase is wrong
				altbase = lsplits[4][0] # for the moment just take the first one
				qualscore = float(lsplits[5])
				if qualscore <= args.qual:
					lowqualsnps += 1
					continue
				if args.AO: # if filtering by AO or DP4
					try: # for freebayes format
						aoval = int(info["AO"])
					except KeyError: # meaning no "AO", possibly from samtools mpileup
						aoval = sum( [ int(i) for i in info["DP4"].split(",")[2:4] ] )
					if aoval < args.AO:
						lowqualsnps += 1
						continue
				if args.DP and int(info["DP"]) < args.DP:
					lowqualsnps += 1
					continue
				snpcounts[scaffold] += 1
				if find_in_interval(position, exonindex[scaffold]):
					snpsinexons[scaffold] += 1
					exoncovdict[int(info["DP"])] += 1
					if find_in_interval(position, irindex[scaffold]): # check subset of retained introns
						snpsinretint[scaffold] += 1
					elif find_in_interval(position, skippedindex[scaffold]): # check subset of skipped exons
						snpsinskexs[scaffold] += 1
				elif find_in_interval(position, intronindex[scaffold]):
					snpsinintrons[scaffold] += 1
					introncovdict[int(info["DP"])] += 1
				else:
					snpsintergenic[scaffold] += 1
					igenecovdict[int(info["DP"])] += 1
				# check if SNP is in CDS
				### TODO needs to deal with double mutations somehow
				### TODO possibly in relation to mapped reads rather than single SNPs
				if args.transdecoder or args.augustus:
					for exbd,tid in genebyRange[scaffold].items():
						if position >= exbd[0] and position <= exbd[1]:
							# means SNP is within gene
							#try:
							reversestrand = strand_from_exon(CDSbyTranscript[tid][0]) # get strand from first exon
							#except IndexError:
							#	sys.stderr.write( CDSbyTranscript[tid] )
							#	continue
							exonorder = sorted(CDSbyTranscript[tid], reverse=reversestrand)
							for cdsbd in CDSbyTranscript[tid]:
								# if reverse strand and within inverted transcript bounds
								# or not reverse strand (forward) and within transcript bounds
								if (reversestrand and position <= cdsbd[0] and position >= cdsbd[1]) or (not reversestrand and position >= cdsbd[0] and position <= cdsbd[1]):
									outlist = [scaffold, tid]
									statlist, snppos = get_codon_stats(position, cdsbd, exonorder, reversestrand, genomedict.get(scaffold,None), refbase, altbase, blosummatrix)
									outlist.extend(statlist)
									codonpositioncounts[snppos] += 1 # count codon positions of mutations
									if snppos==2: # for third codon position, count the mutations
										altbasecounts[(refbase,altbase)] += 1
									print( "\t".join(outlist ) , file=wayout )
	if commentlines:
		sys.stderr.write( "# Counted {} comment lines  {}\n".format(commentlines, time.asctime() ) )
	for k in sorted(typecounts.keys()):
		print( "{}\t{}".format(k, typecounts[k]) , file=sys.stderr)

	# PRINT FINAL OUTPUT
	sys.stderr.write( "Total SNPs: {}  {}\n".format(snptotal, time.asctime() ) )
	validsnps = sum(snpcounts.values())
	if lowqualsnps:
		print( "{} SNPs discarded for QUAL scores of {} or below".format(lowqualsnps, args.qual) , file=sys.stderr)
	print( "{} SNPs on kept contigs".format(validsnps) , file=sys.stderr)
	if contigsum:
		print( "{:.3f} SNPs per kb".format(validsnps * 1000.0/contigsum) , file=sys.stderr)
	exonsnptotal = sum(snpsinexons.values())
	print( "{} SNPs in exons".format(exonsnptotal) , file=sys.stderr)
	exonsum = exonsum - irsum - sksum # to prevent double counting
	print( "{:.3f} SNPs per kb of exon".format(exonsnptotal * 1000.0/exonsum) , file=sys.stderr)
	if irsum: # should be 0 if no exon classes are given
		retintsnptotal = sum(snpsinretint.values())
		print( "{} SNPs in retained introns".format(retintsnptotal) , file=sys.stderr)
		print( "{:.3f} SNPs per kb of retained introns".format(retintsnptotal * 1000.0/irsum) , file=sys.stderr)
	if sksum: # also should be 0 if no exon classes are given
		skexsnptotal = sum(snpsinskexs.values())
		print( "{} SNPs in skipped exons".format(skexsnptotal) , file=sys.stderr)
		print( "{:.3f} SNPs per kb of skipped exons".format(skexsnptotal * 1000.0/sksum) , file=sys.stderr)
	intronsnptotal = sum(snpsinintrons.values())
	if genomedict:
		cdssnptotal = sum(codonpositioncounts.values())
		print( "{} SNPs in coding sequences across all transcripts, including double counts".format(cdssnptotal) , file=sys.stderr)
		for k in sorted(codonpositioncounts.keys()):
			print( "base {}\t{}".format(k, codonpositioncounts[k]) , file=sys.stderr)
		for k in sorted(altbasecounts.keys()):
			if altbasecounts[k] > args.ref_count:
				print( "ref {1}\tto {2}\t{0}".format(altbasecounts[k], *k) , file=sys.stderr)
	print( "{} SNPs in introns".format(intronsnptotal) , file=sys.stderr)
	print( "{:.3f} SNPs per kb of intron".format(intronsnptotal * 1000.0/intronsum) , file=sys.stderr)
	intergenicsnptotal = sum(snpsintergenic.values())
	print( "{} SNPs from intergenic regions".format(intergenicsnptotal) , file=sys.stderr)
	if args.genome:
		intergenicsum = contigsum-genesum
		print( "{:.3f} SNPs per kb of intergenic".format(intergenicsnptotal * 1000.0/intergenicsum) , file=sys.stderr)
	if excludesnps:
		print( "{} SNPs on excluded contigs".format(excludesnps) , file=sys.stderr)
	if args.histogram:
		sys.stderr.write( "# Writing coverage histogram to {}  {}\n".format(args.histogram, time.asctime() ) )
		with open(args.histogram, 'w') as histo:
			print( "Cov\tExons\tIntrons\tIntergenic" , file=histo )
			maxcov = max(chain(exoncovdict.keys(), introncovdict.keys(), igenecovdict.keys()))
			for i in range(maxcov):
				print( "{}\t{}\t{}\t{}".format(i+1, exoncovdict[i+1], introncovdict[i+1], igenecovdict[i+1]) , file=histo )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
