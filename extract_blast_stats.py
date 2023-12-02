#!/usr/bin/env python
#
# extract_blast_stats.py  created 2020-04-02
# v1 2021-12-26
# v1.1 2023-01-31 report subject counts

'''extract_blast_stats.py  v1 last modified 2023-01-31
    extract basic stats about blast results, number of hits, etc

extract_blast_stats.py -q query_prots.fasta -d ref_prots.fasta -b blast_output.tab

    produces tabular output of 7 columns:
qSeqid  qLength  numHits  bestHit  bitScore  sLength  numMatched

    qSeqid = query sequence ID/ FASTA header
    qLength = length of query, as either bases or amino acids
    numHits = total unique subject hits (counts by seq, not by HSPs)
    bestHit = name of highest scoring target
    bitScore = bitscore
    sLength = length of highest scoring target
    numMatched = number of queries that have that target as a hit

'''

import sys
import argparse
import gzip
from collections import defaultdict
from Bio import SeqIO


def get_fasta_seqs(fastafile):
	'''from a file name, read fasta sequences and return a dict where key is sequence ID'''
	if fastafile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading sequences from {} as gzipped\n".format(fastafile) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading sequences from {}\n".format(fastafile) )
	seqdict = SeqIO.to_dict(SeqIO.parse(opentype(fastafile,'rt'),"fasta"))
	sys.stderr.write("# Counted {} sequences from {}\n".format(len(seqdict), fastafile) )
	return seqdict

def parse_tabular_blast(blastfile, query_seqdict, sub_seqdict, is_verbose=False):
	'''from a file name, read blast results and return two dicts, one for each query and subject'''

	# this will therefore not include queries that do not have hits
	query_hits = defaultdict(list) # key is qseqid, value is list of all hits sseqid
	subject_hits = defaultdict(list) # key is sseqid, value is list of all hits qseqid

	best_hit_score = {} # key is seqid, value is summed bitscore of best match, including all HSPs
	best_hit_match = {} # key is seqid, value is seqid of best match

	if blastfile.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		opentype = gzip.open
		sys.stderr.write("# Reading blast results from {} as gzipped\n".format(blastfile) )
	else: # otherwise assume normal open for fasta format
		opentype = open
		sys.stderr.write("# Reading blast results annotations from {}\n".format(blastfile) )

	for line in opentype(blastfile,'rt'):
		line = line.strip()
		if line and line[0]!="#":
			lsplits = line.split("\t")
			qseqid = lsplits[0]
			sseqid = lsplits[1]
			bitscore = float(lsplits[11])

			query_hits[qseqid].append(sseqid)
			subject_hits[sseqid].append(qseqid)

			#TODO this should be changed to sum scores for all hits, not just the best hit
			if qseqid in best_hit_score:
				if sseqid == best_hit_match.get(qseqid,""): # sum bitscores from HPSs
					best_hit_score[qseqid] = best_hit_score.get(qseqid) + bitscore
				else:
					if bitscore > best_hit_score.get(qseqid):
						best_hit_score[qseqid] = bitscore
						best_hit_match[qseqid] = sseqid
					else:
						continue
			else: # meaning no entries yet
				best_hit_score[qseqid] = bitscore
				best_hit_match[qseqid] = sseqid

	sys.stderr.write("# {} of {} queries have hits\n".format( len(query_hits), len(query_seqdict) ) )

	headerline = "#qSeqid\tqLength\tnumHits\tbestHit\tbitScore\tsLength\tnumMatched\n"
	sys.stdout.write( headerline )

	# for queries with hits
	# print one tab delimited line for each query with its best hit
	for qseqid, hitlist in query_hits.items():
		best_score = best_hit_score.get(qseqid, 0.0)
		best_match = best_hit_match.get(qseqid, "None")
		num_hits = len(set(hitlist))
		if best_match != "None":
			best_match_length = len(sub_seqdict[best_match].seq)
		else:
			best_match_length = 0
		num_matched = len(set(subject_hits.get(best_match)))
		outline = "{}\t{}\t{}\t{}\t{:.2f}\t{}\t{}\n".format( qseqid, len(query_seqdict[qseqid].seq), num_hits, best_match, best_score, best_match_length, num_matched)
		sys.stdout.write( outline )

	# for queries without any hits
	for qseqid in query_seqdict.keys():
		if qseqid in query_hits:
			continue
		best_score = best_hit_score.get(qseqid, "None")
		best_match = best_hit_match.get(qseqid, "None")
		num_hits = 0
		best_match_length = 0
		num_matched = 0
		outline = "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format( qseqid, len(query_seqdict[qseqid].seq), num_hits, best_match, best_score, best_match_length, num_matched)
		sys.stdout.write( outline )

	sys.stderr.write("# {} of {} subjects were hit\n".format( len(subject_hits), len(sub_seqdict) ) )
	# no return


def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b','--blast-output', help="tabular blast output of query vs db")
	parser.add_argument('-d','--db', help="db nucleotides or proteins, as fasta")
	parser.add_argument('-q','--query', help="query nucleotides or proteins, as fasta")
	parser.add_argument('-r','--reference-blast', help="tabular blast output of query against some reference sequences")
	parser.add_argument('-v','--verbose', action="store_true", help="verbose output")
	args = parser.parse_args(argv)

	# read proteins from all genomes
	query_seqdict = get_fasta_seqs(args.query)
	sub_seqdict = get_fasta_seqs(args.db)

	# process blast hits
	parse_tabular_blast(args.blast_output, query_seqdict, sub_seqdict, args.verbose)


if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)

