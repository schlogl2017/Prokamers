#!/usr/bin/env python
#
# prints output file name for each file 2018-09-06
# v1.01 some notes on usage 2014-02-12
# degapper.py v1 2014-02-07
#
# script to remove gaps from fasta sequences (alignments)

'''degapper.py  last modified 2023-02-08
remove gaps from many sequence files or alignments

degapper.py proteins1.fasta proteins2.fasta ...
degapper.py *.fasta

generates output files as proteins1.fasta.nogaps ...
'''

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def main(argv, wayout):
	if len(argv) < 2:
		sys.exit(__doc__)
	else:
		exclusion_list = []
		# for example:
		#exclusion_list = ["Pleurobrachia_bachei_AUG"]
		if exclusion_list:
			for exclude_name in exclusion_list:
				sys.stderr.write( "# Excluding taxon: {}\n".format( exclude_name ) )
		filelist = argv[1:]
		filecount = 0
		seqcount = 0
		for fastafile in filelist:
			filecount += 1
			outputfilename = str(fastafile)+".nogaps"
			sys.stderr.write( "# Removing gaps from {}, writing to {}\n".format( fastafile, outputfilename ) )
			with open(outputfilename,'w') as nogapfile:
				for seqrec in SeqIO.parse(fastafile,'fasta'):
					if seqrec.id in exclusion_list: # also remove some taxa
						continue
					seqcount += 1
					gappedseq = str(seqrec.seq)
					degappedseq = Seq(gappedseq.replace("-","").replace("X",""))
					seqrec.seq = degappedseq
					nogapfile.write("%s" % seqrec.format("fasta"))

		sys.stderr.write( "# read {} files\n".format(filecount) )
		sys.stderr.write( "# wrote {} sequences\n".format(seqcount) )


if __name__ == "__main__":
	main(sys.argv, sys.stdout)
