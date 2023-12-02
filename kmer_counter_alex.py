#!/usr/bin/env python

import sys
import subprocess
import itertools

k = 7
chr = 'chrY'
fastaFile = '%s.fa' % (chr)
kmerCmd = 'kmer-counter --fasta --no-rc --k=%d %s' % (k, fastaFile)

try:
    output = subprocess.check_output(kmerCmd, shell=True)
    result = {}
    for line in output.splitlines(): 
        (header, counts) = line.strip().split('\t')
        header = header[1:]
        kmers = dict((key,int(val)) for (key,val) in [d.split(':') for d in counts.split(' ')])
        result[header] = kmers
except subprocess.CalledProcessError as error:
    sys.stderr.write("%s\n" % (str(error)))

kmers = result[chr]
comp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
for kmerList in itertools.product('ACGT', repeat=k):
    kmerKey = ''.join(kmerList)
    kmerCompKey = ''.join(reversed([comp.get(b,b) for b in kmerList]))
    if kmerKey not in kmers and kmerCompKey not in kmers:
        kmers[kmerKey] = 0

for key, val in sorted(kmers.iteritems(), key=lambda (key,val):(val,key)):
    sys.stdout.write("%s\t%s\n" % (key, val))
