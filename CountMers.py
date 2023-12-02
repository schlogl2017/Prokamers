import re, sys
from Bio import SeqIO
import random
import time

k = 9
n = 10
kmers = re.compile("(?=(\w{%s}))" % k)

def timer(start, end):
    hours, rem = divmod(end-start, 3600)
    minutes, seconds = divmod(rem, 60)
    return hours, minutes, seconds

start = time.time()

for rec in SeqIO.parse(sys.argv[1], 'fasta'):
    s = time.time()
    kseqs = kmers.findall(str(rec.seq))
    randomsel = random.sample(kseqs, n)
    e = time.time()
    h, m, s = timer(s, e)
    print("10 random kmers from {} [Elapsed time: {:0>2}:{:0>2}:{:05.2f}]\n{}".formatrec.id, int(h), int(m), s, str(randomsel)))

end = time.time()
hours, minutes, seconds = timer(start, end)
