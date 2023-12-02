#!usr/bin/env python

import numpy
import matplotlib.pyplot as plt


def delta(x, y):
    return 0 if x == y else 1


def matrix(seq1, seq2, i, j, k):
    return sum(delta(x, y) for x, y in zip(seq1[i:i+k], seq2[j:j+k]))


def make_matrix(seq1, seq2, k):
    n = len(seq1)
    m = len(seq2)
    return [[matrix(seq1, 
                    seq2, 
                    i, 
                    j, 
                    k) for j in range(m - k + 1)] for i in range(n -k + 1)]


seqx = "ACCTGAGCTCACCTGAGTTA"
seqy = "ACCTGAGCTCACCTGAGTTA"


dotplot = plt.imshow(numpy.array(make_matrix(seqx, seqy, 1)))

xt = plt.xticks(numpy.arange(len(list(seqx))),list(seqx))
yt = plt.yticks(numpy.arange(len(list(seqx))),list(seqx))
plt.show()
