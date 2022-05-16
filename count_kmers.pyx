#!usr/bin/env python

from itertools import product
from collections import defaultdict
from alphabet import iupac_dna


def get_all_possible_kmers(alphabet, kmin, kmax):
    """Returns a list of all possible combinations of k-mers of
    length k from a input alphabet.

    Inputs:

        alphabet - a alphabet (strings characters) that compound the string sequence
        kmin - minimum DNA kmer length (int)
        kmax - maximum DNA kmer length (int)

    Outputs:

        kmers - list of all possible combinations of k-mers of length k with length
                between kmin and kmax.

    """
    kmers = [''.join(letters) for n in range(kmin, kmax + 1)
             for letters in product(alphabet, repeat=n)]
    return kmers


def count_kmers(str string, int kmin, int kmax, counter=None):
    """
    Count occurrence of kmers in a given string.
    """
    if counter is None:
        counter = defaultdict(int, [(km, 0) for km in get_all_possible_kmers(iupac_dna,
                                                                             kmin,
                                                                             kmax)])
    cdef int i
    cdef int j
    cdef int N = len(string)
    for k in range(kmin, kmax + 1):
        for i in range(N - k + 1):
            kmer = string[i:i+k]
            counter[kmer] = counter.get(kmer, 0) + 1
    return counter
