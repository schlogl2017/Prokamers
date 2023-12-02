#! usr/bin/env python

from array import array
from collections import Counter

DNA_NUCLEOBASES = 'ACTG'

def numerify(s, k):
    '''
    Summarizes consecutive k-length substrings of s as as list of numbers.
    s is a string consisting of letters from DNA_NUCLEOBASES.
    k is an integer, 1 <= k <= 16.
    '''
    TRANS = dict(map(lambda c: (c, DNA_NUCLEOBASES.index(c)), DNA_NUCLEOBASES))
    mask = 4 ** k - 1
    # Depending on k, we might need a byte, short, or long.
    array_type = '!BBBBHHHHLLLLLLLL'[k]
    result = array(array_type, [0] * (len(s) - (k - 1)))
    v = 0
    for i in range(len(s)):
        v = (v << 2 & mask) | TRANS[s[i]]
        result[i - (k - 1)] = v
    return result.tolist()

def stringify(n, k):
    '''
    Converts a number from numerify() back into a k-length string of
    DNA_NUCLEOBASES.
    '''
    result = [DNA_NUCLEOBASES[(n >> 2 * i) & 3] for i in range(k)][::-1]
    return ''.join(result)

def find_clumps(s, k):
    '''
    Lists all k-length substrings of s in descending frequency.
    Returns a list of tuples of the substring and a list of positions at which
    they occur.
    '''
    clumps = numerify(s, k)
    result = []
    for clump, occurrences in Counter(clumps).most_common():
        if occurrences <= 1:
            break
        positions = []
        i = -1
        for _ in range(occurrences):
            i = clumps.index(clump, i + 1)
            positions.append(i)
        result.append((stringify(clump, k), positions))
    return result


s = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
k = 5
for clump, positions in find_clumps(s, k):
    print("%s occurs at %s" % (clump, positions))

