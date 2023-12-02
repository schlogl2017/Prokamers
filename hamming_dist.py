#!usr/bin/env python

def hamming_distance(seq1, seq2):
    assert(len(seq1) == len(seq2))
    hd = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 != c2:
            hd += 1
    return hd

def test():
    hd = hamming_distance("ACTTTGTT", "AGTTTCTT")
    assert(hd == 2)
    print(hd)

test()
