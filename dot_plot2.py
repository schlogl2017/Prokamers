#!usr/bin/env python



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


def plot_matrix(matrix, t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1, matrix):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)


def dotplot(seq1, seq2, nonblank, blank, k = 1, t = 1):
    matrix = make_matrix(seq1, seq2, k)
    #experiment with character choice
    plot_matrix(matrix, t, seq1, seq2, nonblank, blank)


def test():
    seqx = "ACCTGAGCTCACCTGAGTTA"
    seqy = "ACCTGAGCTCACCTGAGTTA"
    dotplot(seqx, seqy, nonblank =".", blank=" ")

test()









