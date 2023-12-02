#!/usr/bin/python

try:
    infile = open('input.txt', 'r')
except IOError, e:
    raise SystemExit, "Hey! Where's the input file?"

A = infile.readline().strip()
B = infile.readline().strip()
w = int(infile.readline().strip())
k = int(infile.readline().strip())
infile.close()

m = len(A)
n = len(B)

# Create a (m - w + 1)x(n - w + 1) matrix initialized with 0s.
dotplot = [[0 for j in range(n - w + 1)] for i in range(m - w + 1)]

# For each position (i,j) in the matrix, count the
# number of matches in a window of length w.
for i in range(m - w + 1):
    for j in range(n - w + 1):
        count = 0
        for l in range(w):
            if A[i + l] == B[j + l]: count += 1
        dotplot[i][j] = (count >= k and 1 or 0)

# Print the result.
for row in dotplot:
    print ' '.join(map(str, row))

