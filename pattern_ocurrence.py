#! usr/bin/env python
from collections import defaultdict
from fasta_parser import parse_fasta

def pattern_occurrences(pattern, text):
    """Generate the indices of the (possibly overlapping) occurrences of
    pattern in text. For example:

    >>> list(pattern_occurrences('ATAT', 'GATATATGCATATACTT'))
    [1, 3, 9]

    """
    position = -1
    while True:
        position = text.find(pattern, position + 1)
        if position == -1:
            return
        yield position


# find substrings at least kmer long that occur at least min_clumpsize
# times in a window of windowsize in genome

def get_substrings(g, k):
    """
Take the input genome window 'g', and produce a list of unique 
substrings of length 'k' contained within it. 
    """
    substrings = list()

    # Start from first character, split into 'k' size chunks
    # Move along one character and repeat. No sense carrying on beyond
    # a starting point of 'k' since that will be the first iteration again.
    for i in range(k):
        line = g[i:]
        substrings += [line[i:i + k]
                       for i in range(0, len(line), k) if i + k <= len(line)]

    # Using collections.Counter increases the runtime by about 3 seconds,
    # during testing.
    results = defaultdict(int)
    for s in substrings:
        results[s] += 1
    return results


def find_clumps(genome, kmer, windowsize, clumpsize):
    """
In a given genome, examines each windowsize section for strings of length kmer
that occur at least clumpsize times. 

Input: 
genome: text string to search
kmer:  length of string to search for
windowsize: size of the genome section to consider for clumping
clumpsize: the kmer length strings must occur at least this many times

Returns: a list of the strings that clump
    """
    window = genome[0:windowsize]

    # Initialise our counter, because the main algorithm can't 
    # start from scratch.
    patterns = get_substrings(window, kmer)

    # Using a dictionary not a list because the lookups are faster
    # once the size of the object becomes large
    relevant = {p: 1 for p in patterns if patterns[p] >= clumpsize}

    starting_string = genome[0:kmer]

    for i in range(windowsize, len(genome)):
        # Move the window along one character
        window = window[1:]
        window += genome[i]

        # This is the only string that can decrease if we've moved
        # one character
        patterns[starting_string] -= 1
        starting_string = window[0:kmer]

        # This is the only string that can increase if we've moved
        # one character
        ending_string = window[-kmer:]
        patterns[ending_string] += 1

        # if there are enough matches of the string at the end, 
        # add it to matches.
        if patterns[ending_string] >= clumpsize and ending_string not in relevant:
            relevant[ending_string] = 1
    return list(relevant)


if __name__ == "__main__":
    for Id, seq in parse_fasta('Genomes_new/Archaea/Asgard_group/C_Heimdallarchaeota/CHR/GCA_020348965.1_chr.fna'):
        genome = seq
    
    kmer = 9  # Length of the k-mer
    windowsize = 500  # genome substring length to register clumps in
    min_clumpsize = 3  # minimum number of repetitions of the k-mer
    clumps = find_clumps(genome, kmer, windowsize, min_clumpsize)
    print("Total: {}".format(len(clumps)))
    print(clumps)


