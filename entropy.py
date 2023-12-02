#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 
import collections
import math
 
 
def estimate_shannon_entropy(dna_sequence):
    m = len(dna_sequence)
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
 
    shannon_entropy_value = 0
    for base in bases:
        # number of residues
        n_i = bases[base]
        # n_i (# residues type i) / M (# residues in column)
        p_i = n_i / float(m)
        entropy_i = p_i * (math.log(p_i, 2))
        shannon_entropy_value += entropy_i
 
    return shannon_entropy_value * -1

def calcRelativeEntropy(seq, resCodes):
    """Calculate a relative entropy value for the residues in a
    sequence compared to a uniform null hypothesis."""
    N = float(len(seq))
    base = 1.0/len(resCodes)
    prop = {}
    for r in resCodes:
        prop[r] = 0
    for r in seq:
        prop[r] += 1
    for r in resCodes:
        prop[r] /= N
    H = 0 
    for r in resCodes:
        if prop[r] != 0.0:
            h = prop[r]* math.log(prop[r]/base, 2.0)
            H += h
    H /= math.log(base, 2.0)
    return H



import collections
 
from scipy.stats import entropy
  
def estimate_shannon_entropy2(dna_sequence):
    bases = collections.Counter([tmp_base for tmp_base in dna_sequence])
    # define distribution
    dist = [x/sum(bases.values()) for x in bases.values()]
 
    # use scipy to calculate entropy
    entropy_value = entropy(dist, base=2)
 
    return entropy_value
