#!usr/bin/env python
# -*- coding: utf-8 -*-

import re
import os
import time
import gzip
import csv
import multiprocessing as mp
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import itertools
from itertools import product
from toolz import groupby
from more_itertools import windowed
from Bio import SeqIO
import alphabet


def get_gc_content(sequence, as_percent=False):
    """
    Finction to calculate the the gc content of a sequence.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.
    
    Outputs:
    
        gc - a float representing the of (g + c) content of a sequence.
    
    """
    seq_len = len(sequence)
    if not sequence:
        return 0.0
    # count all gs and cs
    c_plus_g = sequence.count('C') + sequence.count('G')
    # and divide by the sequence length
    gc = round(c_plus_g / seq_len, 4)
    if as_percent:
         return gc * 100
    return gc
