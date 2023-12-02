#!/usr/bin/env python
# coding: utf-8

# In[9]:


#!/usr/bin/env python
import os
import glob
import math
import numpy as np
import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from fasta_parser import parse_fasta
from sequence_utils import get_sequence_chunks, get_GC_by_slide_window,     get_gc_content, difference_gc, get_chromosomal_gc_variation


# In[11]:


fasta = "Genomes/Test/Chromosome/test_all_chr.fna"


# In[3]:


def get_gc_window_content(filename, window, step):
    gc_window = defaultdict(dict)
    for Id, sequence in parse_fasta(filename):
        geno_id = "_".join(Id.split("_")[:2])
        gc_window[geno_id] = get_GC_by_slide_window(sequence,
                                                get_gc_content,
                                                window,
                                                step)
    return gc_window


# In[6]:


gc_win = get_gc_window_content(fasta, 100, 100)
gc_window = pd.DataFrame(
    gc_win
    ).reset_index().rename(columns={'index': 'pos'})


# In[7]:


gc_window


# In[ ]:


plt.scatter(x=gc_window['pos'], y)
plt.show()


# In[8]:


def at_gc_data(filename):
    basic_data = defaultdict(list)
    for Id, seq in parse_fasta(filename):
        geno_id = "_".join(Id.split("_")[:2])
        gc = get_gc_content(seq)
        at = 1 - gc
        r = at/gc
        basic_data[geno_id] += [gc, at, r]
    return basic_data


# In[10]:


gc_at_r = at_gc_data(fasta)


# In[11]:


gc_at = pd.DataFrame(gc_at_r, 
                     index=['gc', 
                            'at', 
                            'ratio']).reset_index().rename(
    columns={'index': 'stats'})


# In[12]:


gc_at


# In[14]:


y = gc_at[gc_at['stats'] == 'gc'].drop('stats', axis=1).values
x = list(range(3))


# In[17]:


plt.hist(y)
plt.show()


# In[18]:


def get_chr_gc_dif_content(filename, window, step):
    gc_diff = defaultdict(dict)
    for Id, sequence in parse_fasta(filename):
        geno_id = "_".join(Id.split("_")[:2])
        gc = get_gc_content(sequence)
        gc_win = get_GC_by_slide_window(sequence,
                                        get_gc_content,
                                        window,
                                        step)
        gc_diff[geno_id] = difference_gc(gc, gc_win)
    return gc_diff


# In[20]:


gc_chr_dif = get_chr_gc_dif_content(fasta, 100, 100)


# In[21]:


gc_dif = pd.DataFrame(gc_chr_dif).reset_index().rename(columns={'index': 
                                                                'pos'})


# In[22]:


gc_dif['NC_038132.1']


# In[23]:


gc_dif


# In[24]:


def get_chromosomal_gc_variation(difference_dic):
    """
    Calculates the chromosomal GC variation defined as the log-transformed average
    of the absolute value of the difference between the mean GC content of each
    non-overlapping sliding window i and mean chromosomal GC content.
    chromosomal_gc_variantion = log(1/N*sum(|Di|), where N is the maximum number of
    non-overlapping window size in bp in a sliding windows.
    """
    # get the number of chunks
    n = len(difference_dic)
    # Create a new 1-dimensional array from an iterable object
    arr = np.fromiter(
        difference_dic.values(),
        dtype=np.float32,
        count=len(difference_dic),
    )
    var = np.log(np.sum(np.abs(arr)) / n)
    return var


# In[175]:


c = {0:-0.0394,100:-0.0294, 200: -0.0594, 400: -0.0294, 600:0.0106}


# In[25]:


def get_gc_chr_var(filename, window, step):
    gcVar = defaultdict(float)
    for Id, sequence in parse_fasta(filename):
        geno_id = "_".join(Id.split("_")[:2])
        gc = get_gc_content(sequence)
        gc_win = get_GC_by_slide_window(sequence,
                                        get_gc_content,
                                        window,
                                        step)
        gc_dif = difference_gc(gc, gc_win)
        gcVar[geno_id] = get_chromosomal_gc_variation(gc_dif)
    return gcVar


# In[26]:


gc_dif.head()


# In[151]:


gc_dif["NC_038132.1"]


# In[154]:


sns.histplot(data=gc_dif, x="NC_011943.1", kde=True)
plt.xlim((-0.4, 0.4))
plt.show()


# In[45]:


gcvar = get_gc_chr_var(filename, 100, 100)


# In[48]:


pd.DataFrame(gcvar.items(), columns=['geno_id', 'gc_var'])


# In[49]:


def get_chromosomal_gc_variation(difference_dic):
    """
    Calculates the chromosomal GC variation defined as the log-transformed average
    of the absolute value of the difference between the mean GC content of each
    non-overlapping sliding window i and mean chromosomal GC content.
    chromosomal_gc_variantion = log(1/N*sum(|Di|), where N is the maximum number of
    non-overlapping window size in bp in a sliding windows.
    """
    # get the number of chunks
    n = len(difference_dic)
    # Create a new 1-dimensional array from an iterable object
    arr = np.fromiter(
        difference_dic.values(),
        dtype=np.float32,
        count=len(difference_dic),
    )
    var = np.log(np.sum(np.abs(arr)) / n)
    return var


# In[51]:


get_chromosomal_gc_variation(gc_chr_dif['NC_038132.1'])


# In[55]:


gcvar


# In[30]:


GCVAR = get_chromosomal_gc_variation(gc_dif)


# In[34]:


sns.histplot(data=gcvar, x="gc_var", kde=True)
plt.xlim((-0.4, 0.4))
plt.show()


# In[33]:


sns.kdeplot(data=gcvar, x="gc_var")
plt.xlim((-0.4, 0.4))
plt.show()


# In[49]:


dd = {'name':{'a': 3.0, 'b':1.0, 'c': 0.34}}


# In[50]:


dd


# In[44]:


pd.DataFrame(d, index=["gc", "at", "at_gc_ratio"]) #.T.rename(columns={"gc":0})


# In[79]:


def check_paths(full_path_list):
    for p in full_path_list:
        if os.path.isfile(p):
            pass
        elif not os.path.exists(p):
            os.makedirs(p)


# In[81]:


check_paths(['Genomes/Archaea', 'Genomes/tests_empty/text.txt'])


# In[53]:


ss = dd
ss


# In[55]:


ddd = defaultdict(dict)


# In[56]:


ddd.update(ss)


# In[57]:


ddd


# In[29]:


len(s1)


# 05/10/2022

# In[21]:


filename="Genomic_data/Archaea/Asgard/GCA_008000775.1/Chromosome/GCA_008000775.1_Candidatus_Prometheoarchaeum_chr.fna"


# In[22]:


for name, sequence in parse_fasta(filename):
    seq = sequence
    name = name


# In[ ]:





# In[ ]:





# In[ ]:





# In[10]:


# def gc_content(sequence):
#     """
#     Finction to calculate the the gc content of a sequence.
    
#     Inputs:
    
#         sequence - a string representing a DNA sequence.
    
#     Outputs:
    
#         gc - a float representing the of (g + c) content of a sequence.
    
#     """
#     # get the sequence length and 
#     # make all the sequence characters upper case
#     seq_len = len(sequence)
#     all_bases = np.array([sequence.count(i) for i in 'ACGT'])
#     cg = np.sum(np.array([sequence.count(i) for i in 'CG']))

#     return dict(zip('ACGT', all_bases)), round((cg / seq_len) * 100, 4)


# In[23]:


gc_content(seq)


# In[24]:


get_ipython().run_line_magic('timeit', 'gc_content(seq)')


# In[11]:


def gc_content(seq, percent = False):
    """
    Calculate G+C content, return percentage (as float between 0 and 100).
    Copes mixed case sequences, and with the ambiguous nucleotide S (G or C)
    when counting the G and C content.  The percentage is calculated against
    the full length.

    Inputs:
        sequence - a string representing a sequence (DNA/RNA/Protein)

    Outputs:
        gc - a float number or a percent value representing the gc content of the sequence.

    """
    seq = seq.upper()
    seq_len = len(seq)
    g_c = seq.count("G") + seq.count("C")
    gc = g_c / seq_len
    if percent:
        return gc * 100
    return gc


# In[25]:


gc_content(seq, percent = False)


# In[26]:


get_ipython().run_line_magic('timeit', 'gc_content(seq, percent = False)')


# In[12]:


def get_gc_content_slide_window(sequence, window, step, percent =False):
    """
    Calculate the GC (G+C) content along a sequence. Returns a dictionary-like
    object mapping the sequence window to the gc value (floats).
    Returns 0 for windows without any G/C by handling zero division errors, and
    does NOT look at any ambiguous nucleotides.

    Inputs:
        sequence - a string representing a sequence (DNA/RNA/Protein)
        window_size - a integer representing the length of sequence to be analyzed
                   at each step in the full sequence.
        step - a integer representing the length of overlapping sequence allowed.

    Outputs:
        gc- a dictionary-like object mapping the window size sequence to it gc values
            as floating numbers
    """
    gc = defaultdict(float)
    for start, seq in get_sequence_chunks(sequence, window, step):
        gc[start] += gc_content(seq, percent = percent)
    return gc


# In[27]:


get_gc_content_slide_window(sequence, 2000, 2000, percent =False)


# In[28]:


# non overlapping chunks
get_ipython().run_line_magic('timeit', 'get_gc_content_slide_window(sequence, 2000, 2000, percent =False)')


# In[29]:


# overlapping chunks
get_ipython().run_line_magic('timeit', 'get_gc_content_slide_window(sequence, 2000, 500, percent =False)')


# In[13]:


def difference_gc(total_gc, gc_dict):
    """
    Calculates the difference between the mean GC content of window i, and the
    mean chromosomal GC content, as Di = GC i − GC
    """
    # iterate through all keys, gc calculated values for the
    # genome chunk
    # get the difference between the chromosomal mean and the chunks
    # add the difference to the appropriate chunk
    d_i = {
        chunk: gc_cnt - total_gc
        for chunk, gc_cnt in gc_dict.items()
    }
    return d_i


# In[30]:


gctot = gc_content(seq, percent = False)
gc_dict = get_gc_content_slide_window(sequence, 2000, 2000, percent =False)


# In[31]:


difference_gc(gctot, gc_dict)


# In[32]:


get_ipython().run_line_magic('timeit', 'difference_gc(gctot, gc_dict)')


# In[33]:


diff_gc = difference_gc(gctot, gc_dict)


# In[ ]:





# In[61]:


sns.histplot(data=diff_gc,kde=True)
plt.xlim((-0.4, 0.4))
plt.show()


# In[34]:


def get_chromosomal_gc_variation(difference_dic):
    """
    Calculates the chromosomal GC variation defined as the log-transformed average
    of the absolute value of the difference between the mean GC content of each
    non-overlapping sliding window i and mean chromosomal GC content.
    chromosomal_gc_variantion = log(1/N*sum(|Di|), where N is the maximum number of
    non-overlapping window size in bp in a sliding windows.
    """
    # get the number of chunks
    n = len(difference_dic)
    arr = np.fromiter(
        difference_dic.values(),
        dtype=np.float32,
        count=len(difference_dic),
    )
    var = np.log(np.sum(np.abs(arr)) / n)
    return var


# In[35]:


get_chromosomal_gc_variation(diff_gc)


# In[36]:


get_ipython().run_line_magic('timeit', 'get_chromosomal_gc_variation(diff_gc)')


# In[15]:


def get_sequence_skew(sequence):
    """
    Calculates the difference between the total number of
    occurrences of G and the total number of occurrences of C in
    the first i elements of the sequence. 

    Inputs:
        sequence - string representing the sequence     
    
    Outputs:
    
        skew - an array-like object that represents the GC skew of the sequence.
    
    Ex: 
    > get_sequence_skew('ACAACGTAGCAGTAGCAGTAGT')
    [0, 0, -1, -1, -1, -2, -1, -1, -1, 0, -1, -1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 2, 2]
    
    """
    # make the sequence upper case
    sequence = sequence.upper()
    # start the array
    skew = [0]
    # iterates to the sequence elements and it indexes
    for idx, element in enumerate(sequence):
        # check if element[i] is a G
        # if so add 1
        if sequence[idx] == 'G':
            skew.append(skew[idx] + 1)
        # if the element[i] is a C
        # add to the array -1
        elif sequence[idx] == 'C':
            skew.append(skew[idx] - 1)
        else:
            # if it is not G or C add 0
            skew.append(skew[idx])
    return skew


# In[37]:


get_sequence_skew(seq)


# In[38]:


get_ipython().run_line_magic('timeit', 'get_sequence_skew(seq)')


# In[39]:


skgc = get_sequence_skew(seq)


# In[16]:


# def get_minimum_skew(sequence):
#     """
#     Calculates a position in a sequence minimizing the skew.
    
#     Inputs:
#         sequence - string representing the sequence.    
    
#     Outputs:
    
#         min_skew - an array-like object that represents the pistion 
#                    where the GC skew is the minimized in the sequence.    
    
#     Example:
#     get_minimum_skew('ACAACGTAGCAGTAGCAGTAGT')
#     [5]
#     """
#     min_skew = []
#     # calculates the sequence gc skew
#     skew = get_sequence_skew(sequence)
#     m_skew = min(skew)
#     # to get the index positions
#     for idx in range(len(sequence) + 1):
#         # if the position i has the same value 
#         # as the minimum appende to the array
#         if skew[idx] == m_skew:
#             min_skew.append(idx)
#     return min_skew


# In[40]:


get_minimum_skew(seq)


# In[41]:


get_ipython().run_line_magic('timeit', 'get_minimum_skew(seq)')


# In[17]:


def get_min_skew(seq_skew):
    """
    Calculates a position in a sequence minimizing the skew.
    
    Inputs:
        seq_skew - a list representing the sequence skew.    
    
    Outputs:
    
        min_skew - an array-like object that represents the pistion 
                   where the GC skew is the minimized in the sequence.    
    
    Example:
    get_minimum_skew('ACAACGTAGCAGTAGCAGTAGT')
    [5]
    """
    min_skew = []
    skew = seq_skew
    # get the minimized skew values
    m_skew = min(skew)
    # iterates to the length of the sequence
    # to get the index positions
    for idx in range(len(sequence) + 1):
        # if the position i has the same value 
        # as the minimum appende to the array
        if skew[idx] == m_skew:
            min_skew.append(idx)
    return min_skew


# In[42]:


get_min_skew(skgc)


# In[43]:


get_ipython().run_line_magic('timeit', 'get_min_skew(skgc)')


# In[44]:


minskwgc = get_min_skew(skgc)


# In[45]:


def get_gc_strand_difference(sequence, start):
    seq = sequence.upper()
    seq_len = len(seq)
    half = seq_len // 2
    ter = start + half
    g, c = 0, 0
    if ter > seq_len:
        ter = ter - seq_len + 1
    elif ter > start:
        g += 2 * seq[start:ter].count('G') - seq.count('G')
        c += 2 * seq[start:ter].count('C') - seq.count('C')
    else:
        g += seq.count('G') - 2 * seq[start:ter].count('G')
        c += seq.count('C') - 2 * seq[start:ter].count('C')
    return g - c


# In[46]:


get_gc_strand_difference(seq, 0)


# In[48]:


get_gc_strand_difference(seq, 36)


# In[57]:


# def kmer_map(sequence, k):
#     """
#     Function to find the kmer positions in a sequence.
    
#     Inputs:
#         sequence - a string representing a sequence
#         k - a integer representing the length of the
#             substrings.
    
#     Outputs:
#         kmermap - a dictionary-like object mapping the
#                   substrings of k length and they positions
#                   in the sequence.
#     Ex:
#     >> kmer_map('ATTGATTATTG', 3)
#     defaultdict(list,
#             {'ATT': [0, 4, 7],
#              'TTG': [1, 8],
#              'TGA': [2],
#              'GAT': [3],
#              'TTA': [5],
#              'TAT': [6]})
#     """
#     # get the length through iterate over
#     seq_len = len(sequence) - k + 1
#     # get the container fro the data
#     kmermap = defaultdict(list)
#     # iterates through the sequence length
#     for i in range(seq_len):
#         # get the kmer as keys and positions
#         # as values
#         kmermap[sequence[i:i+k]].append(i)
#     return kmermap    


# In[58]:


kmer_map(seq, 4)


# In[59]:


get_ipython().run_line_magic('timeit', 'kmer_map(seq, 4)')


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[97]:


def plot_slide_window_data(data_dic, x_label, y_label, title, as_save=False):
    """
    Plot the data calculated in a slide window of a given sequence.
    
    Inputs:
        data_dic - a dictionary-like object representing the data calculated in a
                   sequence in a window size sample.
        x_label - a string representing the x labels to x axis of the plot.
        y_label - a string representing the y labels to y axis of the plot.
        title - a string representing the title of the plot.
        as_save - a boolean that is default to false and not save any figure.
    
    Outputs:
        a plot figures to be save is as_save or just show the plot.
    """
    plt.figure(num=None, figsize=(15, 7), dpi=100)
    plt.plot(data_dic.keys(), data_dic.values(), color='#1f77b4', linestyle='--', marker='.', alpha=0.5)
    plt.tight_layout()
    plt.xlabel(f'{x_label} (kbp)')
    plt.ylabel(f'{y_label}')
    plt.title(f'{title}')
    plt.xticks(rotation=90)
    plt.grid()
    if as_save:
        plt.savefig(f'{title}.png')
    else:
        plt.show()


# In[102]:


def test() -> None:
    s = (
        "TAGTTGTGAAGAAAATATGGATAAACAGGACGACGAATGCTTTCACCGATAAGGACAACTTTCCATAACTCAGTAAATATAGTGCAGAGTTCACCCTGTC"
    )

    seq_gc_non_overlap = get_gc_content_slide_window(s, 20, 20, True)
    expected = {
        0: 30,
        20: 45,
        40: 40,
        60: 25,
        80: 55,
    }
    assert expected.keys() == seq_gc_non_overlap.keys()
    for k, v in expected.items():
        assert math.isclose(seq_gc_non_overlap[k], v)

    seq_gc_total = gc_content(s, percent=True)
    assert math.isclose(seq_gc_total, 39)

    seq_dif_gctot_gc_slidewindow = difference_gc(seq_gc_total, seq_gc_non_overlap)
    print(seq_dif_gctot_gc_slidewindow)
    expected = {
        0: -9,
        20: 6,
        40: 1,
        60: -14,
        80: 16,
    }
    assert expected.keys() == seq_dif_gctot_gc_slidewindow.keys()
    for k, v in expected.items():
        assert math.isclose(seq_dif_gctot_gc_slidewindow[k], v)

    chromosome_gc_variation = get_chromosomal_gc_variation(seq_dif_gctot_gc_slidewindow)
    assert math.isclose(chromosome_gc_variation, 2.2192034840549946)


# In[103]:


test()


# In[27]:


fasta_files = glob.glob('Data/Genomes_splitted/*/Chromosomes/*_chr.fna.gz')


# In[112]:


def get_gc_data_from_genomes(filename, window, step, percent):
    gc_di = defaultdict(dict)
    gcvar = defaultdict(dict)
    for Id, seq in parse_fasta(filename):
        geno_id = "_".join(Id.split("_")[:2])
        gc_chr = gc_content(seq)
        gc_pb = get_gc_content_slide_window(seq, window, step, percent)
        gc_diff = difference_gc(gc_chr, gc_pb)
        gc_var = get_chromosomal_gc_variation(gc_diff)
        gc_di[geno_id] = gc_diff
        gcvar[geno_id] = gc_var
    return gc_di, gcvar    


# In[164]:


di, var = get_gc_data_from_genomes(fasta, 100, 100, percent = True)


# In[189]:


di


# In[165]:


dif = pd.DataFrame(di).reset_index().rename(columns={'index': 'pos'})


# In[196]:


ndif = dif.drop(columns=['pos'])


# In[198]:


ndif.plot()


# In[168]:


x = dif['pos'].values


# In[169]:


y = dif['NC_038132.1'].values


# In[ ]:





# In[170]:


g =sns.scatterplot(x=x, y=y);


# In[ ]:





# In[199]:


sns.scatterplot(x=x, y=dif['NC_021932.1'].values);


# In[ ]:





# In[171]:


sns.histplot(data=dif, x="NC_038132.1", kde=True)
plt.show()


# In[172]:


sns.kdeplot(data=dif, x="NC_038132.1")
plt.show()


# In[156]:


def gc_content_sequence_window(sequence, as_overlap=False, k=20):
    """GC Content in a DNA/RNA sub-sequence length k. In
    overlapp windows of lenght k.
    
    Inputs:
    
        sequence - a string representing a DNA sequence.    
        as_overlap - boolean that represents if overlap is needed.
        k - a integer reppresenting the lengths of overlappig bases.
            Default is 20.
    
    Outputs:
    
        gc_content - an array-like object with 
    
    
    """
    # make sequence upper case and getting the length of it
    sequence, seq_len = sequence.upper(), len(sequence)
    # the array-like object to collect the data
    gc_content = []
    # non overlap sequence length
    non_overlap = range(0, len(sequence) - k + 1, k)
    # overlap sequence length
    overlap = range(0, seq_len - k + 1)
    # overlap is needed
    if as_overlap:
        # iterates to the overlap region
        for i in overlap:
            # creates the substring to count the gc_content
            subseq = sequence[i:i + k]
            # count and sum up the Gs and Cs counts
            g_c = subseq.count('C') + subseq.count('G')
            # collect the data in the array container
            gc_content.append(round(g_c / len(subseq), 4) * 100)
    # if non overlap is choosed
    else:
        # iterates to the mon overlap region
        for j in non_overlap:
            # creates the substring to count the gc_content
            subseq = sequence[j:j + k]
            # count and sum up the Gs and Cs counts
            g_c = subseq.count('C') + subseq.count('G')
            # collect the data in the array container
            gc_content.append(round(g_c / len(subseq), 4) * 100)
    return gc_content


# In[160]:


for Id, seq in parse_fasta("Genomes/Test/GCF_000019705.1_ASM1970v1_genomic.fna.gz"):
    seq = seq


# In[162]:


gcw = gc_content_sequence_window(seq, as_overlap=False, k=100)


# In[157]:


def gc_var(sequence, as_overlap=False, k=20):
    """
    Calculates the gc content variance in a sequence according to a 
    window of length k.
    
    Inputs:
        sequence - a string representing a DNA sequence.
        k - integer representing the length of the search window.
            default is 20.
    
    Outputs:
    
        log of the gc variantion in the sequence in a window space of
        length k.
    
    """
    # calculates the percent of gc content
    gc = get_gc_content(sequence) * 100
    # get the gc content in the window space as an array
    gc_i = np.array(gc_content_sequence_window(sequence, as_overlap, k=k))
    # get the len of the gc content in the window space
    len_gc_i = np.shape(gc_i)[0]
    # check the difference of each point 
    dif = gc_i - gc
    return np.log((1 / len_gc_i) * sum(abs(dif)))


# In[163]:


gc_var(seq, as_overlap=False, k=100)


# In[176]:


def get_chromosomal_gc_variation(difference_dic):
    """
    Calculates the chromosomal GC variation defined as the log-transformed average
    of the absolute value of the difference between the mean GC content of each
    non-overlapping sliding window i and mean chromosomal GC content.
    chromosomal_gc_variantion = log(1/N*sum(|Di|), where N is the maximum number of
    non-overlapping window size in bp in a sliding windows.
    """
    # get the number of chunks
    n = len(difference_dic)
    # Create a new 1-dimensional array from an iterable object
    arr = np.fromiter(
        difference_dic.values(),
        dtype=np.float32,
        count=len(difference_dic),
    )
    print()
    var = np.log(np.sum(np.abs(arr)) / n)
    return var


# In[ ]:


c = {0:-0.0394,100:-0.0294, 200: -0.0594, 400: -0.0294, 600:0.0106}


# In[177]:


get_chromosomal_gc_variation(c)


# In[179]:


arr = np.fromiter(
        c.values(),
        dtype=np.float32,
        count=len(c),
    )


# In[180]:


arr


# In[182]:


np.sum(np.abs(arr)) / len(c)


# In[183]:


np.log(np.sum(np.abs(arr)) / len(c))


# In[188]:


math.log(0.03364000022411347)


# In[184]:


ar = arr * 100


# In[185]:


ar


# In[186]:


np.sum(np.abs(ar)) / len(c)


# In[187]:


np.log(np.sum(np.abs(ar)) / len(c))


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




