#!/usr/bin/env python
import os
import math
import glob
from time import time
import argparse
from termcolor import colored
from system_utils import make_me_a_folder
import fasta_parser
import count_kmers
from sequence_utils import base_stats, get_all_possible_kmers


def signif(x, digits=6):
    """ 
    It rounds the values in its first argument to the specified number 
    of significant digits. 
    """
    # test if x is a positive number and is finite
    # if not
    if x == 0 or not math.isfinite(x):
        return x
    # if yes
    digits -= math.ceil(math.log10(abs(x)))
    return round(x, digits)


def genome_length(filename, name):
    genome_len = pd.read_csv(filename, names=['name', 'len'])
    gen_len = genome_len['len'].loc[genome_len['name'] == name].values[0]
    return gen_len


def get_kmer_counts(sequence, k, as_freq=False):
    # get kmer counts
    kmc = count_kmers.count_kmers(sequence, k, k)
    total = len(sequence) - k + 1
    # if freq is needed
    if as_freq:
        return {k: v/total for k, v in kmc.items()}
    return kmc    
    
    
def kmer_expected_zero_markov(base_freqs,  kmer_list, seq_len):
    """
    Calculates the 0-order markov given the background nucleotide distribution.
    
    Inputs:
        base_freqs - a dictionary-like obeject mapping the nucleotides to its calculated
                     frequence in a given sequence.
        kmer_list - a list of substrings of length k
        seq_len - a integer representing the length of the original sequence where the
                  calculated bases frequencies came from.
    
    Outputs:
        expected - a dictionary-like obeject mapping the expected kmer values
                   calculated as, (ex = math.pow(A, a)*math.pow(C, c)*math.pow(G, g)\
                   *math.pow(T, t)*(seq_len-k+1), where A,C,G,T are the frequency of each 
                   nucleotide in the genome, a,c,g,t are the number of each nucleotide in the k-mer,
                   and seq_len is the length of the kmer sequence.
    """
    # get kmer length
    k = len(kmer_list[0])
    # create the countainer and it keys from a kmer-list
    expected = defaultdict(float, [(km, 0.0) for km in kmer_list])
    # initiate the variables to receive the base frequencies
    A, C, G, T = base_freqs['A'], base_freqs['C'], base_freqs['G'], base_freqs['T']
    # iterate through each kmer in the list
    for kmer in kmer_list:
        # count the number of bases tha compound the kmer sequence
        a, c, g, t = kmer.count('A'), kmer.count('C'), kmer.count('G'), kmer.count('T')
        # calculates the expected values for taht kmer in the sequence accord to 
        # its base composition and add the value to the respectives keys in the dictionary
        ex = math.pow(A, a) * math.pow(C, c) * math.pow(G, g) * math.pow(T, t) * (seq_len - k + 1)
        expected[kmer] = signif(ex)
    # return the results
    return expected    
    
    
def get_oligo_nucleotide_usage_deviation_normalized(kmer_count, expected_zero):
    """
    The normalized value for a word W is calculated by dividing the
    observed counts by the expected counts. This is the usage deviation
    vector for a genome. It is calculated as TUD(kmer) = Observed(kmer)/Expected(kmer).
    
    Input:
        kmer_count - a dictionary mapping the kmers found in a sequence to its counts (Observed).
        expected_zero - a dictionary mapping the kmers to its calculated expected values (Expected).
    
    Ouputs:
        tud_n - a dictionary mapping the kmers to its tetranucleotide usage deviation (TUD).
    """
    # create a list of kmer from the count dictionary
    kmer_list = list(kmer_count.keys())
    # create the keys and the countainer 
    tud_n = defaultdict(float,[(km, 0.0) for km in kmer_list])
    # for each kmer in the list
    for km in kmer_list:
        # calculate the ratio of observed/expected values
        norm = kmer_count[km]/expected_zero[km]
        # add the values to it respective key
        tud_n[km] = norm
    # return the result
    return tud_n    
    
    
def get_kmer_cnts_from_csv(filename, k):
    """
    Read a csv  file with kmer and counts from a given
    sequence.
    Inputs:
        filename - string representing a path to the csv file.

    Outputs:
        kmer_counts -  a dictionary-like object mapping the kmer
                       of length k to the number of the kmer found/
                       counted in a particular genome.

    """
    # create the dictionary object
    kmer_counts = defaultdict(int)
    # opn the file
    with open(filename, 'r') as file:
        # create a csv reader object
        csv_reader_object = csv.reader(file)
        # iterates through the rows of the csv obj
        for row in csv_reader_object:
            kmer, cnt = row[0], float(row[1])
            # add the kmer and it counts to the counter
            if len(kmer) == k:
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + int(cnt)
    return kmer_counts    


def get_oligo_usage_deviation_normalized(kmer_count, expected_zero):
    """
    The normalized value for a word W is calculated by dividing the
    observed counts by the expected counts. This is the usage deviation
    vector for a genome. It is calculated as TUD(kmer) = Observed(kmer)/Expected(kmer).
    
    Input:
        kmer_count - a dictionary mapping the kmers found in a sequence to its counts (Observed).
        expected_zero - a dictionary mapping the kmers to its calculated expected values (Expected).
    
    Ouputs:
        tud_n - a dictionary mapping the kmers to its tetranucleotide usage deviation (TUD).
    """
    # create a list of kmer from the count dictionary
    kmer_list = list(kmer_count.keys())
    # create the keys and the countainer 
    oud_n = defaultdict(float,[(km, 0.0) for km in kmer_list])
    # for each kmer in the list
    for km in kmer_list:
        # calculate the ratio of observed/expected values
        norm = kmer_count[km]/expected_zero[km]
        # add the values to it respective key
        oud_n[km] = signif(norm)
    # return the result
    return oud_n
    
    
def get_oligo_usage_deviation_zscores(diff_by_window_chunk):
    # diff_by_window_chunk is a list of values calculated by
    # oligo_nucleotide_usage_deviation_slide_window
    len_lst = len(diff_by_window_chunk)
    # get the mean
    mean = np.mean(diff_by_window_chunk)
    # get the standard deviation
    sigma = np.std(diff_by_window_chunk)
    # get z scores
    z_scores = [(diff_by_window_chunk[i] - mean) / sigma for i in range(len_lst)]
    return z_scores    
    
    
def oligo_nucleotide_usage_deviation(sequence, window, step, k):
    # oligo nucleotide usage
    seq_len = len(sequence)
    # get all possible kmers of length k
    kmer_list = get_all_possible_kmers('ACGT', k, k)
    # get the bases/nucleotides frequencies
    base_freqs = base_stats(sequence, 'ACGT', as_count=False, as_dict=True)
    # count the kmer pf length k in the sequence
    obs = get_kmer_counts(sequence, k, as_freq=False)
    # calculate the expected values for each kmer of k length
    exp = kmer_expected_zero_markov(base_freqs,  kmer_list, seq_len)
    # calculates the oligo nucleotide deviation normalized
    oud_n = get_oligo_usage_deviation_normalized(obs, exp)
    # slide along the genome in windows defined by window length 
    # by steps defined by step length
    start = 0
    end = window
    w_idx = []
    diff_by_window = []
    # while theres is sequence
    while end < seq_len:
        if end > seq_len:
            end = seq_len
        # chunk the sequence 
        chunk = sequence[start:end]
        #record list of windows
        w_idx.append((start,end))
        # calculate observed and expected in the window
        # returns a dictionary 
        chunk_len = len(chunk)
        chunk_bases = base_stats(chunk, 'ACGT', as_count=False, as_dict=True)
        obs_chunk = get_kmer_counts(chunk, k, as_freq=False)
        exp_chunk = kmer_expected_zero_markov(chunk_bases,  kmer_list, chunk_len)

        #compute difference sum for all kmers in the chunks
        # using the usage deviation from the complete sequence 
        diff_sum = 0
        for kmer in obs.keys():
            diff_sum += abs((obs_chunk[kmer] / exp_chunk[kmer]) - oud_n[kmer])
        diff_by_window.append(diff_sum)
        #slide along the window by step length
        start += step
        end += step
    # compute Z-score as Z=(x-mu)/sigma
    zscores = get_oligo_usage_deviaiton_zscore(diff_by_window)
    # get all the start positions of the chunks
    starts = [x for x,y in w_idx]
    return starts, zscores 
    
    
def parse_arguments():
    """Parse the command line arguments to the oligo_nucleotide_usage_deviation script.
    Sets up the argparse command-line parser and calls it. These args can be accessed
    using args.args.
    """
    parser = argparse.ArgumentParser(
        description="""A script to calculate the oligo nucleotide usage deviation
         from bacterial genomes/plasmids.""",
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-di',
                        '--dir_in',
                        metavar='path',
                        type=str,
                        required=True,
                        dest='dir_in',
                        help='Directory root. In my case the name is conjugated with a subdir')

    parser.add_argument('-do',
                        '--dir_out',
                        type=str,
                        dest='dir_out',
                        help='directory name for output files.')

    parser.add_argument('-sd',
                        '--sub_dir',
                        type=str,
                        dest='sub_dir',
                        help='Name for a subdirectory, ex., Chromosomes.')

    parser.add_argument('-sw',
                        '--slide_window',
                        action='store_true',
                        dest='window',
                        help='Integer representing the size of the sub sequence')

    parser.add_argument('-w',
                        '--window',
                        type=int,
                        action="store",
                        dest='window',
                        help='Integer representing the size of the sub sequence')

    parser.add_argument('-k',
                        '--length',
                        type=int,
                        action="store",
                        dest='k',
                        help='Integer representing the size of the sub sequence')

    parser.add_argument('-s',
                        '--step',
                        type=int,
                        action="store",
                        dest='step',
                        help='Integer representing the size of the overlap sequence')

    parser.add_argument('-t',
                        '--seq_type',
                        type=str,
                        dest='seq_type',
                        help='Type of sequence, ex, chromosome or plasmid')
    
    return parser.parse_args()


def main():
    # starting count the staring time of the script
    start = time()
    # checking the current directory and printing it
    cwd = os.getcwd()
    print(colored(f'\nThe working directory: {cwd}\n',
                  'green',
                  attrs=['bold']))
    # passing the arguments to the script
    args = parse_arguments()
    # name of the input directory, ex. Data/Genomes_splitted
    dir_in = args.dir_in
    # name of the sub directory to save the final result
    # Chromosomes/Plasmids
    sub_dir = args.sub_dir
    # name of the root directory to save the final result
    dir_out = args.dir_out
    # kmer length
    k = args.k
    # window length
    window = args.window
    # size of the step (length overlap)
    step = args.step
    # if comparing within genome
    slide_window = args.slide_window
    # type of sequence (chr/plsm
    seq_type = args.seq_type
    # get the fasta files
    # glob.glob(f'Data/Genomes_splitted/*/Chromosomes/*_chr.fna.gz')
    fasta_files = glob.glob(f'{dir_in}/*/{sub_dir}/*_chr.fna.gz')
    # csv files
    csv_files = glob.glob(f'{dir_out}/kmer_counts/*/*_k2_8_{seq_type}_all.csv'')
    # check if the output directory exist other wise create it
    if os.path.exists(dir_out):
        print(colored('The directory to save the files already exists!',
                      'red',
                      attrs=['bold']))
        pass
    else:
        make_me_a_folder(dir_out)

    # analysing the oligo deviation within genomes
    if slide_window:
        for filename in fasta_files:
            name = filename.split('/')[2]
            for Id, sequence in fasta_parser.parse_fasta(filename):
                print(f'Calculating the oligonucleotide deviation from genus {name} sequence {Id}')
                oud = oligo_nucleotide_usage_deviation(sequence, window, step, k)
                
                
                
    else:
        # analysing oud in a complete genome
        for csvfile in csv_files:
                name = filename.split('/')[2]
                # get the count from all genomes (it is a mean count)
                kmc = get_kmer_cnts_from_csv(filename, k)
                # get the length of the sequence
                seq_len = genome_length('Results/Length/All_Chromosomes_length.csv', name)
                # calculates the expected values from couunt with zero order markov
                exp = kmer_expected_zero_markov(base_freqs,  kmer_list, seq_len)
                # calculates the olio usage deviation
                print(f'Calculating the oligonucleotide deviation from genus {name} sequence {Id}')                
                # returns a tuple of initial positions and the ratio obs/exp
                oud = get_oligo_nucleotide_usage_deviation_normalized(kmc, exp)
                
                full_path = os.path.join()
                if os.path.exists(full_path):
                    pass
                else:
                    os.makedirs(full_path)
            df_gc_slw_dif_slw.to_csv(f'{full_path}/{csv1}', index=False)
            df_gc_var_gc_total.to_csv(f'{full_path}/{csv2}', index=False)
    
    end = time()
    # print some info
    print(colored(f"Total number of genus/species analyzed: {len(filenames)}\n.",
                  attrs=['bold']))
    print(colored(f'Total time for the script finishes: {round(end - start, 2)}.',
                  'red',
                  attrs=['bold']))
    print(colored('Done!',
                  'green',
                  attrs=['bold']))


if __name__ == "__main__":
    sys.exit(main())
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
