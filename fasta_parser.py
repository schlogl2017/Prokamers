import gzip
import re
from Bio import SeqIO
from collections import defaultdict
from itertools import groupby


def is_header(line):
    return line[0] == '>'


def parse_fasta(filename):
    if filename.endswith('.gz'):
        opener = lambda filename: gzip.open(filename, 'rt')
    else:
        opener = lambda filename: open(filename, 'r')

    with opener(filename) as f:
        fasta_iter = (it[1] for it in groupby(f, is_header))
        for name in fasta_iter:
            name = name.__next__()[1:].strip()
            sequences = ''.join(seq.strip() for seq in fasta_iter.__next__())
            yield name, sequences.upper()


def get_name_sequence_biopython(filename):
    for rec in SeqIO.parse(gzip.open(filename, "rt"), "fasta"):
        name = rec.id
        seq = str(rec.seq)
        yield name, seq


def get_sequence_from_scaffolds(filename):
    return ''.join([seq for name, seq in parse_fasta(filename)])


def fasta_item_counter(filename):
    """It opens and check the number of items in the fasta file."""
    if filename.endswith('.gz'):
        count = sum(g for g, _ in groupby(gzip.open(filename, 'rt'), key=is_header))
    else:
        count = sum(g for g, _ in groupby(open(filename, 'r'), key=is_header))
    return count


def punctuation_strip(string):
    return re.sub(r'[^\w\d\s]+', '', string)


def join_dict_string_data(dic):
    return {k: ''.join(v) for k, v in dic.items()}


def get_name(header, genome_id):
    header_list = header.split(" ")
    name = "_".join(header_list[1:3] + [header_list[0]])
    return f"{genome_id}_{name}"


def get_name_v2(header, genome_id):
    header_list = header.split(" ")
    name = "_".join(header_list[1:3])
    return f"{genome_id}_{name}"


def get_fasta_headers(filename):
    headers = []
    for name, _ in parse_fasta(filename):
        headers.append(name)
    if len(headers) >= 1:
        return headers


def get_plasmids_chromosome_sequences_from_complet_fasta(filename):
    plasmids = defaultdict(list)
    chromosome = defaultdict(list)
    for header, seq in parse_fasta(filename):
        if 'plasmid' in header:
            plasmids[header] = plasmids.get(header, [])
            plasmids[header].append(seq)
        else:
            chromosome[header] = plasmids.get(header, [])
            chromosome[header].append(seq)
    for name, seq in plasmids.items():
        plasmids[name] = ''.join(seq)
    for name, seq in chromosome.items():
        chromosome[name] = ''.join(seq)
    return plasmids, chromosome


def split_sequences_from_fasta_file(filename):
    """
    Function that get a multifasta file and split it in to single fasta files sequence.
    Input:
    
        filename - a full path tho the multi fasta file. Ex; a assembly or a reference
                   ncbi file of a complete bacterial genome, with chromosome and
                   plasmid sequences.
    
    Output:
    
        plasmids, chromosome - a tuple of dictionary-like objects with chromosome and
                               plasmids splitted from a input multi fasta file.                   
    """
    plasmids = defaultdict(list)
    chromosome = defaultdict(list)
    # get sequence ID from the filename. Ex. 
    # ['Data/Assemblies/Acidiphilium/GCF_000202835.1.fna.gz']
    # GCF_000202835.1
    genome_id = filename.split('/')[-1][:15]
    # read the file and get the headers and sequences
    for header, seq in parse_fasta(filename):
        header_clean = punctuation_strip(header)
        # checking for plasmid sequences in the file
        # and get them in the map
        if 'plasmid' in header:
            name = get_name_v2(header_clean, genome_id) + '_plas'
            plasmids[name] = plasmids.get(name, [])
            plasmids[name].append(seq)
        else:
            name = get_name_v2(header_clean, genome_id) + '_chr'
            chromosome[name] = chromosome.get(name, [])
            chromosome[name].append(seq)
    # joint the sequence list
    chro = join_dict_string_data(chromosome)
    # join the chromossome sequences
    plas = join_dict_string_data(plasmids)
    return plas, chro


def write_fasta(seq, name, out_file, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    sequence: sequence of strings -> seq
    out_file : str
        Path to write the sequences to.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    with open(out_file, 'w') as fout:
        fout.write(f'>{name}\n')
        for i in range(0, len(seq), wrap):
            fout.write(f'{seq[i:i + wrap]}\n')


def write_fasta_file(sequence_dict, out_file, wrap=80):
    """Write sequences to a fasta file.

    Parameters
    ----------
    sequence_dict : dict[seq_id] -> seq
        Sequences indexed by sequence id.
    out_file : str
        Path to write the sequences to.
    wrap: int
        Number of AA/NT before the line is wrapped.
    """
    with open(out_file, 'w') as fout:
        for name, seq in sequence_dict.items():
            fout.write(f'>{name}\n')
            for i in range(0, len(seq), wrap):
                fout.write(f'{seq[i:i + wrap]}\n')
