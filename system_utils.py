#!usr/bin/env python
# -*- coding: utf-8 -*-


import os
import glob
from collections import defaultdict
import pandas as pd


def get_files(dir_name):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_lengths(fasta_dict):
    seq_len = defaultdict(int)
    file_names = list(fasta_dict.values())[0]
    names = [os.path.splitext(os.path.split(name)[-1])[0] for name in file_names]
    for name, filename in zip(names, file_names):
        for n, seq in parse_fasta(filename):
            seq_len[name] = seq_len.get(name, 0) + len(seq)
    return seq_len


def get_files_paths(dir_name, sub_dir_name):
    """Returns a dictionary like object using baterial genera
    as keys and a list of path to fasta fileas as values.
    Inputs:
        dir_name: directory name
        sub_dir_name : sub directory name
    Outputs:
        dictionary: with genera name and a list of 
                    files.
    
    Example:
    dir is like: Data/bacteria_splitted/Mycolicibacillus/chromosomes
                 Data/bacteria_splitted/Mycolicibacillus/plasmids (if the case)
    dirname: 'Data/bacteria_splitted'
    sub_dir: 'chromosomes'
    fasta_dicts = get_files(dir_name, 'chromosomes')
    fasta_dicts['Mycolicibacillus']
    ['NZ_AP022594.1_Mycolicibacillus_koreensis_strain_JCM_19956.fna.gz']
    """
    # create a list of file and sub directories
    # names in the given directory
    spc_names = sorted(os.listdir(dir_name))
    all_files = defaultdict(list)
    for name in spc_names:
        all_files[name] = all_files.get(name, [])
        paths = os.path.join(dir_name, name, sub_dir_name)
        for file in os.listdir(paths):
            full_paths = os.path.join(paths, file)
            all_files[name].append(full_paths)
    return all_files


def get_names(path):
    return [(file.split('/')[3][:-3]) for file in get_files(path)]


def get_dir_names(filenames):
    dirn = set()
    subd = set()
    subsub = set()
    for filename in filenames:
        dir_splt = filename.split('/')
        dirn.add(dir_splt[0])
        subd.add(dir_splt[1])
        subsub.add(dir_splt[2])
    return next(iter(dirn)), next(iter(subd)), next(iter(subsub))


def get_dir_name(filename):
    subd, subsub = set(), set()
    names = filename.split('/')
    for i in range(len(names)):
        subd.add(names[1])
        subsub.add(names[2])
    return next(iter(subd)), next(iter(subsub))


def make_me_a_folder(folder_name):
    os.getcwd()
    try:
        os.makedirs(folder_name)
    except OSError:
        print("Directory already exists")
    return os.path.exists(folder_name)


def get_full_name(dir_out, sud_dir, ssub):
    return os.path.join(dir_out, sud_dir, ssub)


def get_fasta_files(dir_name):
    infiles = []
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            input_files = os.path.join(path, name)
            if input_files.endswith('fa.gz') or input_files.endswith('.fa') \
                    or input_files.endswith('.fasta') or input_files.endswith('.fa.gz') \
                    or input_files.endswith('.fna') or input_files.endswith('.fna.gz'):
                infiles.append(input_files)
    return infiles


def get_files2(dir_name, subdir):
    # create a list of file and sub directories
    # names in the given directory
    files = os.listdir(dir_name)
    all_files = []
    # Iterate over all the entries
    for entry in files:
        # Create full path
        full_path = os.path.join(dir_name, entry, subdir)
        # If entry is a directory then get the list of files in this directory
        if os.path.isdir(full_path):
            all_files += get_files(full_path)
        else:
            all_files.append(full_path)
    return all_files


def get_list_paths(spc_names, path, ssub):
    pwd = defaultdict()
    for name in spc_names:
        pwd[name] = get_full_name(path, name, ssub)
    return pwd


def get_fasta_files_paths(dict_paths):
    dict_fasta_paths = defaultdict(list)
    for name, path in dict_paths.items():
        dict_fasta_paths[name] = dict_fasta_paths.get(name, []) + glob.glob(path)
    return dict_fasta_paths


def merge_two_dicts(x, y):
    """Given two dicts, merge them into a new dict as a shallow copy."""
    z = x.copy()
    z.update(y)
    return z


def get_fasta_files_dict(dir_name):
    infiles = defaultdict(list)
    for path, subdirs, files in os.walk(dir_name):
        for name in files:
            infiles[name] = infiles.get(name, [])
            input_files = os.path.join(path, name)
            if input_files.endswith('fa.gz') or input_files.endswith('.fa') \
                    or input_files.endswith('.fasta') or input_files.endswith('.fa.gz') \
                    or input_files.endswith('.fna') or input_files.endswith('.fna.gz'):
                infiles[name].append(input_files)
    return infiles


def fasta_file_paths(fasta_dirs):
    fasta_paths = defaultdict(list)
    exts = tuple(['fa.gz', '.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz'])
    for name, pwd in fasta_dirs.items():
        fasta_paths[name] = fasta_paths.get(name, [])
        for path, _, files in os.walk(pwd):
            for n in files:
                input_files = os.path.join(path, n)
                if input_files.endswith(exts):
                    fasta_paths[name].append(input_files)
    return fasta_paths


def get_species_name(filenames):
    species = []
    with open(filenames, 'r') as fh:
        for line in fh:
            sp = line.strip()
            species.append(sp)
    return species


def find_csv_filenames(dir_name, spc_name, sub_dir, ssub, num_k):
    name = os.path.join(dir_name, spc_name, sub_dir)
    path_files = [os.path.join(name, f'{ssub}{str(i)}') for i in range(2, num_k)]
    files_dict = defaultdict(list)
    files_dict[spc_name] = files_dict.get(spc_name, [])
    for path in path_files:
        filenames = ''.join(os.listdir(path))
        files_dict[spc_name].append(os.path.join(path, filenames))
    return files_dict


def get_csvs_to_df_concatenation(spc_name, filenames_dict):
    """Function tha receives a species name and a dictionary with
    all csvs files to be concatenated.
    """
    dfs = []
    for filename in filenames_dict[spc_name]:
        dfs.append(pd.read_csv(filename))
    dfs = [df.set_index("kmers", drop=True) for df in dfs]
    concat = pd.concat(dfs, axis=0, copy=False).reset_index()
    return concat, len(dfs)


def csv_filenames(dir_name, spc_name, sub_dir, ssub):
    name = os.path.join(dir_name, spc_name, sub_dir)
    path_files = [os.path.join(name, ssub)]
    files_dict = defaultdict(list)
    files_dict[spc_name] = files_dict.get(spc_name, [])
    for path in path_files:
        filenames = ''.join(os.listdir(path))
        files_dict[spc_name].append(os.path.join(path, filenames))
    return files_dict


def get_fasta(spc_list, dir_in, sub_sub_dir):
    ends = ('fa.gz', '.fa', '.fasta', '.fna', '.fna.gz')
    infiles = []
    for dir_name in full_path_for_fasta(spc_list, dir_in, sub_sub_dir):
        for path, subdirs, files in os.walk(dir_name):
            for name in files:
                input_files = os.path.join(path, name)
                if input_files.endswith(ends):
                    infiles.append(input_files)
    for file in infiles:
        yield file


def full_path_for_fasta(spc_list, dir_in, sub_sub_dir):
    for spc in spc_list:
        path = os.path.join(dir_in, spc, sub_sub_dir)
        yield path
