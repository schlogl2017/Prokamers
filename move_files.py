#! usr/bin/env python
#  python move_files.py 5000 Archaea Euryarchaeota Archaeoglobi
#  python move_files.py 10000 Archaea Euryarchaeota Archaeoglobi
# python move_files.py 20000 Archaea Euryarchaeota Archaeoglobi


import os
import sys
import shutil
import glob



window = sys.argv[1] # 5000
supergroup = sys.argv[2] # Archaea 
group = sys.argv[3] # Asgard
subgroup = sys.argv[4] # Cand_Baldrarchaeota


ext = f"_chr.fna_nuc_{window}.txt"
dir_in = f"Metadata/{supergroup}/Metadata_{group}"
dir_out = f"Results/GC_sldw/{supergroup}/{group}/{subgroup}"
dir_files = f"Results/GC_sldw/{supergroup}/{group}/{subgroup}/*{ext}"

filenames =  glob.glob(f"Results/GC_sldw/{supergroup}/{group}/*{ext}")
#print(filenames)

with open(f"{dir_in}/{subgroup}.tsv", "r") as fh:
    header = fh.readline()
    for Id in fh:
        name = Id.strip().split("\t")[2]
        #print("Name", name)
        #name = Id.strip().split("\t")[4]
        print("Name", name)
        for genome in filenames:
            #print(genome)
            filename = genome.split("/")[-1].replace(ext, "")
            #print("filename", filename)
            if name == filename:
                shutil.move(genome, dir_out)
                print(f"Data tranfered to {dir_out}\n{genome}\n")







