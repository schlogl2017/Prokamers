#!/usr/bin/env python

import sys
import glob
import shutil


def get_id_number(data, ext):
    d = []
    for fasta in data:
        name = fasta.split("/")[-1].replace(f"{ext}", "")
        d.append(name)
    return d


# Bacteria
supergroup = sys.argv[1]
# Deferribacteres
group = sys.argv[2]
# Deferribacteres
subgroup = sys.argv[3]


dir_cds = f"Genomes/{supergroup}/{group}/{subgroup}/CDS"
dir_chr = f"Genomes/{supergroup}/{group}/{subgroup}/Chr"

ext_cds = "_cds.fna.gz"
ext_chr = "_chr.fna.gz"

cds = glob.glob(f"{dir_cds}/*{ext_cds}")
chrs = glob.glob(f"{dir_chr}/*{ext_chr}")

print(f"Number CDS: {len(cds)}")
print(f"Number Genomes: {len(chrs)}")


set_chr = get_id_number(chrs, ext_chr)
set_cds = get_id_number(cds, ext_cds)

print(f"Number CDS Ids: {len(cds)}")
print(f"Number Genomes Ids: {len(chrs)}")

diff = set(set_chr).difference(set(set_cds))
print(f"Number of files to copy {len(diff)}")


count_files = 0
for Id in diff:
    count_files += 1
    src_path = f"{dir_chr}/{Id}{ext_chr}"
    # Temp/supergroup/group/subgroup
    dst_path = f"Temp/{supergroup}/{group}/{subgroup}/{Id}{ext_chr}"
    shutil.copy(src_path, dst_path)

print(f"Number of copied: {count_files}")
print("Done!!!")












