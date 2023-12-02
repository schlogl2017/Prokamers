#! usr/bin/env python

import gzip
from Bio import SeqIO
gbk_filename = ""
faa_filename = ""
input_handle  = gzip.open(gbk_filename, "rt")
output_handle = open(faa_filename, "w")

​

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print(f"Dealing with GenBank record {seq_record.id}")
    output_handle.write(f">{seq_record.id}\n{seq_record.seq}\n")

​

input_handle.close()
output_handle.close()
