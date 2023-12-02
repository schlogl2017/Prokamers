from pysam import FastxFile
 
 
def fasta_to_hash(fasta_file):
    """Parse FASTA file into a dict.
 
    Args:
        fasta_file (str): Path to FASTA file.
 
    Returns:
        dict: key is the contig id and value is the sequence.
 
    """
    contigs = {}
    with FastxFile(fasta_file) as fh:
        for entry in fh:
            contigs[entry.name] = entry.sequence
    return contigs
 
 
def extract_phage_sequences(prophage_tbl, input_fasta, output_fasta):
    """Extract phage sequences from the bacterial genome.
 
    Args:
        prophage_tbl (str): Path to Phispy's prophage.tbl output.
        input_fasta (str): Path to input file in FASTA.
        output_fasta (str): Path to output file in FASTA.
 
    """
    contigs = fasta_to_hash(input_fasta)
    with open(prophage_tbl) as my_file, open(output_fasta, "w+") as output:
        for row in my_file:
            sequence_id, coords = row.strip().split()
 
            contig_id, prophage_start, prophage_end = coords.split("_")
 
            prophage_genome = contigs[contig_id][int(prophage_start) - 1:int(prophage_end)]
 
            output.write(">{}_{}_start_{}_end_{}\n{}\n".format(sequence_id, contig_id, prophage_start, prophage_end,
                                                               prophage_genome))
 
    print("Check for: {}\nDone :)".format(output_fasta))
 
 
extract_phage_sequences("U00096.2_annoation_phages/prophage.tbl", "U00096.2.fasta", "U00096.2_prophages_genomes.fasta")
