#! usr/bin/env bash
# Script to get the number of sequences in fasta files

# Genomes/Bacteria/Acidobacteria/Acidobacteria_incertae_sedis/[Chr/CDS]
INPUTDIR="Genomes/$1/$2/$3/CDS"
OUTPUTDIR="Results/Num_genes/$1/$2/$3"
type_seq=$4

# check the folders
if [[ ! -e $OUTPUTDIR ]]; then
    mkdir -p "${OUTPUTDIR}"
fi

for fasta in "${INPUTDIR}"/*_"$type_seq".fna.gz; # *_CDS.fasta.gz
do
  if test -f "$fasta"; then
    echo
  fi
  name=$(basename "${fasta}" _"$type_seq".fna.gz) #_CDS.fasta.gz)
  echo -ne "${name}    "
  zgrep -c "^>" "${fasta}"
done

