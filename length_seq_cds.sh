#! usr/bin/env bash
# Script to get the number of sequences in fasta files
# bash length_seq_cds.sh Archaea Unclassified_viruses
# Genomes/Archaea/Cand_Thermoplasmota/Cand_DHVE2_Group/CDS/GCA_000025665.1_cds.fna

OUTPUTDIR="Results/Length_genes/$1/$2/$3"
echo "${OUTPUTDIR}"
# check the folders
#if [[ ! -e $OUTPUTDIR ]]; then
    #mkdir -p "${OUTPUTDIR}"
#fi

for fasta in Genomes/"$1"/"$2"/"$3"/CDS/*.fna; # *_CDS.fasta.gz
do
  echo "${fasta}"
  #if test -f "$fasta"; then
    #echo "${fasta}"
  #fi
  name=$(basename "${fasta}" _cds.fna) #_CDS.fasta.gz)
  echo "${name}"
  bioawk -c fastx '{ print $name, length($seq) }' "${fasta}" > "${OUTPUTDIR}"/"${name}"_len_cds.tsv
done
