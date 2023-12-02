#! usr/bin/env bash
# Script to get the number of sequences in fasta files
# 
# for name in $(ls Genomes/Bacteria/Thermotogae);  do bash gc_percent_bioawk.sh Bacteria Thermotogae "${name}";  done

OUTPUTDIR="Results/GC/Total//$1/$2/$3"
echo "${OUTPUTDIR}"

if [[ ! -e ${OUTPUTDIR} ]]; then
    mkdir -p "${OUTPUTDIR}"
elif [[ ! -d ${OUTPUTDIR} ]]; then
    echo "${OUTPUTDIR} already exists but is not a directory" 1>&2
fi
# Chr
# Mito/Plastids *.fasta.gz
#Bacteria/Archaea *_chr.fna.gz
#Viruses *.fna.gz
for fasta in Genomes/"$1"/"$2"/"$3"/CDS/*_cds.fna; 
do
  #if test -f "$fasta"; then
    #echo "${fasta}"
  #fi
  name=$(basename "${fasta}" _cds.fna) #_CDS.fasta.gz)
  echo "${name}"
  echo "${fasta}"
  bioawk -c fastx '{ print $name, gc($seq) }' "${fasta}" > "${OUTPUTDIR}"/"${name}"_gc_cds.tsv
done
