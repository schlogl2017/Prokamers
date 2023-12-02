#! usr/bin/env bash
# Calculates the GC content in genomes in slidewindows

set -e

# Genomes/$1=Archaea/$2=Asgard/$3=Cand_Baldrarchaeota/Chr
results="Results/GC_sldw/$1/$2/$3"
infolder="Temp/$2/$3"
echo "$infolder"
# size of the windowz
window=$4

if [[ ! -e ${results} ]]; then
    mkdir -p "${results}"
elif [[ ! -d ${results} ]]; then
    echo "${results} already exists but is not a directory" 1>&2
fi

for file in Temp/"$2"/"$3"/*_chr.fna
do
  echo $(ls Temp/"$2"/"$3"/*_chr.fna)
  echo "${file}"
  #name=$(basename "$file")
  #echo "Calculating the GC content in a window of size $window from genome $name"
  # generate a fasta index (including sizes in second column)
  #samtools faidx "${filename}"
  # extract sequence names (chromosomes) and their lengths to a new file
  #cut -f 1,2 "${filename}".fai > "$infolder"/"${name}".sizes
  # create a BED file with windows of wished width GCA_001940655.1_chr.fna5000.bps.bed
  #bedtools makewindows -g "$infolder"/"${name}".sizes -w "${window}" > "$infolder"/"${name}"."${window}".bps.bed
  #bedtools nuc -fi "${filename}" -bed "${name}"."${window}".bps.bed > "${results}/${name}_nuc_${window}".txt
done


