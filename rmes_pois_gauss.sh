#! usr/bin/env bash

# Genomes_new/Bacteria/Pseudomonadota/Acidithiobacillia/Acidithiobacillaceae_bacterium/GCA_016936895.1_chr.fna
dir_in="Genomes_new/$1/$2/$3"
dir_out="Rmes/$1/$2/$3"
mkdir -p "${dir_out}"

# organelles *.fasta
# bac/arc *_chr.fna
#viruses *.fna
for file in "${dir_in}"/*.fna
do
  #echo "$file"
  name=$(basename "$file" )
  name=${name%.chr}
  #echo "Calculating the compound poisson from ${name} file"
  #rmes --compoundpoisson -l "$3" --max --dna -s "${file}" -o "${dir_out}/${name}_k$3"_comp_pois
  echo "Calculating the gaussian data from ${file} file"
  rmes --gauss -l "$3" --max --dna -s "${file}" -o "${dir_out}/${name}_k$3_gauss"
done

