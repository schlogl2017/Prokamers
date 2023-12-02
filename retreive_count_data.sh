#!usr/bin/env bash
# selecting data from kmer counts

#Inputs
# Genomes_new/Archaea/Asgard_group/Asgard_group_names.txt
# Counts/Archaea/Asgard_group C_Heimdallarchaeota CHR 2

g=$1 # Archaea
sg=$2 # Asgard_group
ssg=$3 # C_Heimdallarchaeota
sqtype=$4 # CHR/PLAS
k=$5 # length kmer 
fileaid=Genomes_new/"$g"/"$sg"/"$ssg"/"$sqtype"/"$ssg".ids
dirin=Counts/"$g"/"$sg"/"$ssg"/"$sqtype"


for id in $(cat $fileaid); 
do
    echo "$id"
    #echo "$dirin"
    #echo "$id"_"$k".counts
    awk -F, -v awkVar="$k" '{if(length($1) == awkVar ) print }' "$dirin"/"$id"_*_kmer.counts > "$dirin"/"$id"_"$k"_kmer.counts
done
