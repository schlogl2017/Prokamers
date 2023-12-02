#! usr/bin/env bash
# $1 Archaea/Bacteria
# $2 Asgard/Terrabacteria
# $3 Cand_Baldrarchaeota 
# $4 Ch/Plas
# Genomes/Archaea/Asgard/Cand_Baldrarchaeota/Chr/GCA_019056805.1_chr.fna.gz
# Genomes/Archaea/Asgard/Cand_Heimdallarchaeota/Chr/GCA_003144275.1_chr.fna.gz
# bash seq_lengths_bioawk.sh Bacteria Nitrospinae_Tectomicrobia Chr > Results/Lengths/Bacteria/Nitrospinae_Tectomicrobia_chr_lengths.tsv

outfolder="Results/Lengths/$1/$2/$3"

infolder="Genomes/$1/$2/$3/$4/CDS"

if [[ ! -e $outfolder ]]; then
    mkdir -p "$outfolder"
fi
# for bac/arch/viruses *_chr.fna.gz *_cds_from_genomic.fna.gz *_pls.fna.gz
# for Mito/Plasmids *.fasta.gz
for fasta in "${infolder}"/*_cds_from_genomic.fna.gz
do
    bioawk -c fastx '{ print $name, length($seq) }' "${fasta}"
done
