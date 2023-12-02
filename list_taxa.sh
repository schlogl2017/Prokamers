#!usr/bin/env bash

IN="Genomes_new/$1/$2/"
for name in $(ls -1 "$IN");
do
  echo  "$name" >> "$IN"/"$2"_names.txt
done
