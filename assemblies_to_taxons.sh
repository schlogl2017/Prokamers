#! usr/bin/env bash
INPUT=$1
for line in $(cat $INPUT); 
do
  echo -n "${line} "
  esearch -db assembly -query "${line}" | esummary | xtract -pattern DocumentSummary -element Taxid,SpeciesName,ScientificName;
done

