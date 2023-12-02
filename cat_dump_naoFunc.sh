# for i in {1..10}; do cat "GCA_000020825.1_$i.dump" >> GCA_000020825.1_k1_10.txt; done

GROUP=$1
SUBGR=$2
SPC=$3
IN_DATA=Jelly_count/"${GROUP}"/"${SUBGR}"/"${SPC}"
IDS="jelly"
for id in $(cat "${IN_DATA}"/"${IDS}");
do
  #echo "${id}"
  find "${IN_DATA}"/CHR -type f -name "${id}"_*.dump -exec cat {} + >> "${id}"_kmer.counts
  awk '{ print length($1), $0 | "sort -n" }' "${id}"_kmer.counts > "${IN_DATA}"/CHR/"${id}"_kmer_sort.counts
  rm "${id}"_kmer.counts
done

