cat Metadata/Archaea/Metadata_Asgard/Cand_Baldrarchaeota.tsv | while read -r acc ; do
    esearch -db assembly -query $acc </dev/null \
        | esummary \
        | xtract -pattern DocumentSummary -element FtpPath_GenBank \
        | while read -r url ; do
            fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
            wget -nc "$url/$fname" ;
        done ;
    done
