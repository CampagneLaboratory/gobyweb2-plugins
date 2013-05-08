# Installation script for Pathogen Data
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'FUNGI' )

            #${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/" "fungi.*.rna.fna.gz"
            #${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/" "fungi.*.rna.fna.gz"
            ${JOB_DIR}/fetch_url_pattern "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/" "fungi.*.rna.fna.gz"
            ${JOB_DIR}/fetch_url_pattern "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/fungi/" "fungi.*.genomic.fna.gz"

            echo working dir `pwd`
            #Concatenate them
            cat *.fna.gz > joined.fa.gz
            #Extract them to one big fasta file
            cat *.fna.gz > fungi.all.fna.gz
            gunzip fungi.all.fna.gz
            #Grep fasta file for lines starting with '>' and write the output to another file
            grep -a ">" fungi.all.fna > fungi.all.names.map
            #Cut the first character from each of those lines (the >)
            cat fungi.all.names.map | cut -c2- > fungi.all.names.edited.map
            awk '{print "\t"$0}' fungi.all.names.edited.map
            awk '{print $1"\t"$1=""; print $0}' fungi.all.names.edited.map  | more
            awk '{$1=$1"\t"; print $0}' fungi.all.names.edited.map  > fungal-names.map
            lastdb -Q 0 -v fungiref fungi.all.rna.fna

            tar -cvzf fungalref.tar.gz fungalref*

            mkdir -p "${installation_path}/fungal"
            cp fungalref.tar.gz "${installation_path}/fungal"
            cp fungal-names.map "${installation_path}/fungal"

            rm  -rf *.*

            return 0
            ;;

        'MICROBIAL' )
            #${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/microbial/" "microbial.*.rna.fna.gz"
            #${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/microbial/" "microbial.*.rna.fna.gz"
            ${JOB_DIR}/fetch_url_pattern "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/microbial/" "microbial.*.rna.fna.gz"
            ${JOB_DIR}/fetch_url_pattern "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/microbial/" "microbial.*.genomic.fna.gz"

            #Concatenate them
            cat *.fna.gz > joined.fa.gz
            #Extract them to one big fasta file
            cat *.fna.gz > micro.all.fna.gz
            gunzip micro.all.fna.gz
            #Grep fasta file for lines starting with '>' and write the output to another file
            grep -a ">" micro.all.fna > micro.all.names.map
            #Cut the first character from each of those lines (the >)
            cat micro.all.names.map | cut -c2- > micro.all.names.edited.map
            awk '{print "\t"$0}' micro.all.names.edited.map
            awk '{print $1"\t"$1=""; print $0}' micro.all.names.edited.map  | more
            awk '{$1=$1"\t"; print $0}' micro.all.names.edited.map  > micro-names.map
            lastdb -Q 0 -v microref micro.all.rna.fna

            tar -cvzf microref.tar.gz microref*

            mkdir -p "${installation_path}/bacterial"
            cp microref.tar.gz "${installation_path}/bacterial"
            cp micro-names.map "${installation_path}/bacterial"

            rm  -rf *.*

            return 0
            ;;
        'VIRAL' )
              #${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/" "viral.*.rna.fna.gz"
            #${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/" "viral.*.rna.fna.gz"
            ${JOB_DIR}/fetch_url_pattern "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/" "viral.*.rna.fna.gz"
            ${JOB_DIR}/fetch_url_pattern "ftp://ftp.ncbi.nlm.nih.gov/refseq/release/viral/" "viral.*.genomic.fna.gz"
            echo working dir `pwd`
            #Concatenate them
            cat *.fna.gz > joined.fa.gz
            #Extract them to one big fasta file
            cat *.fna.gz > viral.all.fna.gz
            gunzip viral.all.fna.gz
            #Grep fasta file for lines starting with '>' and write the output to another file
            grep -a ">" viral.all.fna > viral.all.names.map
            #Cut the first character from each of those lines (the >)
            cat viral.all.names.map | cut -c2- > viral.all.names.edited.map
            awk '{print "\t"$0}' viral.all.names.edited.map
            awk '{print $1"\t"$1=""; print $0}' viral.all.names.edited.map  | more
            awk '{$1=$1"\t"; print $0}' viral.all.names.edited.map  > viral-names.map
            lastdb -Q 0 -v viralref viral.all.rna.fna

            tar -cvzf viralref.tar.gz viral*

            mkdir -p "${installation_path}/viral"
            cp viralref.tar.gz "${installation_path}/viral"
            cp viral-names.map "${installation_path}/viral"

            rm  -rf *.*

            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}