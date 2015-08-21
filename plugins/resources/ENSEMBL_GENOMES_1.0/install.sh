    # Installation script for ENSEMBL_GENOMES
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


      'TOPLEVEL' )
            ORGANISM=$3
            GENOME_REFERENCE_ID=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID} ENSEMBL_RELEASE=${ENSEMBL_RELEASE}"
            ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`
            # get genome from local NFS, or wget from ensembl servers if not present
            BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`

            if [ "${BUILD_NUMBER}" = "CUSTOM" ]; then
                install_custom_fasta ${GENOME_REFERENCE_ID}
            else
                if [ "${BUILD_NUMBER}" = "GRCh37" -o "${BUILD_NUMBER}" = "1000GENOMES" ]; then
                    # For human, use the compatible 1000g assembly instead of the Ensembl build: (coordinates are compatible),
                    # see http://www.1000genomes.org/category/assembly
                    ${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz human.dna.toplevel.fa.gz
                else
                    ${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ORG_LOWERCASE}/dna/" ".dna.toplevel.fa.gz"

                fi
            fi
            cp *.dna.toplevel.fa.gz ${installation_path}/genome-toplevel.fasta.gz
            if [ -e ${installation_path}/genome-toplevel.fasta.gz ]; then

             exit 0
            else

             exit 1
            fi
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}
function install_custom_fasta() {
   id=$1
   FASTA_FILE_TAG=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
   # TODO: use fileset --fetch with tag to retrieve the file (to be implemented in the plugins-sdk command line).
   # TODO: For now, just copy from a fixed directory.
   cp ~gobyweb/tmp/custom-ref-${FASTA_FILE_TAG}.fa.gz custom-${FASTA_FILE_TAG}.dna.toplevel.fa.gz
}

function get_attribute_values() {

    id=$1
    out=$2

    echo get_attribute_values for ID=${id}

    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
    ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
    echo >>${out} "organism=${ORGANISM}"
    echo >>${out} "reference-build=${BUILD_NUMBER}"
    echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"

    echo "Printing result from ${out}:"
    cat ${out}
    return 0
}