    # Installation script for ENSEMBL_GENOMES
function plugin_install_artifact {

    id=$1
    installation_path=$2
set -x
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
                install_custom_fasta ${ENSEMBL_RELEASE}
            else
                if ["${BUILD_NUMBER}" = "1000GENOMES" ]; then
                 fail "1000GENOMES is not supported for transcripts."
                else
                ftp://ftp.ensembl.org/pub/release-82/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
                    ${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ORG_LOWERCASE}/cdna/" ".dna.toplevel.fa.gz"

                fi
            fi
            cp *.cdna.all.fa.gz ${installation_path}/transcripts-all.fasta.gz
            if [ -e ${installation_path}/transcripts-all.fasta.gz ]; then

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

   FASTA_FILE_TAG=$1
   # TODO: use fileset --fetch with tag to retrieve the file (to be implemented in the plugins-sdk command line).
   # TODO: For now, just copy from a fixed directory.
   cp ~gobyweb/customs/custom-ref-${FASTA_FILE_TAG}.fa.gz custom-${FASTA_FILE_TAG}.dna.toplevel.fa.gz
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