# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'SAMTOOLS_FAI_INDEX' )
            ORG=$3
            BUILD_NUMBER=$4
            ENSEMBL_RELEASE=$5

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

            # Link to the genome in the indexed directory, then put the index next to the link:
            gzip -c -d ${GENOME_DIR}/genome-toplevel.fasta.gz >${installation_path}/genome-toplevel.fasta
            ${RESOURCES_ARTIFACTS_SAMTOOLS_BINARIES}/samtools faidx ${installation_path}/genome-toplevel.fasta

            if [ -e ${installation_path}/genome-toplevel.fasta.fai ]; then
                return 0
            else
                return 1
            fi
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}


function get_attribute_values() {

    id=$1
    out=$2

    echo "get_attribute_values for ID=${id}"

    set -xv
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
    ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
    echo >>${out} "organism=${ORGANISM}"
    echo >>${out} "reference-build=${BUILD_NUMBER}"
    echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"

    echo "Printing result from ${out}:"
    cat ${out}
    exit 0
}