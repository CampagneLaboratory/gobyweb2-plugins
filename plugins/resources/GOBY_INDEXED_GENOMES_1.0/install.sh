# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    . ${SGE_O_WORKDIR}/constants.sh
    . ${SGE_O_WORKDIR}/auto-options.sh
    . ${RESOURCES_GOBY_SHELL_SCRIPT}

    case ${id} in


        'SEQUENCE_CACHE' )
            ORG=$3
            BUILD_NUMBER=$4
            ENSEMBL_RELEASE=$5

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

            run-goby build-sequence-cache -i ${GENOME_DIR}/genome-toplevel.fasta.gz -o random-access-genome
            cp random-access-genome* ${installation_path}/


            if [ -e ${installation_path}/random-access-genome.sizes ]; then
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


function get_attribute_values() {

    id=$1
    out=$2

    echo "get_attribute_values for ID=${id}"

    # get environment variables for GobyWeb job from SGE work directory:
    . ${SGE_O_WORKDIR}/constants.sh
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