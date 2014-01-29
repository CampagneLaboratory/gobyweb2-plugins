# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'INDEX' )
            ORGANISM=$3
            BUILD_NUMBER=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"

            LASTDB=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/bin/lastdb

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_LAST_INDEX_TOPLEVEL_IDS_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

            NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta
            ${LASTDB} -s 25G index ${INPUT_FASTA_NO_GZ}
            if [ $? != 0 ]; then
               return 1;
            fi

            cp -r index* ${installation_path}/
            ls -l
            ls -l ${installation_path}/

            if [ -e ${installation_path}/index.prj ]; then
               return 0
            else
               return 127
            fi

            ;;

        'TOPLEVEL_IDS' )
                    ORGANISM=$3
                    BUILD_NUMBER=$4
                    ENSEMBL_RELEASE=$5
                    echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"

                    . ${RESOURCES_GOBY_SHELL_SCRIPT}

                    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
                    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
                    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

                    INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

                    goby fasta-to-compact ${INDEXED_GENOME_DIR}/genome-toplevel.fasta  --exclude-sequences  --include-identifiers -o ${installation_path}/toplevel-ids.compact-reads

                    if [ -e ${installation_path}/toplevel-ids.compact-reads ]; then
                          return 0
                    else
                          return 127
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
    set +xv

    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
    ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
    echo >>${out} "organism=${ORGANISM}"
    echo >>${out} "reference-build=${BUILD_NUMBER}"
    echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"
    set -xv
    echo "Printing result from ${out}:"
    cat ${out}
    exit 0
}
