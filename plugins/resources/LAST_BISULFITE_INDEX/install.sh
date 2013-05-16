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
            BISULFITE_FORWARD_SEED=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/examples/bisulfite_f.seed
            BISULFITE_REVERSE_SEED=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/examples/bisulfite_r.seed

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})


            NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
            SJDB=gencode.v14.annotation.gtf.sjdb

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta
            ${LASTDB} -u ${BISULFITE_FORWARD_SEED} index_f ${INPUT_FASTA_NO_GZ}
            if [ $? != 0 ]; then
               return 1;
            fi

            ${LASTDB} -u ${BISULFITE_REVERSE_SEED} index_r ${INPUT_FASTA_NO_GZ}
            if [ $? != 0 ]; then
               return 1;
            fi
            cp -r index_f* ${installation_path}/
            cp -r index_r* ${installation_path}/
            ls -l
            ls -l ${installation_path}/

            if [ -e ${installation_path}/index ]; then
               return 0
            else
               return 127
            fi
            return 127
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
