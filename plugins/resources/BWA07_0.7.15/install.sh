# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    set -x
    case ${id} in

        'EXECUTABLE' )
           VERISON="0.7.15"
            (${RESOURCES_FETCH_URL_SCRIPT} http://downloads.sourceforge.net/project/bio-bwa/bwa-0.7.15.tar.bz2?r=https%3A%2F%2Fsourceforge.net%2Fprojects%2Fbio-bwa%2Ffiles%2F&ts=1476387328&use_mirror=heanet  bwa-${VERSION}.bzip2
            bunzip2 bwa-${VERSION}.tar.bz2
            tar -xvf bwa-${VERSION}.tar
            cd    bwa-${VERSION}
            make
            cp bwa ${installation_path}/bwa
            chmod +x ${installation_path}/bwa
            )
            if [ -e ${installation_path}/bin/bwa ]; then
                return 0
            else
               return 127
            fi

            return 0
            ;;

        'INDEX' )
            ORGANISM=$3
            BUILD_NUMBER=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${BUILD_NUMBER} ${ENSEMBL_RELEASE} "

            BWA_EXEC=${RESOURCES_ARTIFACTS_BWA07_EXECUTABLE}/bin/bwa

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

            #INPUT_FASTA_NO_GZ=genome.fasta

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta
            if [ "${COLOR_SPACE}" = "true" ]; then
                 SPACE_PARAM="-c"
            fi

            ${BWA_EXEC} index -a bwtsw -p index ${INPUT_FASTA_NO_GZ}
            STATUS=$?
            if [ ${STATUS} != 0 ]; then
             return ${STATUS}
            fi
            ls -ltr
            cp index* ${installation_path}/
            if [ ! -f ${installation_path}/index.bwt ]; then
             exit 1;
            fi
            echo "Finished indexing Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"
            exit 0
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