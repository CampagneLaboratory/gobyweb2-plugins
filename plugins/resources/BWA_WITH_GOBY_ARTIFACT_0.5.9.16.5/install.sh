# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    set -x
    case ${id} in

        'EXECUTABLE' )
            . ${RESOURCES_ARTIFACTS_GOBY_CPP_API_LIBRARIES}/setup.sh

            (${RESOURCES_FETCH_URL_SCRIPT} http://campagnelab.org/files/bwa-0.5.9-goby-2.3.zip  bwa-goby-support.zip

                unzip bwa-goby-support.zip
                cd bwa-goby-support/
                chmod +x autogen.sh
                ./autogen.sh
                ./configure --with-goby --prefix=${installation_path}
                make
                make install
                # cp bwa ${installation_path}/bwa-icb
                #chmod +x ${installation_path}/bwa-icb
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

            BWA_EXEC=${RESOURCES_ARTIFACTS_BWA_WITH_GOBY_ARTIFACT_EXECUTABLE}/bin/bwa

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

            #INPUT_FASTA_NO_GZ=genome.fasta

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta
            if [ "${COLOR_SPACE}" = "true" ]; then
                 SPACE_PARAM="-c"
            fi

            ${BWA_EXEC} index -a bwtsw ${SPACE_PARAM} -p index ${INPUT_FASTA_NO_GZ}
            STATUS=$?
            if [ ${STATUS} != 0 ]; then

             return ${STATUS}
            fi
            ls -ltr
            cp index* ${installation_path}/

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