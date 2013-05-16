# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    set -x
    case ${id} in

        'SCRIPTS' )

            VERSION="0.7.12"
            (${RESOURCES_FETCH_URL_SCRIPT} http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v${VERSION}.tar.gz  bismark_v${VERSION}.tar.gz
                gzip -c -d bismark_v${VERSION}.tar.gz |tar -xf -
                cd bismark_v${VERSION}
                cp * ${installation_path}/
                chmod +x ${installation_path}/bedGraph2cytosine
                chmod +x ${installation_path}/bismark2bedGraph
                chmod +x ${installation_path}/bismark_methylation_extractor
                chmod +x ${installation_path}/bismark_genome_preparation
                chmod +x ${installation_path}/bismark_genome_preparation
                chmod +x ${installation_path}/bismark
            )
            if [ -e ${installation_path}/bismark ]; then
                return 0
            else
               return 127
            fi


            ;;

        'INDEX' )
            ORGANISM=$3
            BUILD_NUMBER=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${BUILD_NUMBER} ${ENSEMBL_RELEASE} "

            BOWTIE2_EXEC=${RESOURCES_ARTIFACTS_BOWTIE2_ARTIFACT_BINARIES}/bowtie2

            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            cp ${FAI_INDEXED_GENOME_DIR}/* ${installation_path}/

            #INPUT_FASTA_NO_GZ=genome.fasta

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta

            ${RESOURCES_ARTIFACTS_BISMARK_ARTIFACT_SCRIPTS}/bismark_genome_preparation \
                --path_to_bowtie ${RESOURCES_ARTIFACTS_BOWTIE2_ARTIFACT_BINARIES}/ \
                --verbose --bowtie2 \
                ${installation_path}/

            STATUS=$?
            if [ ${STATUS} != 0 ]; then

             return ${STATUS}
            fi
            ls -ltr   ${installation_path}/

            echo "Finished indexing Bisulfite genome with bismark Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"
            return 0

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