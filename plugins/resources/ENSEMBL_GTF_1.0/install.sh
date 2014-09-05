    # Installation script for ENSEMBL_GENOMES
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


      'ANNOTATIONS' )
            ORGANISM=$3
            GENOME_REFERENCE_ID=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID} ENSEMBL_RELEASE=${ENSEMBL_RELEASE}"
            ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`
            # get genome from local NFS, or wget from ensembl servers if not present
            BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
            if [ "${BUILD_NUMBER}" = "GRCh37" -o "${BUILD_NUMBER}" = "1000GENOMES" ]; then
                # For human, use the compatible 1000g assembly instead of the Ensembl build: (coordinates are compatible),
                # see http://www.1000genomes.org/category/assembly
                ${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp.ensembl.org/pub/release-74/gtf/homo_sapiens/Homo_sapiens.GRCh37.74.gtf.gz human.dna.toplevel.gtf.gz
            else
                ${RESOURCES_FETCH_URL_SCRIPT_PATTERN} "ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/gtf/${ORG_LOWERCASE}/" '.gtf.gz'
            fi
            gunzip *.gtf.gz
            cp *.gtf ${installation_path}/genome.gtf

            if [ -e ${installation_path}/genome.gtf ]; then
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