# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    set -x
    case ${id} in

        'BINARIES' )

            VERSION="2.1.0"
            (${RESOURCES_FETCH_URL_SCRIPT} http://sourceforge.net/projects/bowtie-bio/files/bowtie2/${VERSION}/bowtie2-${VERSION}-source.zip/download  bowtie2-${VERSION}-source.zip

                unzip bowtie2-${VERSION}-source.zip
                cd bowtie2-${VERSION}/

                make
                cp bowtie2 ${installation_path}/
                cp bowtie2-build ${installation_path}/
                cp bowtie2-align ${installation_path}/
                cp bowtie2-inpect ${installation_path}/
                chmod +x ${installation_path}/*
            )
            if [ -e ${installation_path}/bowtie2 ]; then
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