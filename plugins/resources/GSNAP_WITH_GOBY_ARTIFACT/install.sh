# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    set -x
    case ${id} in

        'EXECUTABLE' )
            . ${RESOURCES_ARTIFACTS_GOBY_CPP_API_LIBRARIES}/setup.sh
            VERSION="2013-06-26"
            (${RESOURCES_FETCH_URL_SCRIPT} http://research-pub.gene.com/gmap/src/gmap-gsnap-${VERSION}.tar.gz  gmap-gsnap.tar.gz

                tar -xvf gmap-gsnap.tar.gz
                cd  gmap-${VERSION}
                ./configure --with-goby=${RESOURCES_ARTIFACTS_GOBY_CPP_API_LIBRARIES} --prefix=${installation_path}
                make
                make install

                ls -ltr

                ls -ltr ${installation_path}/
                # cp bwa ${installation_path}/bwa-icb
                #chmod +x ${installation_path}/bwa-icb
            )
            if [ -e ${installation_path}/bin/gmap_build ]; then
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

            GMAP_BUILD_EXEC=${RESOURCES_ARTIFACTS_GSNAP_WITH_GOBY_ARTIFACT_EXECUTABLE}/bin/gmap_build


            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

            #INPUT_FASTA_NO_GZ=genome.fasta

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta
            mkdir share
            ${GMAP_BUILD_EXEC} -d index -k 15 ${INPUT_FASTA_NO_GZ} -D share
            cp -r share  ${installation_path}/
            STATUS=$?
            ls -ltr
            if [ ${STATUS} != 0 ]; then

             return ${STATUS}
            fi

            if [ -e ${installation_path}/share/index/index.chromosome ]; then
                return  0
            else
                return 127
            fi

            echo "Finished indexing Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"
            exit 11
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