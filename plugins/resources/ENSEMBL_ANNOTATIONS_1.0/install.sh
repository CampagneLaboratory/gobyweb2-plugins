# Installation script for ENSEMBL_ANNOTATIONS

function plugin_install_artifact {

    id=$1
    installation_path=$2

    set -e

    GOBY_JAR=${RESOURCES_GOBYWEB_SERVER_SIDE_GLOBAL_GOBY_JAR}
    case ${id} in

        'ANNOTATIONS' )
            ORGANISM=$3
            GENOME_REFERENCE_ID=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID} ENSEMBL_RELEASE=${ENSEMBL_RELEASE}"
            ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`

            BIOMART_SERVER_PREFIX=useast
            VIRTUAL_SCHEMA_NAME=default

            TOPLEVEL_IDS_FILE=${installation_path}/toplevel-ids.compact-reads
            CHROMOSOME_LIST_FILE=${installation_path}/chromosome.list.txt

            create_chromosome_list_file

            echo "SCRIPT="${RESOURCES_ENSEMBL_ANNOTATIONS_BIOMART_SCRIPT}

            ${RESOURCES_GROOVY_EXECUTABLE} \
            -cp ${RESOURCES_GOBYWEB_SERVER_SIDE_ICB_GROOVY_SUPPORT_JAR}:${RESOURCES_ARTIFACTS_SAM_JDK_JAR}/sam-jdk.jar:${JOB_DIR} \
                        ${RESOURCES_ENSEMBL_ANNOTATIONS_BIOMART_SCRIPT} \
            --virtual-schema-name ${VIRTUAL_SCHEMA_NAME} --url-prefix ${BIOMART_SERVER_PREFIX} \
            --tabix-executable ${RESOURCES_ARTIFACTS_TABIX_BINARIES}/tabix \
            --bgzip-executable ${RESOURCES_ARTIFACTS_TABIX_BINARIES}/bgzip \
            --organism ${ORG_LOWERCASE} \
            --output-folder ${installation_path}   \
            --chromosome-list-file ${CHROMOSOME_LIST_FILE} \
            --exports all
            #exon-annotations, gene-annotations, five-prime-annotations
            #ref-start-end-gene,gene-id-description,ref-start-end-gene-hgnc


            if [ "$?" -eq "0" ]; then
             # Make sure tabix indices have a more recent timestamp than tabix data files:
             sleep 3
             touch   ${installation_path}/*.tbi
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

function create_chromosome_list_file {
CHROMOSOME_LIST_FILE_TEMP=chromosome.list.txt.tmp
GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${GENOME_REFERENCE_ID}_${ENSEMBL_RELEASE}})

if [ ! -e ${TOPLEVEL_IDS_FILE} ]; then
    echo "Creating toplevel ids file"

    java -Xmx6G -jar ${GOBY_JAR} --mode fasta-to-compact -n 1 --include-identifiers  \
        --exclude-sequences -o ${TOPLEVEL_IDS_FILE} ${GENOME_DIR}/genome-toplevel.fasta.gz
    if [ "$?" != "0" ]; then
       exit 1;
    fi
fi
if [ ! -e ${CHROMOSOME_LIST_FILE} ]; then
    CHROMOSOME_LIST_FILE_TEMP=chromosome.list.txt.tmp

    java -Xmx6G -jar ${GOBY_JAR} --mode compact-to-fasta \
        --identifier-to-header \
        -i ${TOPLEVEL_IDS_FILE} -o ${CHROMOSOME_LIST_FILE_TEMP}
        if [ "$?" != "0" ]; then
               exit 1;
            fi
    grep '>' ${CHROMOSOME_LIST_FILE_TEMP} | cut -d '>' -f 2 | sort -f -n > ${CHROMOSOME_LIST_FILE}
    if [ "$?" != "0" ]; then
           exit 1;
    fi
fi


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
    return 0
}