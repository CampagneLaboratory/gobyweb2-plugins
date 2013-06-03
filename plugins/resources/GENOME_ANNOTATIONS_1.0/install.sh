# Installation script for genomic annotations

function plugin_install_artifact {

    id=$1
    installation_path=$2

    set -e

    case ${id} in

        'ANNOTATIONS' )
            ORGANISM=$3
            GENOME_REFERENCE_ID=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID} ENSEMBL_RELEASE=${ENSEMBL_RELEASE}"
            ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`

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
            --exports exon-annotations,ref-start-end-gene,gene-id-description,ref-start-end-gene-hgnc


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
