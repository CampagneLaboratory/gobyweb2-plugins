#!/bin/sh

#
. ${JOB_DIR}/constants.sh

READ_FILES_LIST=""

function plugin_task {
     set -xv
     cd ${TMPDIR}
     echo "fileset command: ${FILESET_COMMAND}"

     ${FILESET_COMMAND} --has-fileset INPUT_VCF
     if [ $? -eq 0 ]; then

        VCF_INPUT=`${FILESET_COMMAND} --fetch INPUT_VCF`
        dieUponError "Failed to fetch input vcf files ${VCF_INPUT}"
        echo "${VCF_INPUT}"
        cp "${VCF_INPUT}" `basename ${VCF_INPUT}`
        export QUEUE_WRITER
        ${RESOURCES_GROOVY_EXECUTABLE} -cp ${GOBY_DIR}:${RESOURCES_GOBYWEB_SERVER_SIDE_GLOBAL_GOBY_JAR}:${RESOURCES_GOBYWEB_SERVER_SIDE_ICB_GROOVY_SUPPORT_JAR} \
        ${RESOURCES_GOBYWEB_SERVER_SIDE_TSV_VCF_TO_SQLITE} \
        --job-start-status "${JOB_START_STATUS}" \
        --queue-writer-prefix-variable QUEUE_WRITER \
        --export-format lucene   *.vcf
        jobDieUponError "failed to convert VCF results to Lucene Table"
     fi

     ${FILESET_COMMAND} --has-fileset INPUT_TSV
     if [ $? -eq 0 ]; then

        TSV_INPUT=`${FILESET_COMMAND} --fetch INPUT_TSV`
        dieUponError "Failed to fetch input tsv files ${TSV_INPUT}"
        echo "${TSV_INPUT}"
        cp "${TSV_INPUT}" `basename ${TSV_INPUT}`
        export QUEUE_WRITER
        ${RESOURCES_GROOVY_EXECUTABLE} -cp ${GOBY_DIR}:${RESOURCES_GOBYWEB_SERVER_SIDE_GLOBAL_GOBY_JAR}:${RESOURCES_GOBYWEB_SERVER_SIDE_ICB_GROOVY_SUPPORT_JAR} \
        ${RESOURCES_GOBYWEB_SERVER_SIDE_TSV_VCF_TO_SQLITE} \
        --job-start-status "${JOB_START_STATUS}" \
        --queue-writer-prefix-variable QUEUE_WRITER \
        --export-format lucene  *.tsv
        jobDieUponError "failed to convert TSV results to Lucene Table"
     fi

     REGISTERED_TAGS=`${FILESET_COMMAND} --push LUCENE_TABLE: *.lucene.index `
     dieUponError "Failed to push back the Lucene Table Result"

     ALL_REGISTERED_TAGS="${REGISTERED_TAGS} "
     info "The following tags were registered by this plugin: ${ALL_REGISTERED_TAGS}"
     exit 0
}


