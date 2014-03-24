#!/bin/sh

#
. constants.sh

READ_FILES_LIST=""

function plugin_task {
     set -xv
     echo "fileset command: ${FILESET_COMMAND}"

     ${FILESET_COMMAND} --has-fileset UPLOADS_FILES
     if [ $? != 0 ]; then
       dieUponError "Input files are not available"

     fi

     ${FILESET_COMMAND} --has-fileset UPLOAD_MERGE_PLAN
     if [ $? != 0 ]; then
            dieUponError "Merge Plan is not available"
     fi

     READ_FILES_LIST=`${FILESET_COMMAND} --fetch UPLOADS_FILES`
     if [ $? != 0 ]; then
        dieUponError "Failed to fetch uploaded files ${READ_FILES_LIST}"
        echo ${READ_FILES_LIST}

     fi

    MERGE_PLAN_FILE=`${FILESET_COMMAND} --fetch UPLOAD_MERGE_PLAN`
    if [ $? != 0 ]; then
     dieUponError "Failed to fetch merge plan ${MERGE_PLAN_FILE}"
     echo ${MERGE_PLAN_FILE}
    fi

    mkdir CONVERTED

    echo "Localized input files: ${READ_FILES_LIST}"
    READS_FIRST_FILE_TAG=`echo ${READ_FILES_LIST} | awk '{ print $1}' `
    CURRENT_RETRY=1
    MAX_RETRIES=4
    RETURN_STATUS=1
    GOBY_JAR_DIR="goby"
    #
    # I (i.e., Kevin Dorff, early versions) have, on occasion, gotten errors from running
    # this such as stale nfs and java launch errors.
    # Here I will repeat up to MAX_RETRIES times before actually failing
    #
    until [ ${RETURN_STATUS} -eq 0 ] || [ ${CURRENT_RETRY} -gt ${MAX_RETRIES} ]; do
        if [ ${CURRENT_RETRY} -gt 1 ]; then
            ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_START_STATUS} --description "Previous attempt to process reads failed. This will be retried up to ${MAX_RETRIES} times." --index 0 --job-type job-part
            # Last execution failed. Wait 60-120 seconds and try again.
            WAIT_TIME=$(( 60 + ($RANDOM % 60) ))
            sleep $WAIT_TIME
        fi
        export QUEUE_WRITER
        chmod +x ${RESOURCES_GROOVY_EXECUTABLE}
        ${RESOURCES_GROOVY_EXECUTABLE} -cp ${GOBY_DIR}:${RESOURCES_GOBYWEB_SERVER_SIDE_GLOBAL_GOBY_JAR}:${RESOURCES_GOBYWEB_SERVER_SIDE_ICB_GROOVY_SUPPORT_JAR} \
            ${RESOURCES_PROCESS_READS_PROCESS_SAMPLES} \
            --jvm-flags "${GRID_JVM_FLAGS}" \
            --goby-jar-dir ${GOBY_JAR_DIR} \
            --cluster-reads-dir ./CONVERTED ${READ_FILES_LIST} \
            --sample-tag ${TAG} \
            --first-file-tag ${PLUGINS_TASK_PROCESS_READS_TASK_TAG} \
            --quality-encoding ${PLUGINS_TASK_PROCESS_READS_TASK_QUALITY_ENCODING} \
            --platform ${PLUGINS_TASK_PROCESS_READS_TASK_READS_PLATFORM} \
            --sample-name ${PLUGINS_TASK_PROCESS_READS_TASK_SAMPLE_NAME} \
            --color-space ${PLUGINS_TASK_PROCESS_READS_TASK_READS_COLOR_SPACE} \
            --ssh-prefix ${WEB_SERVER_SSH_PREFIX} \
            --web-files-dir REMOVE_THIS_OPTION \
            --merge-plan-filename ${MERGE_PLAN_FILE} \
            --queue-writer-prefix-variable QUEUE_WRITER \
            --job-start-status "${JOB_START_STATUS}" \
            --work-dir ${TMPDIR} \
            --output-stats output-stats.properties \
            ${WEB_SAMPLE_FILES}
        RETURN_STATUS=$?
        (( CURRENT_RETRY++ ))
    done

    if [ ! $RETURN_STATUS -eq 0 ]; then
        ls -lat
        ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Processing of sample on cluster failed, failure code is ${RETURN_STATUS}." --index 0 --job-type job
        exit
    fi
    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_COMPLETED_STATUS} --description "Processing of sample on cluster completed" --index 0 --job-type job

     # DO SOMETHING WITH THE FILES

     #push back the generated compact-reads:
     REGISTERED_TAGS=`${FILESET_COMMAND} --push COMPACT_READ_FILES: CONVERTED/*.compact-reads`
     if [ $? != 0 ]; then
        dieUponError "Failed to push back the output TSV file"
     fi
     echo "PROCESS_READS registered the following FileSet instances: ${REGISTERED_TAGS}"

     #push back the output stats:
     REGISTERED_TAGS=`${FILESET_COMMAND} --push OUTPUT_STATS: output-stats.properties `
     if [ $? != 0 ]; then
        dieUponError "Failed to push back the reads statistics properties file"
     fi
     echo "Read statistics registered the following FileSet instances: ${REGISTERED_TAGS}"

     return 0
}


