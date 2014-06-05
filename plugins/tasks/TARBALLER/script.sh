#!/bin/sh

#
. constants.sh

TXT_FILES_LIST=""
JPEG_FILES_LIST=""

OUTPUT_FILE=${JOB_DIR}/output-data/out.tar

function plugin_task {

     # ${FILESET_COMMAND}
     # options:  (-i|--pb-file) <pbfile> [(-d|--download-dir) <downloadDir>] [-q|--has-fileset] [-f|--fetch] [-p|--push] [-h|--help] [-l|--info] fileset1 fileset2 ... filesetN

     ${FILESET_COMMAND} --has-fileset TEXT
     if [ $? == 0 ]; then
       echo Input TEXTs are not available
     else
        TXT_FILES_LIST=`${FILESET_COMMAND} --fetch TEXT`
        if [ $? == 0 ]; then
            echo Failed to fecth TEXT entries
            echo ${TXT_FILES_LIST}
            return 0
         fi
     fi

     ${FILESET_COMMAND} --has-fileset JPEG
     if [ $? == 0 ]; then
       echo Input JPEGs are not available
     else
        JPEG_FILES_LIST=`${FILESET_COMMAND} --fetch JPEG.*`
         if [ $? == 0 ]; then
            echo Failed to fecth JPEG entries
            echo ${JPEG_FILES_LIST}
            return 0
         fi
     fi

     echo "Localized filesets ${TXT_FILES_LIST} ${JPEG_FILES_LIST}"

     mkdir "${JOB_DIR}/output-data"
     tar -cvf ${OUTPUT_FILE} ${TXT_FILES_LIST} ${JPEG_FILES_LIST}

     #will be replaced by a request to file-manager
     #scp ${OUTPUT_FILE} ${WEB_SERVER_SSH_PREFIX}:${RESULTS_WEB_DIR}
}
