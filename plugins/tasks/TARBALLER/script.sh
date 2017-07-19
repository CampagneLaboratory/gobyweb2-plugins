#!/bin/sh

#
. ${JOB_DIR}/constants.sh

TXT_FILES_LIST=""
JPEG_FILES_LIST=""

OUTPUT_FILE=${JOB_DIR}/output-data/out.tar.gz

function plugin_task {

     # ${FILESET_COMMAND}
     # options:  (-i|--pb-file) <pbfile> [(-d|--download-dir) <downloadDir>] [-q|--has-fileset] [-f|--fetch] [-p|--push] [-h|--help] [-l|--info] fileset1 fileset2 ... filesetN
     info "TARBALLER execution started"
     ${FILESET_COMMAND} --has-fileset TEXT
     if [ $? != 0 ]; then
       echo Input TEXTs are not available
     else
        TXT_FILES_LIST=`${FILESET_COMMAND} --fetch TEXT`
        if [ $? != 0 ]; then
            error "Failed to fecth TEXT entries ${TXT_FILES_LIST}"
            return 0
         fi
     fi                                                                     
     debug "Fetched TEXT entries ${TXT_FILES_LIST}"  
     ${FILESET_COMMAND} --has-fileset IMAGE
     if [ $? != 0 ]; then
        error "Input IMAGEs are not available"
     else
        JPEG_FILES_LIST=`${FILESET_COMMAND} --fetch IMAGE`
         if [ $? != 0 ]; then
            error "Failed to fecth IMAGE entries ${JPEG_FILES_LIST} "
            return 0
         fi
     fi
     debug "Fetched IMAGE entries ${JPEG_FILES_LIST}"

     mkdir -p "${JOB_DIR}/output-data"
     tar -zcvf ${OUTPUT_FILE} ${TXT_FILES_LIST} ${JPEG_FILES_LIST}
     REGISTERED_TAGS=`${FILESET_COMMAND} --push ARCHIVE: ${OUTPUT_FILE}`
     dieUponError "Failed to register the tar archive."

     ALL_REGISTERED_TAGS="${ALL_REGISTERED_TAGS} COMPACT_READ_FILES:[${REGISTERED_TAGS}]"
     info "The following tags were registered by this plugin: ${ALL_REGISTERED_TAGS}"
     info "TARBALLER execution completed"
}
