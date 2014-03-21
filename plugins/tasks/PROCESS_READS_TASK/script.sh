#!/bin/sh

#
. constants.sh

READ_FILES_LIST=""

function plugin_task {

     echo "fileset command: ${FILESET_COMMAND}"

     ${FILESET_COMMAND} --has-fileset UPLOADS_FILES
     if [ $? != 0 ]; then
       dieUponError "Input files are not available"

     fi

     READ_FILES_LIST=`${FILESET_COMMAND} --fetch UPLOADS_FILES`
     if [ $? != 0 ]; then
        dieUponError "Failed to fetch uploaded files ${UPLOADS_FILES}"
        echo ${UPLOADS_FILES}

     fi
     echo "Localized input files: ${UPLOADS_FILES}"

     # DO SOMETHING WITH THE FILES

     #push back the generated compact-reads:
     REGISTERED_TAGS=`${FILESET_COMMAND} --push COMPACT_READ_FILES: *.compact-reads`
     if [ $? != 0 ]; then
        dieUponError "Failed to push back the output TSV file"

     fi
     echo "PROCESS_READS registered the following FileSet instances: ${REGISTERED_TAGS}"
     return 0
}


