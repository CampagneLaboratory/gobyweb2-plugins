#!/bin/sh

#
. constants.sh

READ_FILES_LIST=""

function plugin_task {

     echo "fileset command: ${FILESET_COMMAND}"

     ${FILESET_COMMAND} --has-fileset INPUT_READS
     if [ $? != 0 ]; then
       dieUponError "Input compact reads are not available"

     fi

     READ_FILES_LIST=`${FILESET_COMMAND} --fetch INPUT_READS`
     if [ $? != 0 ]; then
        dieUponError "Failed to fecth compact reads ${READ_FILES_LIST}"
        echo ${READ_FILES_LIST}

     fi
     echo "Localized filesets ${READ_FILES_LIST}"

     java  -jar ${RESOURCES_RNASELECT_RNASELECT_TOOL} --output out.tsv ${READ_FILES_LIST}

     #push back the generated tsv
     REGISTERED_TAGS=`${FILESET_COMMAND} --push STATS: *.tsv`
     if [ $? != 0 ]; then
        dieUponError "Failed to push back the output TSV file"

     fi
     echo "RNA-select registered the following FileSet instances: ${REGISTERED_TAGS}"
     return 0
}


