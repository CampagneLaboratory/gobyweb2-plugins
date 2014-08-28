#!/bin/sh

#
. ${JOB_DIR}/constants.sh

READ_FILES_LIST=""

function plugin_task {
     set -xv

     GOBY_JAR_DIR="goby"
     export JAVA_OPTS="${GRID_JVM_FLAGS}"

    # source setup of R and RJAVA:
     . ${RESOURCES_ARTIFACTS_R_BINARIES}/setup.sh
     . ${RESOURCES_ARTIFACTS_RJAVA_BINARIES}/setup.sh
     . ${RESOURCES_GOBY_SHELL_SCRIPT}

     env
     # Note that we override the grid jvm flags to request only 4Gb:
     run_goby ${PLUGIN_NEED_PROCESS_JVM} test-r-connection 2>&1 | tee log.txt
     STATUS=$?

     # push back the output stats:
     LOG_REGISTERED_TAGS=`${FILESET_COMMAND} --push EXECUTION_LOG: log.txt `
     dieUponError "Failed to push back the execution log"

     ALL_REGISTERED_TAGS="${LOG_REGISTERED_TAGS}"
     echo "The following tags were registered by this plugin: ${ALL_REGISTERED_TAGS}"
     return $STATUS
}


