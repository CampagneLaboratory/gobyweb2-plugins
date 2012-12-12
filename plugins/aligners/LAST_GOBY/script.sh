# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true|false
# COLOR_SPACE = true|false
# READS = reads file
# START_POSITION = start index in the reads file
# END_POSITION = end index in the reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search

# ALIGNER_OPTIONS = any Last options the end-user would like to set

# Please note that Goby must be configured with appropriate path to Last aligner executable.
. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_align {

      OUTPUT=$1
      BASENAME=$2

      COLOR_SPACE_OPTION=""
      if [ "${COLOR_SPACE}" == "true" ]; then
          COLOR_SPACE_OPTION=" --color-space "
      fi

      # Extract the reads if a split is needed
      if [ ! -z ${SGE_TASK_ID} ] && [ "${SGE_TASK_ID}" != "undefined" ] && [ "${SGE_TASK_ID}" != "unknown" ]; then
          ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_SPLIT_STATUS} --description "Split, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, starting" --index ${CURRENT_PART} --job-type job-part
          # The reads file to process
          READS_FILE=${READS##*/}


          run-goby ${PLUGIN_NEED_ALIGN_JVM} reformat-compact-reads --output ${READS_FILE} \
              --start-position ${START_POSITION} --end-position ${END_POSITION} ${READS}

          dieUponError "split reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

      fi

       # Override the path to last to use the resource
       dirname=`dirname ${RESOURCES_LAST_EXEC_PATH}`
       echo "executables.path.last = ${dirname}" >>${TMPDIR}/goby.properties

       ALIGNER_OPTIONS_COMPLETE="matchQuality="${PLUGINS_ALIGNER_LAST_GOBY_MATCH_QUALITY}","\
"maxGapsAllowed="${PLUGINS_ALIGNER_LAST_GOBY_MAX_GAPS_ALLOWED}","\
"gapOpeningCost="${PLUGINS_ALIGNER_LAST_GOBY_GAP_EXISTENCE_COST}","\
"gapExtensionCost="${PLUGINS_ALIGNER_LAST_GOBY_GAP_EXTENSION_COST}","\
${ALIGNER_OPTIONS}

       # This Goby wrapper detects automatically if the reads file is paired end:
       run-goby ${PLUGIN_NEED_ALIGN_JVM} align --reference ${REFERENCE} --aligner last ${COLOR_SPACE_OPTION} --search \
           --ambiguity-threshold ${PLUGINS_ALIGNER_LAST_GOBY_AMBIGUITY_THRESHOLD} ${PLUGINS_ALIGNER_LAST_GOBY_ALL_OTHER_OPTIONS} \
           --database-name ${INDEX_PREFIX} --database-directory ${INDEX_DIRECTORY} \
           ${ALIGNER_OPTIONS} --reads ${READS_FILE} --basename ${OUTPUT} --options ${ALIGNER_OPTIONS_COMPLETE}

}
