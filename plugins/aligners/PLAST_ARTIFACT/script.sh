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


       set -x
      # Extract the reads if a split is needed
      if [ ! -z ${SGE_TASK_ID} ] && [ "${SGE_TASK_ID}" != "undefined" ] && [ "${SGE_TASK_ID}" != "unknown" ]; then
          ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_SPLIT_STATUS} --description "Split, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, starting" --index ${CURRENT_PART} --job-type job-part
          # The reads file to process
          READS_FILE=${READS##*/}

          ls -l ${READS_FILE}
          ls -l ${READS}
          run-goby ${PLUGIN_NEED_ALIGN_JVM} reformat-compact-reads --output ${READS_FILE} \
              --start-position ${START_POSITION} --end-position ${END_POSITION} ${READS}

          dieUponError "split reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

      fi

       # Make sure the scripts we use are executable:
      chmod +x ${PLUGINS_ALIGNER_PLAST_ARTIFACT_FILES_ALIGN_BOTH_STRANDS} ${RESOURCES_PLAST_SCRIPT} ${RESOURCES_GOBY_SHELL_SCRIPT}
      export PAIRED_END_ALIGNMENT
      export JOB_DIR
      ${RESOURCES_PLAST_SCRIPT} ${JOB_DIR} ${PLUGINS_ALIGNER_PLAST_ARTIFACT_FILES_ALIGN_BOTH_STRANDS} ${READS_FILE} ${TAG}-tmp-align ${PAIRED_END_ALIGNMENT}
      goby merge-compact-alignments  ${TAG}-tmp-align -o ${OUTPUT}
      dieUponError "Aligning forward and reverse strand results failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"


}
