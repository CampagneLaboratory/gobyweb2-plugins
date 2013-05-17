# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:


# READS = reads file
# START_POSITION = start index in the reads file
# END_POSITION = end index in the reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search
# REFERENCE = Top level id file for reference genome.
# ALIGNER_OPTIONS = any Last options the end-user would like to set

. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_align {

      OUTPUT=$1
      BASENAME=$2

      if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
           dieUponError "Plugin LAST_BISULFITE does not support paired-end read files, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
      fi

      if [ ! -e ${RESOURCES_LAST_BISULFITE_REVERSE_MATRIX} ]; then
       dieUponError "Last Resource is missing a required file (BISULFITE_REVERSE_MATRIX). Aborting."
      fi
      if [ ! -e ${RESOURCES_LAST_BISULFITE_FORWARD_MATRIX} ]; then
       dieUponError "Last Resource is missing a required file (BISULFITE_FORWARD_MATRIX). Aborting."
      fi
      if [ ! -e ${RESOURCES_LAST_MAP_PROBS_EXEC} ]; then
       dieUponError "Last Resource is missing a required file (MAP_PROBS_EXEC). Aborting."
      fi

      # Extract the reads if a split is needed
      if [ ! -z ${SGE_TASK_ID} ] && [ "${SGE_TASK_ID}" != "undefined" ] && [ "${SGE_TASK_ID}" != "unknown" ]; then
          ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_SPLIT_STATUS} --description "Split, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, starting" --index ${CURRENT_PART} --job-type job-part
          # The reads file to process
          READS_FILE=${READS##*/}

          goby reformat-compact-reads ${READS} --output ${READS_FILE} \
              --start-position ${START_POSITION} --end-position ${END_POSITION}

          dieUponError "split reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

      fi
      ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_ALIGN_STATUS} --description "Align, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, starting" --index ${CURRENT_PART} --job-type job-part

      # Make sure the scripts we use are executable:
      chmod +x ${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_FILES_ALIGN_BOTH_STRANDS} ${RESOURCES_PLAST_SCRIPT} ${RESOURCES_GOBY_SHELL_SCRIPT}

      ${RESOURCES_PLAST_SCRIPT} ${JOB_DIR} ${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_FILES_ALIGN_BOTH_STRANDS} ${READS_FILE} ${TAG}-tmp-align
      goby merge-compact-alignments  ${TAG}-tmp-align -o ${OUTPUT}
      dieUponError "Aligning forward and reverse strand results failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"


}
