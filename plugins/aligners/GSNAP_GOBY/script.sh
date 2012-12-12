# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true|false
# READS = reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search

# ${RESOURCES_ILLUMINA_ADAPTERS_FILE_PATH} = path to adapters.txt, obtained from the ILLUMINA_ADAPTERS resource
# ${RESOURCES_GSNAP_WITH_GOBY_EXEC_PATH} = path to gsnap, obtained from the GSNAP_GOBY resource

# ALIGNER_OPTIONS = any GSNAP options the end-user would like to set

. ${RESOURCES_EXTRACT_NONMATCHED_SHELL_SCRIPT}

function plugin_align {

     OUTPUT=$1
     BASENAME=$2
     # set the number of threads to the number of cores available on the server:
     NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
     ALIGNER_OPTIONS="${ALIGNER_OPTIONS} -t ${NUM_THREADS}"

     SPLICED_OPTION=""
     if [ "${PLUGINS_ALIGNER_GSNAP_GOBY_SPLICED_ALIGNMENT}" == "spliced" ]; then
        SPLICED_OPTION="-s ${GSNAP_SPLICE_FILE}"
     fi

     BISULFITE_OPTION=""
     if [ "${BISULFITE_SAMPLE}" == "true" ]; then

         goby reformat-compact-reads  --start-position=${START_POSITION} --end-position=${END_POSITION}  ${READS_FILE} -o small-reads.compact-reads
         dieUponError "reformat reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

         # GSNAP version 2011-03-11 and newer, for older use -C

         STRANDNESS="${PLUGINS_ALIGNER_GSNAP_GOBY_STRANDNESS}"
         BISULFITE_OPTION=" --mode "cmet-${STRANDNESS}" -m 1 -i 100 --terminal-threshold=100    "

         # Trim the reads if they are bisulfite.
         goby trim  -i small-reads.compact-reads -o small-reads-trimmed.compact-reads --complement -a  ${RESOURCES_ILLUMINA_ADAPTERS_FILE_PATH}  --min-left-length 4
         dieUponError "trim reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

         WINDOW_OPTIONS=" "
         READ_FILE_SMALL=small-reads-trimmed.compact-reads
     else
         WINDOW_OPTIONS=" --creads-window-start=${START_POSITION} --creads-window-end=${END_POSITION}  "
         READ_FILE_SMALL=" ${READS_FILE} "
     fi


     if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
         # PAIRED END alignment, native aligner
         nice ${RESOURCES_GSNAP_WITH_GOBY_EXEC_PATH} ${WINDOW_OPTIONS} -B 4 ${SPLICED_OPTION} ${BISULFITE_OPTION} ${ALIGNER_OPTIONS} ${PLUGINS_ALIGNER_GSNAP_GOBY_ALL_OTHER_OPTIONS} -A goby --goby-output="${OUTPUT}" -D ${INDEX_DIRECTORY} -d ${INDEX_PREFIX} -o ${PAIRED_END_DIRECTIONS} ${READ_FILE_SMALL}
         dieUponError "GSNAP alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

     else
         # Single end alignment, native aligner
         nice ${RESOURCES_GSNAP_WITH_GOBY_EXEC_PATH} ${WINDOW_OPTIONS} -B 4 ${SPLICED_OPTION} ${BISULFITE_OPTION} ${ALIGNER_OPTIONS} ${PLUGINS_ALIGNER_GSNAP_GOBY_ALL_OTHER_OPTIONS}  -A goby --goby-output="${OUTPUT}" -D ${INDEX_DIRECTORY} -d ${INDEX_PREFIX} ${READ_FILE_SMALL}
         dieUponError "GSNAP alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
     fi



#extra variables:

#RESULT_DIR= directory on shared filesystem, send output files to $RESULT_DIR/split-results
#CURRENT_PART= unique id associated with this part of the job


#eventually need to change xml to output the final file, but that comes later
#also make it able to handle errors

     if [ "${PLUGINS_ALIGNER_GSNAP_GOBY_NON_MATCHING}" == "true" ]; then
     
     	 extract_unmatched_in_plugin_align ${READS} ${OUTPUT} ${NUMBER_OF_PARTS} ${CURRENT_PART} ${SGE_O_WORKDIR}/split-results
     fi
}

# This function is called after the alignment slices have been combined into one final output.
# It is called with three arguments, the basename of the alignment (present in the directory where this function is
# invoked), the reads filename (full path), the tag useful to create an output.

function plugin_alignment_combine {
    TAG=$1
    READS=$2
    BASENAME=$3

    if [ "${PLUGINS_ALIGNER_GSNAP_GOBY_NON_MATCHING}" == "true" ]; then

        combine_splits  ${SGE_O_WORKDIR}/split-results ${RESULT_DIR}/ ${BASENAME}

    fi

}
