# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true (alignment must be paired-end)
# COLOR_SPACE = true|false
# READS = reads file
# START_POSITION = start index in the reads file
# END_POSITION = end index in the reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search

# BWA_GOBY_EXEC_PATH = path to BWA, obtained from environment.sh
# BWA_GOBY_NUM_THREADS = number of threads to run with, obtained from environment.sh

# ALIGNER_OPTIONS = any BWA options the end-user would like to set
. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_align {

    OUTPUT=$1
    BASENAME=$2

    COLOR_SPACE_OPTION=""
    if [ "${COLOR_SPACE}" == "true" ]; then
        COLOR_SPACE_OPTION="-c"
    fi

    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then

        # quick test that Goby is available:
        run-goby ${PLUGIN_NEED_ALIGN_JVM} version

        # PAIRED END alignment, native aligner
        SAI_FILE_0=${READS##*/}-0.sai
        SAI_FILE_1=${READS##*/}-1.sai
        nice ${BWA_GOBY_EXEC_PATH} aln -w 0 -t ${BWA_GOBY_NUM_THREADS} ${COLOR_SPACE_OPTION} -f ${SAI_FILE_0} -l ${INPUT_READ_LENGTH} ${ALIGNER_OPTIONS} -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIRECTORY}/${INDEX_PREFIX} ${READS_FILE}
        dieUponError "bwa aln failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} samse ${COLOR_SPACE_OPTION} -F goby -f ${OUTPUT}_0 -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIRECTORY}/${INDEX_PREFIX} ${SAI_FILE_0} ${READS_FILE}
        dieUponError "bwa samse failed for _0, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} aln -w 1 -t ${BWA_GOBY_NUM_THREADS} ${COLOR_SPACE_OPTION} -f ${SAI_FILE_1} -l ${INPUT_READ_LENGTH} ${ALIGNER_OPTIONS} -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIRECTORY}/${INDEX_PREFIX} ${READS_FILE}
        dieUponError "bwa aln failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} samse ${COLOR_SPACE_OPTION} -F goby -f ${OUTPUT}_1 -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIRECTORY}/${INDEX_PREFIX} ${SAI_FILE_1} ${READS_FILE}
        dieUponError "bwa samse failed for _1, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        # aln worked, we use the goby merge mode to combine the data into one alignment:
        run-goby ${PLUGIN_NEED_ALIGN_JVM}  merge-compact-alignments --hi-c ${OUTPUT}_0 ${OUTPUT}_1  -o ${OUTPUT}
        RETURN_STATUS=$?
    else
        # Single end alignment, native aligner
        dieUponError "bwa Hi-C requires a paired-end sample."
    fi
}
