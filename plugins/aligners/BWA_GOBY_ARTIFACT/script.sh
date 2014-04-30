# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true|false
# COLOR_SPACE = true|false
# READS = reads file
# START_POSITION = start index in the reads file
# END_POSITION = end index in the reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search
# ALIGNER_OPTIONS = any BWA options the end-user would like to set

function plugin_align {

    OUTPUT=$1
    BASENAME=$2
    . ${RESOURCES_ARTIFACTS_PROTOBUF_CPP_LIBRARIES}/setup.sh
    . ${RESOURCES_ARTIFACTS_GOBY_CPP_API_LIBRARIES}/setup.sh
    COLOR_SPACE_OPTION=""
    if [ "${COLOR_SPACE}" == "true" ]; then
        COLOR_SPACE_OPTION="-c"
    fi
    BWA_GOBY_EXEC_PATH=${RESOURCES_ARTIFACTS_BWA_WITH_GOBY_ARTIFACT_EXECUTABLE}/bin/bwa
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
    SAMPE_SAMSE_OPTIONS="${PLUGINS_ALIGNER_BWA_GOBY_ARTIFACT_SAMPE_SAMSE_OPTIONS}"
    ALN_OPTIONS="${PLUGINS_ALIGNER_BWA_GOBY_ARTIFACT_ALN_OPTIONS}"
    ALL_OTHER_OPTIONS="${PLUGINS_ALIGNER_BWA_GOBY_ARTIFACT_ALL_OTHER_OPTIONS}"
    INDEX_DIR=$(eval echo \${RESOURCES_ARTIFACTS_BWA_WITH_GOBY_ARTIFACT_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})/index
    BWA_GOBY_NUM_THREADS=3
    MAXIMUM_NUMBER_GAP_OPENS="${PLUGINS_ALIGNER_BWA_GOBY_ARTIFACT_MAXIMUM_NUMBER_GAP_OPENS}"
    MAXIMUM_NUMBER_GAP_EXTENSIONS="${PLUGINS_ALIGNER_BWA_GOBY_ARTIFACT_MAXIMUM_NUMBER_GAP_EXTENSIONS}"
    AMBIGUITY_THRESHOLD="${PLUGINS_ALIGNER_BWA_GOBY_ARTIFACT_AMBIGUITY_THRESHOLD}"

    GAP_OPTIONS="-o ${MAXIMUM_NUMBER_GAP_OPENS} -e ${MAXIMUM_NUMBER_GAP_EXTENSIONS} "
    SAMPLE_NAME=`basename ${READS_FILE}`
    PLATFORM_NAME=${READS_PLATFORM}
    READ_GROUPS="@RG\tID:1\tSM:${SAMPLE_NAME}\tPL:${PLATFORM_NAME}\tPU:1"
    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
        # PAIRED END alignment, native aligner
        SAI_FILE_0=${READS##*/}-0.sai
        SAI_FILE_1=${READS##*/}-1.sai

        nice ${BWA_GOBY_EXEC_PATH} aln -w 0 -t ${BWA_GOBY_NUM_THREADS} ${COLOR_SPACE_OPTION} -f ${SAI_FILE_0} -l ${INPUT_READ_LENGTH} ${ALN_OPTIONS} -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIR} ${READS_FILE} ${GAP_OPTIONS}
        dieUponError "primary alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
        nice ${BWA_GOBY_EXEC_PATH} aln -w 1 -t ${BWA_GOBY_NUM_THREADS} ${COLOR_SPACE_OPTION} -f ${SAI_FILE_1}  ${ALN_OPTIONS} -l ${INPUT_READ_LENGTH}  -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIR} ${READS_FILE}
        dieUponError "secondary alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
        # aln worked, let's sampe
        nice ${BWA_GOBY_EXEC_PATH} sampe ${COLOR_SPACE_OPTION}  -F goby -f ${OUTPUT} ${SAMPE_SAMSE_OPTIONS} -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIR} ${SAI_FILE_0} ${SAI_FILE_1} ${READS_FILE} ${READS_FILE} -r ${READ_GROUPS} -n ${AMBIGUITY_THRESHOLD}
        dieUponError "samse failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
    else
        # Single end alignment, native aligner
        SAI_FILE_0=${READS##*/}.sai
        nice ${BWA_GOBY_EXEC_PATH} aln ${COLOR_SPACE_OPTION} -t ${BWA_GOBY_NUM_THREADS} -f ${SAI_FILE_0} -l ${INPUT_READ_LENGTH} ${ALN_OPTIONS} -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIR} ${READS_FILE} ${GAP_OPTIONS}
        dieUponError "alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
        # aln worked, let's samse
        nice ${BWA_GOBY_EXEC_PATH} samse ${COLOR_SPACE_OPTION} ${SAMPE_SAMSE_OPTIONS}  -F goby -f ${OUTPUT} -x ${START_POSITION} -y ${END_POSITION} ${INDEX_DIR} ${SAI_FILE_0} ${READS_FILE} -r ${READ_GROUPS} -n ${AMBIGUITY_THRESHOLD}
        dieUponError "samse failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
    fi
}
