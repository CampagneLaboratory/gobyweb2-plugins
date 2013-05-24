# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true|false
# COLOR_SPACE = true|false
# READS = reads file
# START_POSITION = start index in the reads file
# END_POSITION = end index in the reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search

# BWA_GOBY_EXEC_PATH = path to BWA, obtained from environment.sh
# BWA_GOBY_NUM_THREADS = number of threads to run with, obtained from environment.sh

# ALIGNER_OPTIONS = any BWA options the end-user would like to set

function plugin_align {

    OUTPUT=$1
    BASENAME=$2

    COLOR_SPACE_OPTION=""
    if [ "${COLOR_SPACE}" == "true" ]; then
        COLOR_SPACE_OPTION="-c"
    fi
    BWA_GOBY_EXEC_PATH=${RESOURCES_ARTIFACTS_BWA_WITH_GOBY_ARTIFACT_EXECUTABLE}/bin/bwa
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    INDEX_DIR=$(eval echo \${RESOURCES_ARTIFACTS_BWA_WITH_GOBY_ARTIFACT_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})/index
    ALIGNER_OPTIONS="${PLUGINS_ALIGNER_BWA_BAM_ARTIFACT_ALIGNER_OPTIONS}"
    ALL_OTHER_OPTIONS="${PLUGINS_ALIGNER_BWA_BAM_ARTIFACT_ALL_OTHER_OPTIONS}"

    set -x
    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
        # PAIRED END alignment, native aligner
        SAI_FILE_0=${READS##*/}-0.sai
        SAI_FILE_1=${READS##*/}-1.sai

        nice ${BWA_GOBY_EXEC_PATH} aln -w 0 -t ${BWA_GOBY_NUM_THREADS} ${COLOR_SPACE_OPTION} -f ${SAI_FILE_0} -l ${INPUT_READ_LENGTH} ${ALIGNER_OPTIONS}  ${INDEX_DIR} ${READS}
        RETURN_STATUS=$?
        if [ $RETURN_STATUS -eq 0 ]; then
            nice ${BWA_GOBY_EXEC_PATH} aln -w 1 -t ${BWA_GOBY_NUM_THREADS} ${COLOR_SPACE_OPTION} -f ${SAI_FILE_1} -l ${INPUT_READ_LENGTH} ${ALIGNER_OPTIONS} ${INDEX_DIR} ${READS}
            RETURN_STATUS=$?
            if [ $RETURN_STATUS -eq 0 ]; then
                # aln worked, let's sampe
                nice ${BWA_GOBY_EXEC_PATH} sampe ${COLOR_SPACE_OPTION}  ${ALL_OTHER_OPTIONS}  -F sam -f ${OUTPUT} ${INDEX_DIR} ${SAI_FILE_0} ${SAI_FILE_1} ${READS} ${READS}
                RETURN_STATUS=$?
            fi
        fi
    else
        # Single end alignment, native aligner
        SAI_FILE_0=${READS##*/}.sai
        nice ${BWA_GOBY_EXEC_PATH} aln ${COLOR_SPACE_OPTION} -t ${BWA_GOBY_NUM_THREADS} -f ${SAI_FILE_0} -l ${INPUT_READ_LENGTH} ${ALIGNER_OPTIONS} ${INDEX_DIR} ${READS}
        RETURN_STATUS=$?
        if [ $RETURN_STATUS -eq 0 ]; then
            # aln worked, let's samse
            nice ${BWA_GOBY_EXEC_PATH} samse ${COLOR_SPACE_OPTION}  ${ALL_OTHER_OPTIONS}  -F sam -f ${OUTPUT} ${INDEX_DIR} ${SAI_FILE_0} ${READS}
            RETURN_STATUS=$?
        fi
    fi
     # aln worked, let's convert to BAM and sort on the fly:

    nice ${RESOURCES_SAMTOOLS_EXEC_PATH}  view -uS ${OUTPUT}  | ${RESOURCES_SAMTOOLS_EXEC_PATH}  sort - ${BASENAME}
    dieUponError "samtools view|sort step failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

    # sort worked. We index the BAM file. If this works, the return code will be 0, indicating no problem with plugin_align
    nice ${RESOURCES_SAMTOOLS_EXEC_PATH} index ${BASENAME}.bam
    dieUponError "samtools index step failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

}
