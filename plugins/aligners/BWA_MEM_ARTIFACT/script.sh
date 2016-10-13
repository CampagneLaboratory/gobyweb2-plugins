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
. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_align {

    OUTPUT=$1
    BASENAME=$2
    COLOR_SPACE_OPTION=""
    if [ "${COLOR_SPACE}" == "true" ]; then
     echo "Color space option is not supported"
     exit 1;
    fi
    BWA_GOBY_EXEC_PATH=${RESOURCES_ARTIFACTS_BWA07_EXECUTABLE}/bwa
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

        #run_goby ${PLUGIN_NEED_ALIGN_JVM} reformat-compact-reads --output paired-content.compact-reads \

        run_goby ${PLUGIN_NEED_ALIGN_JVM} compact-to-fasta paired-content.compact-reads -t fastq \
            --start-position ${START_POSITION} --end-position ${END_POSITION} ${READS} \
            -o read1.fq -p read2.fq
        dieUponError "compact read to fastq failed (paired end), sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} mem  -t ${BWA_GOBY_NUM_THREADS}  ${ALN_OPTIONS}  ${INDEX_DIR} read1.fq read2.fq >output.sam
        dieUponError "pair alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

    else
        # Single end alignment, native aligner
       run_goby ${PLUGIN_NEED_ALIGN_JVM} compact-to-fasta paired-content.compact-reads -t fastq \
            --start-position ${START_POSITION} --end-position ${END_POSITION} ${READS} \
            -o read1.fq
        dieUponError "compact read to fastq failed (paired end), sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} mem  -t ${BWA_GOBY_NUM_THREADS}  ${ALN_OPTIONS}  ${INDEX_DIR} read1.fq >output.sam
        dieUponError "SE alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
    fi
    run_goby ${PLUGIN_NEED_ALIGN_JVM} sam-to-compact -i output.sam -o ${OUTPUT}
    dieUponError "SAM conversion to Goby fomat failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
}
