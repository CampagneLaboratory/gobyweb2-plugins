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
. ${RESOURCES_GOBY3_SHELL_SCRIPT}

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

    INDEX_DIR=$(eval echo \${RESOURCES_ARTIFACTS_BWA07_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})/index

    NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
    NUM_THREADS=$((${NUM_THREADS} - 2))
    if [ ${NUM_THREADS} -lt 2 ]; then
       # Make sure we use at least two threads:
       NUM_THREADS=2
    fi
    AMBIGUITY_THRESHOLD="${PLUGINS_ALIGNER_BWA07_AMBIGUITY_THRESHOLD}"
    MEM_OPTIONS="${PLUGINS_ALIGNER_BWA07_MEM_OPTIONS}"

    SAMPLE_NAME=`basename ${READS_FILE}`
    PLATFORM_NAME=${READS_PLATFORM}
    READ_GROUPS="@RG\tID:1\tSM:${SAMPLE_NAME}\tPL:${PLATFORM_NAME}\tPU:1"
    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
        # PAIRED END alignment, native aligner

        run_goby ${PLUGIN_NEED_ALIGN_JVM} compact-to-fasta -i ${READS} -t fastq \
            --start-position ${START_POSITION} --end-position ${END_POSITION}  \
            -o read1.fq -p read2.fq
        dieUponError "compact read to fastq failed (paired end), sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} mem  -t ${NUM_THREADS}  ${MEM_OPTIONS}  ${INDEX_DIR} read1.fq read2.fq >output.sam
        dieUponError "pair alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

    else
        # Single end alignment, native aligner
       run_goby ${PLUGIN_NEED_ALIGN_JVM} compact-to-fasta -i ${READS} -t fastq \
            --start-position ${START_POSITION} --end-position ${END_POSITION}  \
            -o read1.fq
        dieUponError "compact read to fastq failed (paired end), sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        nice ${BWA_GOBY_EXEC_PATH} mem  -t ${NUM_THREADS}  ${MEM_OPTIONS}  ${INDEX_DIR} read1.fq >output.sam
        dieUponError "SE alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
    fi
    set -x
     run_goby ${PLUGIN_NEED_ALIGN_JVM} sam-to-compact -i output.sam -o ${OUTPUT} --genome ${GENOME_DIR}/random-access-genome

     INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
     # ADD MD tags to the sam file with samtools: NB: we sort the BAM file because calmd is terribly slow on non-sorted input
     ${RESOURCES_SAMTOOLS_EXEC_PATH} view -S -b -u Aligned.out.sam |${RESOURCES_SAMTOOLS_EXEC_PATH} sort - sam_sorted
     ${RESOURCES_SAMTOOLS_EXEC_PATH} calmd -u sam_sorted.bam ${INDEXED_GENOME_DIR}/*toplevel.fasta >Aligned.out.bam
    # cp Aligned.out.bam ${JOB_DIR}/output-${CURRENT_PART}.bam

    export GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_GOBY_INDEXED_GENOMES_SEQUENCE_CACHE_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
    run_goby ${PLUGIN_NEED_ALIGN_JVM} sam-to-compact -i Aligned.out.bam -o ${OUTPUT}  --genome ${GENOME_DIR}/random-access-genome --read-names-are-query-indices
    dieUponError "SAM conversion to Goby fomat failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
}
