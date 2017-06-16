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

# override the number of alignment parts to make larger chunks:
NUMBER_OF_ALIGN_PARTS=`echo $(( NUMBER_OF_ALIGN_PARTS / 20 ))`
if [ ${NUMBER_OF_ALIGN_PARTS} -lt 1 ]; then
    NUMBER_OF_ALIGN_PARTS=1
fi

function plugin_align {

      OUTPUT=$1
      BASENAME=$2

      COLOR_SPACE_OPTION=""
      if [ "${COLOR_SPACE}" == "true" ]; then
          COLOR_SPACE_OPTION=" --color-space "
      fi

        NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
        NUM_THREADS=$((${NUM_THREADS} - 2))
        if [ ${NUM_THREADS} -lt 2 ]; then
           # Make sure we use at least two threads:
           NUM_THREADS=2
        fi


    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    INDEX_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_LAST_INDEX_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}}/)
    TOPLEVEL_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_LAST_INDEX_TOPLEVEL_IDS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}}/)

    #make sure index exists
    if [ ! -e ${INDEX_DIRECTORY}/index.prj ]; then
        failThisLine
      #	dieUponError "last index could not be found"
    fi
       set -x
      # Extract the reads if a split is needed
      if [ ! -z ${SGE_TASK_ID} ] && [ "${SGE_TASK_ID}" != "undefined" ] && [ "${SGE_TASK_ID}" != "unknown" ]; then
          ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_SPLIT_STATUS} --description "Split, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, starting" --index ${CURRENT_PART} --job-type job-part
          # The reads file to process
          READS_FILE=${READS##*/}

          ls -l ${READS_FILE}
          ls -l ${READS}
          run_goby ${PLUGIN_NEED_ALIGN_JVM} reformat-compact-reads --output ${READS_FILE} \
              --start-position ${START_POSITION} --end-position ${END_POSITION} ${READS}

          dieUponError "split reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
      fi
    TEMP_FILENAME="last-alignment"
    ${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/bin/lastal -P ${NUM_THREADS} -s2 -Q1 -d${PLUGINS_ALIGNER_LAST_NATIVE_PARALLEL_ARTIFACT_D} \
        -e${PLUGINS_ALIGNER_LAST_NATIVE_PARALLEL_ARTIFACT_E} ${INDEX_DIRECTORY}/index ${READS_FASTQ} > ${TEMP_FILENAME}.maf
    dieUponError "last could not align reads"
    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
        ${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}//bin/lastal -P ${NUM_THREADS} -s2 -Q1 -d${PLUGINS_ALIGNER_LAST_NATIVE_PARALLEL_ARTIFACT_D} \
           -e${PLUGINS_ALIGNER_LAST_NATIVE_PARALLEL_ARTIFACT_E} ${INDEX_DIRECTORY}/index ${PAIRS_FASTQ} > ${TEMP_FILENAME}-pairs.maf
        dieUponError "last could not align paired reads"
        ${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/bin/last-pair-probs ${TEMP_FILENAME}.maf ${TEMP_FILENAME}-pairs.maf > ${TEMP_FILENAME}-2.maf
        if [ $? != 0 ]; then
           cat ${TEMP_FILENAME}.maf ${TEMP_FILENAME}-pairs.maf > ${TEMP_FILENAME}-2.maf
        fi
        dieUponError "last could not last-pair-probs"
    else
        cat ${TEMP_FILENAME}.maf |  ${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/bin/last-map-probs -s${PLUGINS_ALIGNER_LAST_NATIVE_PARALLEL_ARTIFACT_S} > ${TEMP_FILENAME}-2.maf
        dieUponError "last could not last-map-probs"
    fi
    REFERENCE=${TOPLEVEL_DIRECTORY}/toplevel-ids.compact-reads

    java -Xmx${PLUGIN_NEED_ALIGN_JVM} -Dlog4j.debug=true -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                        -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                        -jar ${RESOURCES_GOBY_GOBY_JAR} \
                        --mode last-to-compact -i ${TEMP_FILENAME}-2.maf -o ${OUTPUT} --third-party-input true \
                        --only-maf -q ${READS} -t ${REFERENCE} --quality-filter-parameters threshold=1.0
     dieUponError "Aligning forward and reverse strand results failed. Goby could not convert maf to compact alignment format, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
}
