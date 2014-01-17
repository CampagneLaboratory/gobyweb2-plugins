<<<<<<< HEAD
set -x
FULL_READS_INPUT=$1
PAIRED_END_ALIGNMENT=$2
READS_FASTQ=$3

if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
    PAIRS_FASTQ=$4
    shift 1
fi

TEMP_FILENAME=$4
OUTPUT=$5
JOB_DIR=$6

=======
#!/bin/bash
FULL_READS_INPUT=$1
READS_FASTQ=$2
TEMP_FILENAME=$3
OUTPUT=$4
JOB_DIR=$5
>>>>>>> 1809ee0ac8d84a4e697cfb6f114e4e36f22bd503
# Grab the variables and functions we need:
. ${JOB_DIR}/artifacts.sh
expose_artifact_environment_variables

<<<<<<< HEAD


=======
>>>>>>> 1809ee0ac8d84a4e697cfb6f114e4e36f22bd503
    RESOURCES_LAST_EXEC_PATH=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/bin/lastal
    RESOURCES_LAST_MERGE_BATCHES_EXEC=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/scripts/last-merge-batches.py
    RESOURCES_LAST_MAP_PROBS_EXEC=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/scripts/last-map-probs.py


    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    INDEX_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_LAST_INDEX_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}}/)
    TOPLEVEL_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_LAST_INDEX_TOPLEVEL_IDS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}}/)

    #make sure index exists
    if [ ! -e ${INDEX_DIRECTORY}/index.prj ]; then
        failThisLine
<<<<<<< HEAD
     # 	dieUponError "last index could not be found"
=======
      	dieUponError "last index could not be found"
>>>>>>> 1809ee0ac8d84a4e697cfb6f114e4e36f22bd503
    fi

    ${RESOURCES_LAST_EXEC_PATH} -v -s1 -Q1 -d${PLUGINS_ALIGNER_PLAST_ARTIFACT_D} \
        -e${PLUGINS_ALIGNER_PLAST_ARTIFACT_E} ${INDEX_DIRECTORY}/index ${READS_FASTQ} -o ${TEMP_FILENAME}.maf
<<<<<<< HEAD
    # dieUponError "last could not align reads"

    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
        ${RESOURCES_LAST_EXEC_PATH} -v -s1 -Q1 -d${PLUGINS_ALIGNER_PLAST_ARTIFACT_D} \
           -e${PLUGINS_ALIGNER_PLAST_ARTIFACT_E} ${INDEX_DIRECTORY}/index ${PAIRS_FASTQ} -o ${TEMP_FILENAME}-pairs.maf
      #  dieUponError "last could not align paired reads"
        ${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/scripts/last-pair-probs.py ${TEMP_FILENAME}.maf ${TEMP_FILENAME}-pairs.maf > ${TEMP_FILENAME}-2.maf
      #  dieUponError "last could not last-pair-probs.py"
    else
        cat ${TEMP_FILENAME}.maf | ${RESOURCES_LAST_MAP_PROBS_EXEC} -s${PLUGINS_ALIGNER_PLAST_ARTIFACT_S} > ${TEMP_FILENAME}-2.maf
     #   dieUponError "last could not last-map-probs.py"
    fi

=======

    cat ${TEMP_FILENAME}.maf | ${RESOURCES_LAST_MAP_PROBS_EXEC} -s${PLUGINS_ALIGNER_PLAST_ARTIFACT_S} > ${TEMP_FILENAME}-2.maf
>>>>>>> 1809ee0ac8d84a4e697cfb6f114e4e36f22bd503
    REFERENCE=${TOPLEVEL_DIRECTORY}/toplevel-ids.compact-reads

    java -Xmx${PLUGIN_NEED_ALIGN_JVM} -Dlog4j.debug=true -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                        -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                        -jar ${RESOURCES_GOBY_GOBY_JAR} \
                        --mode last-to-compact -i ${TEMP_FILENAME}-2.maf -o ${OUTPUT} --third-party-input true \
                        --only-maf -q ${FULL_READS_INPUT} -t ${REFERENCE} --quality-filter-parameters threshold=1.0
<<<<<<< HEAD
    #dieUponError "goby could not convert maf to compact alignment format"
=======
>>>>>>> 1809ee0ac8d84a4e697cfb6f114e4e36f22bd503
