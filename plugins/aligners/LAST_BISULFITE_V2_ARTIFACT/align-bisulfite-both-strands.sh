#!/bin/bash

[ $# -eq 6 ] && shift
FULL_READS_INPUT=$1
READS_FASTQ=$2
TEMP_FILENAME1=$3
OUTPUT=$4
JOB_DIR=$5
# Grab the variables and functions we need:
. ${JOB_DIR}/constants.sh
. ${JOB_DIR}/auto-options.sh

    RESOURCES_LAST_EXEC_PATH=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/bin/lastal
    RESOURCES_LAST_MERGE_BATCHES_EXEC=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/scripts/last-merge-batches.py
    RESOURCES_LAST_MAP_PROBS_EXEC=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/scripts/last-map-probs.py
    RESOURCES_LAST_BISULFITE_FORWARD_MATRIX=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/examples/bisulfite_f.mat
    RESOURCES_LAST_BISULFITE_REVERSE_MATRIX=${RESOURCES_ARTIFACTS_LAST_ARTIFACT_BINARIES}/examples/bisulfite_r.mat



    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    INDEX_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_LAST_BISULFITE_INDEX_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}}/)
    TOPLEVEL_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_LAST_BISULFITE_INDEX_TOPLEVEL_IDS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}}/)
    REFERENCE=
    #make sure index exists
    if [ ! -e ${INDEX_DIRECTORY}/index_f.prj ]; then
        failThisLine
      	dieUponError "last-bisulfite index could not be found"
    fi

    # Avoid mapping bias by converting cytosines to thymines (use small t so we can convert back to Cs after mapping).
    perl -pe 'y/C/t/ if $. % 4 == 2' ${READS_FASTQ} >READS-${TEMP_FILENAME1}.fastq

    ${RESOURCES_LAST_EXEC_PATH} -v -p ${RESOURCES_LAST_BISULFITE_FORWARD_MATRIX} -s1 -Q1 -d${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_D} \
        -e${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_E} ${INDEX_DIRECTORY}/index_f READS-${TEMP_FILENAME1}.fastq -o f-${TEMP_FILENAME1}.maf

    ${RESOURCES_LAST_EXEC_PATH} -v -p ${RESOURCES_LAST_BISULFITE_REVERSE_MATRIX} -s0 -Q1 -d${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_D} \
              -e${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_E} ${INDEX_DIRECTORY}/index_r READS-${TEMP_FILENAME1}.fastq -o r-${TEMP_FILENAME1}.maf

    ${RESOURCES_LAST_MERGE_BATCHES_EXEC} f-${TEMP_FILENAME1}.maf r-${TEMP_FILENAME1}.maf | ${RESOURCES_LAST_MAP_PROBS_EXEC} -s${PLUGINS_ALIGNER_LAST_BISULFITE_V2_ARTIFACT_S} > ${TEMP_FILENAME1}-both-strands.maf
    REFERENCE=${TOPLEVEL_DIRECTORY}/toplevel-ids.compact-reads

    java -Xmx${PLUGIN_NEED_ALIGN_JVM} -Dlog4j.debug=true -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                        -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                        -jar ${RESOURCES_GOBY_GOBY_JAR} \
                        --mode last-to-compact -i ${TEMP_FILENAME1}-both-strands.maf -o ${OUTPUT} --third-party-input true \
                        --only-maf -q ${FULL_READS_INPUT} -t ${REFERENCE} --quality-filter-parameters threshold=1.0 \
                        --substitutions t/C,a/G
