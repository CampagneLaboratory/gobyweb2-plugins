
FULL_READS_INPUT=$1
PAIRED_END_ALIGNMENT=$2
READS_FASTQ=$3
TEMP_FILENAME1=$4
OUTPUT=$5
JOB_DIR=$6

# Grab the variables and functions we need:
. ${JOB_DIR}/constants.sh
. ${JOB_DIR}/auto-options.sh

set -x
  # Avoid mapping bias by converting cytosines to thymines (use small t so we can convert back to Cs after mapping).
  perl -pe 'y/C/t/ if $. % 4 == 2' ${READS_FASTQ} >READS-${TEMP_FILENAME1}.fastq

  ${RESOURCES_LAST_EXEC_PATH} -v -p ${RESOURCES_LAST_BISULFITE_FORWARD_MATRIX} -s1 -Q1 -d${PLUGINS_ALIGNER_LAST_BISULFITE_V2_D} \
        -e${PLUGINS_ALIGNER_LAST_BISULFITE_V2_E} ${INDEX_DIRECTORY}/index_f READS-${TEMP_FILENAME1}.fastq -o f-${TEMP_FILENAME1}.maf

  ${RESOURCES_LAST_EXEC_PATH} -v -p ${RESOURCES_LAST_BISULFITE_REVERSE_MATRIX} -s0 -Q1 -d${PLUGINS_ALIGNER_LAST_BISULFITE_V2_D} \
              -e${PLUGINS_ALIGNER_LAST_BISULFITE_V2_E} ${INDEX_DIRECTORY}/index_r READS-${TEMP_FILENAME1}.fastq -o r-${TEMP_FILENAME1}.maf

  ${RESOURCES_LAST_MERGE_BATCHES_EXEC} f-${TEMP_FILENAME1}.maf r-${TEMP_FILENAME1}.maf | ${RESOURCES_LAST_MAP_PROBS_EXEC} -s${PLUGINS_ALIGNER_LAST_BISULFITE_V2_S} > ${TEMP_FILENAME1}-both-strands.maf

  java -Xmx${PLUGIN_NEED_ALIGN_JVM} -Dlog4j.debug=true -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                        -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                        -jar ${RESOURCES_GOBY_GOBY_JAR} \
                        --mode last-to-compact -i ${TEMP_FILENAME1}-both-strands.maf -o ${OUTPUT} --third-party-input true \
                        --only-maf -q ${FULL_READS_INPUT} -t ${REFERENCE} --quality-filter-parameters threshold=1.0 \
                        --substitutions t/C,a/G
