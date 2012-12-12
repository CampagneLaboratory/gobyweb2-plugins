
FULL_READS_INPUT=$1
READS_FASTQ=$2
TEMP_FILENAME1=$3
OUTPUT=$4
JOB_DIR=$5
# Grab the variables and functions we need:
. ${JOB_DIR}/constants.sh
. ${JOB_DIR}/auto-options.sh

set -x
  ${RESOURCES_LAST_EXEC_PATH} -v -p ${RESOURCES_LAST_BISULFITE_FORWARD_MATRIX} -s1 -Q1 -d${PLUGINS_ALIGNER_LAST_BISULFITE_D} \
        -e${PLUGINS_ALIGNER_LAST_BISULFITE_E} ${INDEX_DIRECTORY}/index_f ${READS_FASTQ} -o f-${TEMP_FILENAME1}.maf

  java -Xmx${PLUGIN_NEED_ALIGN_JVM} -Dlog4j.debug=true -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                       -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_GOBY_GOBY_JAR} \
                       --mode  last-to-compact -i f-${TEMP_FILENAME1} -o ${OUTPUT} --third-party-input true --only-maf \
                       -q ${FULL_READS_INPUT} -t ${REFERENCE} --quality-filter-parameters threshold=1.0

