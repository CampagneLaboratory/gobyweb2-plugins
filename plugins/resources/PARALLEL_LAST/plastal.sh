[ $# -eq 6 ] && shift
JOB_DIR=$1
ALIGN_COMMAND=$2
READS_FOR_SPLIT=$3
OUTPUT=$4
PAIRED_END_ALIGNMENT=$5

NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
NUM_THREADS=$((${NUM_THREADS} - 2))
if [ ${NUM_THREADS} -lt 2 ]; then
   # Make sure we use at least two threads:
   NUM_THREADS=2
fi

# Grab the variables and functions we need:
. ${JOB_DIR}/constants.sh
. ${JOB_DIR}/auto-options.sh
. ${RESOURCES_GOBY_SHELL_SCRIPT}
if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
  PAIR="%pair.fastq%"
  PAIRED_ARG="--paired"
else
  PAIR=" "
  PAIRED_ARG=" "
  PAIRED_END_ALIGNMENT="false"
fi

# This script splits a compact-reads file in num-parts (n) and executes ALIGN_COMMAND on each bit. It concatenates alignments that result to produce OUTPUT
java ${GRID_JVM_FLAGS} -Dlog4j.debug=true -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                       -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_GOBY_GOBY_JAR} \
                       --mode run-parallel -i "${READS_FOR_SPLIT}" -n "${NUM_THREADS}" ${PAIRED_ARG} -o "${OUTPUT}"  ${ALIGN_COMMAND} ${READS_FOR_SPLIT} ${PAIRED_END_ALIGNMENT} %read.fastq% ${PAIR} %tmp1% %output% ${JOB_DIR}
