#!/bin/sh

#
. ${JOB_DIR}/constants.sh

READ_FILES_LIST=""

function plugin_task {
     set -xv
     cd ${TMPDIR}
     echo "fileset command: ${FILESET_COMMAND}"

     ${FILESET_COMMAND} --has-fileset INPUT_VCF
     dieUponError "Input files are not available"

     VCF_INPUT=`${FILESET_COMMAND} --fetch INPUT_VCF`
     dieUponError "Failed to fetch input vcf files ${VCF_INPUT}"
     echo "${VCF_INPUT}"


     mv ${VCF_INPUT} ./
     # vcf-annotate requires uncompressed vcf:
     gunzip *.vcf.gz

     export JAVA_OPTS="${GRID_JVM_FLAGS}"

     . ${RESOURCES_ANNOTATE_VCF_EXEC_PATH}
     annotate_vep true *.vcf out-vep.vcf ${PLUGINS_TASK_ANNOTATE_WITH_VEP_ONLY_NON_SYNONYMOUS} 2>&1  | tee ${JOB_DIR}/log.txt

     STATUS=$?

     # push back the output stats:
     LOG_REGISTERED_TAGS=`${FILESET_COMMAND} --push EXECUTION_LOG: ${JOB_DIR}/log.txt `
     dieUponError "Failed to push back the execution log"

     VCF_REGISTERED_TAGS=`${FILESET_COMMAND} --push ANNOTATED_VCF: out-vep.vcf `
     dieUponError "Failed to push back the annotated VCF"

     ALL_REGISTERED_TAGS="${LOG_REGISTERED_TAGS} ${VCF_REGISTERED_TAGS}"
     echo "The following tags were registered by this plugin: ${ALL_REGISTERED_TAGS}"
     exit $STATUS
}


