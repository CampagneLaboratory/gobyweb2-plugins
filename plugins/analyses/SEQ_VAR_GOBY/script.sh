# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# IS_TRANSCRIPT = whether alignments were done against a cDNA reference.
# GROUPS_DEFINITION = description of the groups, in the format group-1=sample_i,sample_j/group-2=sample_k,..
# COMPARE_DEFINITION
# ANNOTATION_FILE = file describing annotations in the Goby annotation format.
# ANNOTATION_TYPES = gene|exon|other, specifies the kind of annotations to calculate counts for.
# USE_WEIGHTS_DIRECTIVE = optional, command line flags to have Goby annotation-to-counts adjust counts with weigths.

# All output files must be created in the directory where the analysis script is run.
# STATS_OUTPUT = name of the statistics file produced by the analysis. Format can be tsv, or VCF. If the file is VCF,
# the filename points to the vcf.gz file, and a secondary index file vcf.gz.tbi must also be produced by the analysis.
# IMAGE_OUTPUT_PNG = name of an optional image file output (must be written in PNG format)

# OTHER_ALIGNMENT_ANALYSIS_OPTIONS = any options defined by the end-user or assembled with the auto-format mechanism.

. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_alignment_analysis_split {

  NUMBER_OF_PARTS=$1
  SPLICING_PLAN_RESULT=$2
  shift
  shift
  goby suggest-position-slices \
          --number-of-slices ${NUMBER_OF_PARTS} \
          --output ${SPLICING_PLAN_RESULT} \
          $*
}

# This function return the number of parts in the slicing plan. It returns zero if the alignments could not be split.
function plugin_alignment_analysis_num_parts {
   SPLICING_PLAN_FILE=$1

   if [ $? -eq 0 ]; then

        return `grep -v targetIdStart ${SPLICING_PLAN_FILE} | wc -l `
   fi

   return 0
}

function plugin_alignment_analysis_process {

   SLICING_PLAN_FILENAME=$1
   ARRAY_JOB_INDEX=$2
   shift
   shift
   MINIMUM_VARIATION_SUPPORT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_MINIMUM_VARIATION_SUPPORT}
   THRESHOLD_DISTINCT_READ_INDICES=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_THRESHOLD_DISTINCT_READ_INDICES}
   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_OUTPUT_FORMAT}
   ANNOTATIONS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_ANNOTATIONS}
   if [ ! "${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_ANNOTATIONS}" == "NONE" ]; then
     ANNOTATION_OPTION=" -x MethylationRegionsOutputFormat:annotations=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_ANNOTATIONS} "
   else
     ANNOTATION_OPTION=" "
   fi
    set -x
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    SEQUENCE_CACHE_DIR=$(eval echo \${RESOURCES_ARTIFACTS_GOBY_INDEXED_GENOMES_SEQUENCE_CACHE_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})


   # These variables are defined: SLICING_PLAN_FILENAME
     echo "Processing run_single_alignment_analysis_process for part ${SGE_TASK_ID}"

     WINDOW_LIMITS=`awk -v arrayJobIndex=${ARRAY_JOB_INDEX} '{ if (lineNumber==arrayJobIndex) print " -s "$3" -e "$6; lineNumber++; }' ${SLICING_PLAN_FILENAME}`
     STAT2_FILENAME=${SGE_O_WORKDIR}/results/${TAG}-variations-stats2.tsv

     echo "Discovering sequence variants for window limits: ${WINDOW_LIMITS} and statsFilename: ${STAT2_FILENAME}"

     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Start discover-sequence-variations for part # ${CURRENT_PART}." --index ${CURRENT_PART} --job-type job-part
     REALIGNMENT_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_REALIGN_AROUND_INDELS}
     if [ "${REALIGNMENT_OPTION}" == "true" ]; then

            REALIGNMENT_ARGS=" --processor realign_near_indels "
     else
            REALIGNMENT_ARGS="  "
     fi
     CALL_INDELS_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_CALL_INDELS}
     FORCE_DIPLOID=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_FORCE_DIPLOID}

     # Note that we override the grid jvm flags to request only 4Gb:
     run-goby ${PLUGIN_NEED_PROCESS_JVM} discover-sequence-variants \
           ${WINDOW_LIMITS} \
           --groups ${GROUPS_DEFINITION} \
           --compare ${COMPARE_DEFINITION} \
           --format ${OUTPUT_FORMAT} \
           --eval filter \
           ${REALIGNMENT_ARGS} \
           --genome ${SEQUENCE_CACHE_DIR}/random-access-genome \
           --minimum-variation-support ${MINIMUM_VARIATION_SUPPORT} \
           --threshold-distinct-read-indices ${THRESHOLD_DISTINCT_READ_INDICES} \
           --output ${TAG}-dsv-${ARRAY_JOB_INDEX}.vcf  \
           --call-indels ${CALL_INDELS_OPTION} \
           --diploid ${FORCE_DIPLOID}  \
           ${ANNOTATION_OPTION} \
           ${ENTRIES_FILES}

      dieUponError  "Compare sequence variations part, sub-task ${CURRENT_PART} failed."

      ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "End discover-sequence-variations for part # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part

     . ${RESOURCES_ANNOTATE_VCF_EXEC_PATH}
     annotate_vep ${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_ANNOTATE_VARIATIONS} \
            ${TAG}-dsv-${ARRAY_JOB_INDEX}.vcf ${TAG}-dsv-${ARRAY_JOB_INDEX}-vep.vcf

     annotate_ensembl_genes ${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_ANNOTATE_VARIATIONS} \
                 ${TAG}-dsv-${ARRAY_JOB_INDEX}-vep.vcf ${TAG}-discover-sequence-variants-output-${ARRAY_JOB_INDEX}.vcf.gz


}

function plugin_alignment_analysis_combine {

   RESULT_FILE=stats.vcf.gz
   shift
   PART_RESULT_FILES=$*

   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_OUTPUT_FORMAT}
   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_NUM_TOP_HITS}
   COLUMNS=" "
   for groupName in {1..${NUM_GROUPS}}
   do
      echo -n "${groupName} "
   done
    #GROUP1_PAIR=A/B
    #GROUP2_NAME=A/C
    #NUM_COMPARISON_PAIRS=2

    for ((i=1; i <= NUM_COMPARISON_PAIRS ; i++))
    do
     GROUP_PAIR=`eval echo "$""GROUP"$i"_COMPARISON_PAIR"`

     # More than one group, some P-values may need adjusting:
     if [ "${OUTPUT_FORMAT}" == "ALLELE_FREQUENCIES" ]; then

        COLUMNS="${COLUMNS} --column P[${GROUP_PAIR}]"

       elif [ "${OUTPUT_FORMAT}" == "COMPARE_GROUPS" ]; then

         COLUMNS="${COLUMNS} --column FisherP[${GROUP_PAIR}]"

       elif [ "${OUTPUT_FORMAT}" == "METHYLATION" ]; then

         COLUMNS="${COLUMNS} --column FisherP[${GROUP_PAIR}]"
       fi
    done

   echo "Adjusting P-value columns: $COLUMNS"
   if [ "${OUTPUT_FORMAT}" == "GENOTYPES" -o ${NUM_GROUPS} == 1 ]; then

        # Do not attempt FDR adjustment when there is no p-value, just concat the split files and sort:

        ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-concat ${PART_RESULT_FILES} | \
        ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/bin/vcf-sort | \
        ${RESOURCES_TABIX_BGZIP_EXEC_PATH} -c > ${RESULT_FILE}

   else
       Q_VALUE_THRESHOLD=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_Q_VALUE_THRESHOLD}

        run-goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
          --vcf \
          --q-threshold ${Q_VALUE_THRESHOLD} \
          --top-hits ${NUM_TOP_HITS} \
          ${PART_RESULT_FILES}  \
          ${COLUMNS} \
          --output ${TMPDIR}/${TAG}-pre.vcf.gz

        gunzip -c -d ${TMPDIR}/${TAG}-pre.vcf.gz | ${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/vcf-sort | ${RESOURCES_TABIX_BGZIP_EXEC_PATH} -c > ${RESULT_FILE}
   fi

   ${TABIX_EXEC_PATH} -f -p vcf ${RESULT_FILE}

}