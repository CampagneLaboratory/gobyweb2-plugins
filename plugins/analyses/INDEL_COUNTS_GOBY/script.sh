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
          --number-of-bytes 50000000 \
          --output ${SPLICING_PLAN_RESULT} \
          $*
}


# This function return the number of parts in the slicing plan. It returns zero if the alignments could not be split.
function plugin_alignment_analysis_num_parts {
   SPLICING_PLAN_FILE=$1

   if [ $? -eq 0 ]; then

        echo `grep -v targetIdStart ${SPLICING_PLAN_FILE} | wc -l `
   fi

   echo 0
}

function plugin_alignment_analysis_process {
   SLICING_PLAN_FILENAME=$1
   ARRAY_JOB_INDEX=$2
   shift
   shift
   MINIMUM_VARIATION_SUPPORT=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_MINIMUM_VARIATION_SUPPORT}
   THRESHOLD_DISTINCT_READ_INDICES=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_THRESHOLD_DISTINCT_READ_INDICES}
   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_OUTPUT_FORMAT}

   # These variables are defined: SLICING_PLAN_FILENAME
     echo "Processing run_single_alignment_analysis_process for part ${SGE_TASK_ID}"

     WINDOW_LIMITS=`awk -v arrayJobIndex=${ARRAY_JOB_INDEX} '{ if (lineNumber==arrayJobIndex) print " -s "$3" -e "$6; lineNumber++; }' ${SLICING_PLAN_FILENAME}`
     STAT2_FILENAME=${SGE_O_WORKDIR}/results/${TAG}-variations-stats2.tsv

     echo "Discovering sequence variants for window limits: ${WINDOW_LIMITS} and statsFilename: ${STAT2_FILENAME}"

     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Start discover-sequence-variations for part # ${CURRENT_PART}." --index ${CURRENT_PART} --job-type job-part
     REALIGNMENT_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_REALIGN_AROUND_INDELS}
     if [ "${REALIGNMENT_OPTION}" == "true" ]; then

            REALIGNMENT_ARGS=" --processor realign_near_indels "
     else
            REALIGNMENT_ARGS="  "
     fi

     # Note that we override the grid jvm flags to request only 4Gb:
     run-goby ${PLUGIN_NEED_PROCESS_JVM} discover-sequence-variants \
           ${WINDOW_LIMITS} \
           --groups ${GROUPS_DEFINITION} \
           --compare ${COMPARE_DEFINITION} \
           --format ${OUTPUT_FORMAT} \
           --eval filter \
           ${REALIGNMENT_ARGS} \
           --genome ${REFERENCE_DIRECTORY}/random-access-genome \
           --minimum-variation-support ${MINIMUM_VARIATION_SUPPORT} \
           --threshold-distinct-read-indices ${THRESHOLD_DISTINCT_READ_INDICES} \
           --output indel-counts.tsv  \
           --call-indels true \
           ${ENTRIES_FILES}

     dieUponError  "Failed to count indels in alignment files, sub-task ${CURRENT_PART} failed."
     cp indel-counts.tsv  ${TAG}-ic-${ARRAY_JOB_INDEX}.tsv
     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "End discover-sequence-variations for part # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part

}

function plugin_alignment_analysis_combine {

   RESULT_FILE=stats.tsv
   shift
   PART_RESULT_FILES=$*


   run-goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
          ${PART_RESULT_FILES}  \
          --output ${TMPDIR}/${TAG}-indel-counts.tsv
   dieUponError  "Failed to merge split results for indel counts."

   cp ${TMPDIR}/${TAG}-indel-counts.tsv ${RESULT_FILE}

}

function plugin_alignment_analysis_sequential {
   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_OUTPUT_FORMAT}
   NUM_TOP_HITS= ${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_NUM_TOP_HITS}

   ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Start counting called indels with Goby" --index ${CURRENT_PART} --job-type job-part

   MINIMUM_VARIATION_SUPPORT=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_MINIMUM_VARIATION_SUPPORT}
   THRESHOLD_DISTINCT_READ_INDICES=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_THRESHOLD_DISTINCT_READ_INDICES}
   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_OUTPUT_FORMAT}

   REALIGNMENT_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_INDEL_COUNTS_GOBY_REALIGN_AROUND_INDELS}
   if [ "${REALIGNMENT_OPTION}" == "true" ]; then

           REALIGNMENT_ARGS=" --processor realign_near_indels "
    else
           REALIGNMENT_ARGS="  "
   fi

     # Note that we override the grid jvm flags to request only 4Gb:
   run-goby ${PLUGIN_NEED_PROCESS_JVM} discover-sequence-variants \
           --groups ${GROUPS_DEFINITION} \
           --compare ${COMPARE_DEFINITION} \
           --format ${OUTPUT_FORMAT} \
           --eval filter \
           ${REALIGNMENT_ARGS} \
           --genome ${REFERENCE_DIRECTORY}/random-access-genome \
           --minimum-variation-support ${MINIMUM_VARIATION_SUPPORT} \
           --threshold-distinct-read-indices ${THRESHOLD_DISTINCT_READ_INDICES} \
           --output indel-counts.tsv  \
           ${ENTRIES_FILES}

   dieUponError  "Failed to count indels in alignment files, sub-task ${CURRENT_PART} failed."
   RESULT_FILE=stats.tsv
   cp indel-counts.tsv ${RESULT_FILE}
   ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "End discover-sequence-variations for part # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part


}
