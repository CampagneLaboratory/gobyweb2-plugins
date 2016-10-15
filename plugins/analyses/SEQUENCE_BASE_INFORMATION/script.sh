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

. ${RESOURCES_GOBY3_SHELL_SCRIPT}

function plugin_alignment_analysis_split {

  NUMBER_OF_PARTS=$1
  SPLICING_PLAN_RESULT=$2
  shift
  shift
  goby suggest-position-slices \
          --number-of-bytes 5000000 \
          --output ${SPLICING_PLAN_RESULT} \
          $*
}

# This function return the number of parts in the slicing plan. It returns zero if the alignments could not be split.
function plugin_alignment_analysis_num_parts {
   SPLICING_PLAN_FILE=$1

   if [ $? -eq 0 ]; then

        echo `grep -v targetIdStart ${SPLICING_PLAN_FILE} | wc -l `

   else
        echo 0
   fi


}

function plugin_alignment_analysis_process {

   SLICING_PLAN_FILENAME=$1
   ARRAY_JOB_INDEX=$2
   shift
   shift
   MINIMUM_VARIATION_SUPPORT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_MINIMUM_VARIATION_SUPPORT}
   THRESHOLD_DISTINCT_READ_INDICES=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_THRESHOLD_DISTINCT_READ_INDICES}
   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_OUTPUT_FORMAT}

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
     REALIGNMENT_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_REALIGN_AROUND_INDELS}
     if [ "${REALIGNMENT_OPTION}" == "true" ]; then

            REALIGNMENT_ARGS=" --processor realign_near_indels "
     else
            REALIGNMENT_ARGS="  "
     fi
     CALL_INDELS_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_CALL_INDELS}
     FORCE_DIPLOID=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_FORCE_DIPLOID}

     # Note that we override the grid jvm flags to request only 4Gb:
     run_goby ${PLUGIN_NEED_PROCESS_JVM} discover-sequence-variants \
           ${WINDOW_LIMITS} \
           --groups ${GROUPS_DEFINITION} \
           --compare ${COMPARE_DEFINITION} \
           --format ${OUTPUT_FORMAT} \
           --eval filter \
           ${COVARIATES_OPTION} \
           ${REALIGNMENT_ARGS} \
           --genome ${SEQUENCE_CACHE_DIR}/random-access-genome \
           --minimum-variation-support ${MINIMUM_VARIATION_SUPPORT} \
           --threshold-distinct-read-indices ${THRESHOLD_DISTINCT_READ_INDICES} \
           --output ${TAG}-out-${CURRENT_PART}  \
           --call-indels ${CALL_INDELS_OPTION} \
           --diploid ${FORCE_DIPLOID}  \
           ${ENTRIES_FILES}

      dieUponError  "Compare sequence variations part, sub-task ${CURRENT_PART} failed."
      mkdir -p ${JOB_DIR}/split-results
      dieUponError  "cannot create split-results directory. sub-task ${CURRENT_PART} failed."
      cp ${TAG}-out-${CURRENT_PART}.sbi  ${JOB_DIR}/split-results/
      cp ${TAG}-out-${CURRENT_PART}.sbip  ${JOB_DIR}/split-results/
      ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "End discover-sequence-variations for part # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part
      # Create an empty TSV file
      touch ${TAG}-out-${ARRAY_JOB_INDEX}.tsv

}

function plugin_alignment_analysis_combine {

   RESULT_FILE=stats.vcf.gz
   shift
   PART_RESULT_FILES=$*

   java -cp ${RESOURCES_DLVARIATION_JAR}/model-training-bin.jar  -Xmx${PLUGIN_NEED_COMBINE_JVM}  \
                                        org.campagnelab.dl.varanalysis.intermediaries.Randomizer2 \
                                                    ${JOB_DIR}/split-results/*.sbi ${TMPDIR}/${TAG}-out

   mkdir -p ${JOB_DIR}/results
   cp ${TMPDIR}/${TAG}-out.sbi* ${JOB_DIR}/results
    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Result written to JOB_DIR/results (will not be pushed to web server) # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part
   echo "Result written to JOB_DIR/results (will not be pushed to web server)">  ${TMPDIR}/${TAG}-out.tsv
}