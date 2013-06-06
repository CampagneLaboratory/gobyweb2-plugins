# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# IS_TRANSCRIPT = whether alignments were done against a cDNA reference.
# GROUPS_DEFINITION = description of the groups, in the format group-1=sample_i,sample_j/group-2=sample_k,..
# COMPARE_DEFINITION
# ANNOTATION_FILE = file describing annotations in the Goby annotation format.
# ANNOTATION_SOURCE = gene|exon|other, specifies the kind of annotations to calculate counts for.
# USE_WEIGHTS_DIRECTIVE = optional, command line flags to have Goby annotation-to-counts adjust counts with weigths.

# All output files must be created in the directory where the analysis script is run.
# STATS_OUTPUT = name of the statistics file produced by the analysis. Format can be tsv, or VCF. If the file is VCF,
# the filename points to the vcf.gz file, and a secondary index file vcf.gz.tbi must also be produced by the analysis.
# IMAGE_OUTPUT_PNG = name of an optional image file output (must be written in PNG format)

# OTHER_ALIGNMENT_ANALYSIS_OPTIONS = any options defined by the end-user or assembled with the auto-format mechanism.

. ${RESOURCES_GOBY_SHELL_SCRIPT}
. ${RESOURCES_IGVTOOLS_SHELL_SCRIPT}

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
   else

    echo 0
   fi
}


function setupAnnotationSource {
    . ${JOB_DIR}/artifacts.sh
    expose_artifact_environment_variables
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
    ANNOTATION_PATH=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_ANNOTATIONS_ANNOTATIONS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
    ANNOTATION_CLASS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_ANNOTATION_TYPE}

    local ANNOTATION_SOURCE=""
    case "${ANNOTATION_CLASS}" in
                  "ENSEMBL_PROMOTER")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/promoter-annotations.tsv"
                  ;;
                  "GENES")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/gene-annotations.tsv"
                        ;;
                  "EXONS")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/exon-annotations.tsv"
                  ;;
                  "INTRONS")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/intron-annotations.tsv"
                  ;;
                  "5_UTR")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/5_utrs-annotations.tsv"
                  ;;
                  "3_UTR")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/3_utr-annotations.tsv"
                  ;;
                  "CPG_ISLAND")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/promoter-annotations.tsv"
                  ;;
                  "INTERGENIC")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/intergenic-annotations.tsv"
                  ;;
                  "1KB_TILE")
                         ANNOTATION_SOURCE="${ANNOTATION_PATH}/tile-annotations.tsv"
                  ;;
                  *)
                          echo "ERROR: Invalid type of annotation class: ${ANNOTATION_CLASS}"
    esac

    echo ${ANNOTATION_SOURCE}
}

function setupGenomeCache{
    set -x
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    SEQUENCE_CACHE_DIR=$(eval echo \${RESOURCES_ARTIFACTS_GOBY_INDEXED_GENOMES_SEQUENCE_CACHE_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

    echo ${SEQUENCE_CACHE_DIR}
}
function plugin_alignment_analysis_process {
   SLICING_PLAN_FILENAME=$1
   ARRAY_JOB_INDEX=$2
   shift
   shift
   MINIMUM_VARIATION_SUPPORT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_MINIMUM_VARIATION_SUPPORT}
   THRESHOLD_DISTINCT_READ_INDICES=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_THRESHOLD_DISTINCT_READ_INDICES}
   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_OUTPUT_FORMAT}

    ANNOTATION_SOURCE=$(setupAnnotationSource)

   if [ ! "${ANNOTATION_TYPE}" == "NONE" ]; then
     ANNOTATION_OPTION=" -x MethylationRegionsOutputFormat:annotations=${ANNOTATION_SOURCE} "
   else
     ANNOTATION_OPTION=" "
   fi

   # These variables are defined: SLICING_PLAN_FILENAME
     echo "Processing run_single_alignment_analysis_process for part ${SGE_TASK_ID}"

     WINDOW_LIMITS=`awk -v arrayJobIndex=${ARRAY_JOB_INDEX} '{ if (lineNumber==arrayJobIndex) print " -s "$3" -e "$6; lineNumber++; }' ${SLICING_PLAN_FILENAME}`
     STAT2_FILENAME=${SGE_O_WORKDIR}/results/${TAG}-variations-stats2.tsv

     echo "Discovering sequence variants for window limits: ${WINDOW_LIMITS} and statsFilename: ${STAT2_FILENAME}"

     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Start discover-sequence-variations for part # ${CURRENT_PART}." --index ${CURRENT_PART} --job-type job-part
     REALIGNMENT_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_REALIGN_AROUND_INDELS}
     if [ "${REALIGNMENT_OPTION}" == "true" ]; then

            REALIGNMENT_ARGS=" --processor realign_near_indels "
     else
            REALIGNMENT_ARGS="  "
     fi
     CALL_INDELS_OPTION=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_CALL_INDELS}
     INDEL_RATE=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_INDEL_RATE}
     FORCE_DIPLOID=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_FORCE_DIPLOID}
     WRITE_COUNTS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_WRITE_COUNTS}
     ESTIMATE_DENSITY=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_ESTIMATE_INTRA_GROUP_DIFFERENCE_DENSITY}
     COMBINATOR=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_PVALUE_COMBINATOR}
     EXTRA_ARGS=" "
     CONTEXTS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_CONTEXTS}

     if [ "${INDEL_RATE}" == "true" ]; then
         CALL_INDELS_OPTION="true"
         EXTRA_ARGS=" -x MethylationRegionsOutputFormat:do-indel-rate=true "
     fi


     if [ "${ESTIMATE_DENSITY}" == "true" ]; then
         run_methyl_regions ${TAG}-intra-group-differences-estimate-${ARRAY_JOB_INDEX}.bin -x AnnotationAveragingWriter:estimate-intra-group-differences=${ESTIMATE_DENSITY} -x AnnotationAveragingWriter:estimate-empirical-P=false -x AnnotationAveragingWriter:binning-strategy=fastslog10

         dieUponError  "Estimating density failed for part ${CURRENT_PART}."
         mkdir -p ${SGE_O_WORKDIR}/split-results/
         cp ${TAG}-intra-group-differences-estimate-${ARRAY_JOB_INDEX}.bin ${SGE_O_WORKDIR}/split-results/
         #cp *-null-observations.tsv ${SGE_O_WORKDIR}/split-results/

         dieUponError  "Could not copy estimated density to result directory for part ${CURRENT_PART}."
         EXTRA_ARGS=" -x AnnotationAveragingWriter:estimate-empirical-P=true -x AnnotationAveragingWriter:estimate-intra-group-differences=false -x AnnotationAveragingWriter:serialized-estimator-filename=${TAG}-intra-group-differences-estimate-${ARRAY_JOB_INDEX}.bin -x AnnotationAveragingWriter:combinator=${COMBINATOR} "
     fi

     run_methyl_regions ${TAG}-mr-${ARRAY_JOB_INDEX}.tsv
     dieUponError  "Compare methylation region part, sub-task ${CURRENT_PART} failed."
     #cp *-test-observations.tsv ${SGE_O_WORKDIR}/split-results/

     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "End discover-sequence-variations for part # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part


}
function run_methyl_regions {
    output="$1"
    shift

    SEQUENCE_CACHE_DIR=$(setupGenomeCache)

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
           --output ${output}  \
           --call-indels ${CALL_INDELS_OPTION} \
           ${ANNOTATION_OPTION} \
           --diploid ${FORCE_DIPLOID} \
           -x AnnotationAveragingWriter:contexts:${CONTEXTS} \
           -x AnnotationAveragingWriter:write-counts=${WRITE_COUNTS} \
           ${EXTRA_ARGS} \
           ${ENTRIES_FILES} $*
}

function plugin_alignment_analysis_combine {

   RESULT_FILE=stats.tsv
   shift
   PART_RESULT_FILES=$*

   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_OUTPUT_FORMAT}
   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_NUM_TOP_HITS}
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
     # only adjust fisher p values, because empirical p-values are already controlling FDR.
        COLUMNS="${COLUMNS} --column-selection-filter fisherP[${GROUP_PAIR}]"
    done

   echo "Adjusting P-value columns: $COLUMNS"
   
   Q_VALUE_THRESHOLD=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_METHYLATION_REGIONS_ARTIFACT_Q_VALUE_THRESHOLD}
   ADJUSTMENT_OPTIONS=" --q-threshold ${Q_VALUE_THRESHOLD}  --top-hits ${NUM_TOP_HITS} "

   if [ "${NUM_GROUPS}" == "1" ]; then

       ADJUSTMENT_OPTIONS=" "
   fi


   run-goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
          ${ADJUSTMENT_OPTIONS} \
          ${PART_RESULT_FILES}  \
          ${COLUMNS} \
          --output ${RESULT_FILE}
   cp ${RESULT_FILE} raw.igv

   igvtools ${PLUGIN_NEED_COMBINE_JVM} sort raw.igv methyl-regions.igv
   # We can't convert to tdf at this time because we don't know what genome build id to use in the following command line
   # igvtools ${PLUGIN_NEED_COMBINE_JVM} tile methyl-regions.igv methyl-regions.tdf hg19
}