# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:


# GROUPS_DEFINITION = description of the groups, in the format group-1=sample_i,sample_j/group-2=sample_k,..
# COMPARE_DEFINITION

# All output files must be created in the directory where the analysis script is run.
# the script generates one TSV file with the statistics, as well as images for the scatter plots:
# GENE.png
# EXON.png
# OTHER.png
# TRANSCRIPT.png

#set -o errexit
#set -o nounset

function evaluate {
EVAL=raw-counts
}

function setupWeights {

   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_WEIGHT_ADJUSTMENT}" == "NONE" ]; then

       USE_WEIGHTS_DIRECTIVE=" "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_WEIGHT_ADJUSTMENT}" == "GC_CONTENT" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights gc --adjust-gc-bias ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_BIAS_ADJUSTMENT_FORMULA} "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_WEIGHT_ADJUSTMENT}" == "HEPTAMERS" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights heptamers "
   else
     dieUponError "weight adjustment  not supported: ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_WEIGHT_ADJUSTMENT}"
   fi

}

function setupAnnotationTypes {
   ANNOTATION_TYPES=""
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_ESTIMATE_COUNTS_GENE}" == "true" ]; then

       ANNOTATION_TYPES="${ANNOTATION_TYPES}gene"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_ESTIMATE_COUNTS_EXON}" == "true" ]; then
       if [ "${ANNOTATION_TYPES}" != "" ]; then
          ANNOTATION_TYPES="${ANNOTATION_TYPES},"
       fi
       ANNOTATION_TYPES="${ANNOTATION_TYPES}exon"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_ESTIMATE_COUNTS_OTHER}" == "true" ]; then
       if [ "${ANNOTATION_TYPES}" != "" ]; then
          ANNOTATION_TYPES="${ANNOTATION_TYPES},"
       fi
       ANNOTATION_TYPES="${ANNOTATION_TYPES}other"
   fi

   if [ "${ANNOTATION_TYPES}" == "" ]; then
     dieUponError "At least one annotation type must be selected to run a differential analysis."
   fi

}

function setupAnnotationSource {
    . ${JOB_DIR}/artifacts.sh
    expose_artifact_environment_variables
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
    ANNOTATION_PATH=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_ANNOTATIONS_ANNOTATIONS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

  if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_ANNOTATION_SOURCE}" == "GENE_EXON_OTHER" ]; then
    # gene exon annotation file.
    ANNOTATION_SOURCE="${ANNOTATION_PATH}/exon-annotations.tsv"
  else
    # CNV annotation file.
    ANNOTATION_SOURCE="${ANNOTATION_PATH}/cnv-annotations.tsv"
  fi
}

. ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_FILES_PARALLEL_SCRIPT}

function plugin_alignment_analysis_combine {
   set -x
   set -T
   RESULT_FILE=$1
   shift
   PART_RESULT_FILES=$*

   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_NUM_TOP_HITS}
   Q_VALUE_THRESHOLD=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_Q_VALUE_THRESHOLD}

   run_fdr

   # Estimate stats on complete file
   NORMALIZATION_METHOD="${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_NORMALIZATION_METHOD}"
   if [ -z "${NORMALIZATION_METHOD}" ]; then
        NORMALIZATION_METHOD="aligned-count"
   fi

   setupWeights

   run_goby ${PLUGIN_NEED_COMBINE_JVM}  stats --info info.xml \
          ${OUT_FILENAME} \
          --parallel \
          --groups ${GROUPS_DEFINITION} \
          --compare ${COMPARE_DEFINITION} ${USE_WEIGHTS_DIRECTIVE} \
          --normalization-methods ${NORMALIZATION_METHOD} \
          -o stats.tsv

   dieUponError "statistics evaluation failed."

   if [ $RETURN_STATUS -eq 0 ]; then
            IMAGE_OUTPUT_PNG=
            R -f ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ARTIFACT_FILES_R_SCRIPT} --slave --quiet --no-restore --no-save --no-readline --args input=stats.tsv graphOutput=.png
   fi

}
