# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# IS_TRANSCRIPT = whether alignments were done against a cDNA reference.
# GROUPS_DEFINITION = description of the groups, in the format group-1=sample_i,sample_j/group-2=sample_k,..
# COMPARE_DEFINITION
# ANNOTATION_FILE = file describing annotations in the Goby annotation format.
# ANNOTATION_TYPES = gene|exon|other, specifies the kind of annotations to calculate counts for.
# USE_WEIGHTS_DIRECTIVE = optional, command line flags to have Goby annotation-to-counts adjust counts with weights.

# All output files must be created in the directory where the analysis script is run.
# the script generates one TSV file with the statistics, as well as images for the scatter plots:
# GENE.png
# EXON.png
# OTHER.png
# TRANSCRIPT.png

function eval {
EVAL=raw-counts
}

function setupWeights {

   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_WEIGHT_ADJUSTMENT}" == "NONE" ]; then

       USE_WEIGHTS_DIRECTIVE=" "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_WEIGHT_ADJUSTMENT}" == "GC_CONTENT" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights gc --adjust-gc-bias ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_BIAS_ADJUSTMENT_FORMULA} "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_WEIGHT_ADJUSTMENT}" == "HEPTAMERS" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights heptamers "
   else
     dieUponError "weight adjustment  not supported: ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_WEIGHT_ADJUSTMENT}"
   fi

}

function setupAnnotationTypes {
   ANNOTATION_TYPES=""
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ESTIMATE_COUNTS_GENE}" == "true" ]; then

       ANNOTATION_TYPES="${ANNOTATION_TYPES}gene"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ESTIMATE_COUNTS_EXON}" == "true" ]; then
       if [ "${ANNOTATION_TYPES}" != "" ]; then
          ANNOTATION_TYPES="${ANNOTATION_TYPES},"
       fi
       ANNOTATION_TYPES="${ANNOTATION_TYPES}exon"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ESTIMATE_COUNTS_OTHER}" == "true" ]; then
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

  ANNOTATION_SOURCE=""
  if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_ANNOTATION_SOURCE}" == "GENE_EXON_OTHER" ]; then
    # gene exon annotation file.
    ANNOTATION_SOURCE="${REFERENCE_DIRECTORY}/exon-annotations.tsv"
  else
    # CNV annotation file.
    ANNOTATION_SOURCE="${REFERENCE_DIRECTORY}/cnv-annotations.tsv"
  fi
}

. ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_FILES_PARALLEL_SCRIPT}

function plugin_alignment_analysis_combine {
   set -x
   set -T
   RESULT_FILE=$1
   shift
   PART_RESULT_FILES=$*

   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_NUM_TOP_HITS}
   Q_VALUE_THRESHOLD=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_Q_VALUE_THRESHOLD}

   run_fdr

   # Estimate stats on complete file
   NORMALIZATION_METHOD="${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_NORMALIZATION_METHOD}"
   if [ -z "${NORMALIZATION_METHOD}" ]; then
        NORMALIZATION_METHOD="aligned-count"
   fi

   setupWeights

   run-goby ${PLUGIN_NEED_COMBINE_JVM}  stats --info info.xml \
          ${OUT_FILENAME} \
          --parallel \
          --groups ${GROUPS_DEFINITION} \
          --compare ${COMPARE_DEFINITION} ${USE_WEIGHTS_DIRECTIVE} \
          --normalization-methods ${NORMALIZATION_METHOD} \
          -o stats.tsv

   dieUponError "statistics evaluation failed."

   if [ $RETURN_STATUS -eq 0 ]; then
            IMAGE_OUTPUT_PNG=
            R -f ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_GOBY_FILES_R_SCRIPT} --slave --quiet --no-restore --no-save --no-readline --args input=stats.tsv graphOutput=.png
   fi

}
