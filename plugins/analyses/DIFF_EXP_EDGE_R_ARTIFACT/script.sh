# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# IS_TRANSCRIPT = whether alignments were done against a cDNA reference.
# GROUPS_DEFINITION = description of the groups, in the format group-1=sample_i,sample_j/group-2=sample_k,..
# COMPARE_DEFINITION
# ANNOTATION_FILE = file describing annotations in the Goby annotation format.
# ANNOTATION_TYPES = gene|exon|other, specifies the kind of annotations to calculate counts for.
# USE_WEIGHTS_DIRECTIVE = optional, command line flags to have Goby annotation-to-counts adjust counts with weigths.

function evaluate {
EVAL=raw-counts
}

function setupWeights {

   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_WEIGHT_ADJUSTMENT}" == "NONE" ]; then

       USE_WEIGHTS_DIRECTIVE=" "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_WEIGHT_ADJUSTMENT}" == "GC_CONTENT" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights gc --adjust-gc-bias ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_BIAS_ADJUSTMENT_FORMULA} "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_WEIGHT_ADJUSTMENT}" == "HEPTAMERS" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights heptamers "
   else
     dieUponError "weight adjustment  not supported: ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_WEIGHT_ADJUSTMENT}"
   fi

}

function setupAnnotationTypes {
   ANNOTATION_TYPES=""
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_ESTIMATE_COUNTS_GENE}" == "true" ]; then

       ANNOTATION_TYPES="${ANNOTATION_TYPES}gene"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_ESTIMATE_COUNTS_EXON}" == "true" ]; then
       if [ "${ANNOTATION_TYPES}" != "" ]; then
          ANNOTATION_TYPES="${ANNOTATION_TYPES},"
       fi
       ANNOTATION_TYPES="${ANNOTATION_TYPES}exon"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_ESTIMATE_COUNTS_OTHER}" == "true" ]; then
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
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `

    ANNOTATION_PATH=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_ANNOTATIONS_ANNOTATIONS_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

  ANNOTATION_SOURCE=""
  if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_ANNOTATION_SOURCE}" == "GENE_EXON_OTHER" ]; then
    # gene exon annotation file.
    ANNOTATION_SOURCE="${ANNOTATION_PATH}/exon-annotations.tsv"
  else
    # CNV annotation file.
    ANNOTATION_SOURCE="${ANNOTATION_PATH}/cnv-annotations.tsv"
  fi
}

. ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_FILES_PARALLEL_SCRIPT}
. ${RESOURCES_R_SHELL_SCRIPT}
. ${RESOURCES_EDGE_R_SCRIPT_SETUP}

function plugin_alignment_analysis_combine {
   set -x
   set -T
   RESULT_FILE=$1
   shift
   PART_RESULT_FILES=$*

   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_NUM_TOP_HITS}
   Q_VALUE_THRESHOLD=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_Q_VALUE_THRESHOLD}
   NORMALIZATION_FACTORS_METHOD=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_NORMALIZATION_FACTORS_METHOD}
   DISPERSION_METHOD=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_DISPERSION_METHOD}
   FILTERING=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_EDGE_R_ARTIFACT_FILTERING}


   # following function from DIFF_EXP_GOBY/parallel.sh
   run_fdr

   if [ $? -eq 0 ]; then
      GENE_OUT_FILENAME=gene-counts.stats.tsv
      EXON_OUT_FILENAME=exon-counts-stats.tsv

      # Extract gene part of file:
      head -1 ${OUT_FILENAME} >${GENE_OUT_FILENAME}
      grep GENE ${OUT_FILENAME} |cat >>${GENE_OUT_FILENAME}
      wc -l ${GENE_OUT_FILENAME}

      HAS_GENES=`cat ${GENE_OUT_FILENAME} | wc -l`


      # Extract exon part of file:
      head -1 ${OUT_FILENAME} >${EXON_OUT_FILENAME}
      grep EXON ${OUT_FILENAME} 1>>${EXON_OUT_FILENAME}
      ls -lat
      HAS_EXONS=`cat ${EXON_OUT_FILENAME} | wc -l `

    cp ${SGE_O_WORKDIR}/sampleToGroups.tsv .
    SAMPLE_GROUP_MAPPING="sampleGroupMapping=sampleToGroups.tsv"

      # Run EdgeR on gene and exon input:

      if [ "$HAS_GENES" != "1" ]; then

        EDGE_R_OUTPUT="output=gene-stats.tsv mdsPlotOutput=mds.png smearPlotOutput=smear.png"
        run_R -f ${RESOURCES_EDGE_R_SCRIPT_R_SCRIPT} --slave --quiet --no-restore --no-save --no-readline --args input=${GENE_OUT_FILENAME} ${EDGE_R_OUTPUT} ${SAMPLE_GROUP_MAPPING} elementType=GENE normalizationMethod=${NORMALIZATION_FACTORS_METHOD} dispersionMethod=${DISPERSION_METHOD} filterFlag=${FILTERING}
        cp ${GENE_OUT_FILENAME} counts-table.tsv
      fi
      if [ "$HAS_EXONS" != "1" ]; then

        EDGE_R_OUTPUT="output=exon-stats.tsv mdsPlotOutput=mds.png smearPlotOutput=smear.png"
        run_R -f ${RESOURCES_EDGE_R_SCRIPT_R_SCRIPT} --slave --quiet --no-rest/ore --no-save --no-readline --args input=${EXON_OUT_FILENAME} ${EDGE_R_OUTPUT} ${SAMPLE_GROUP_MAPPING} elementType=EXON normalizationMethod=${NORMALIZATION_FACTORS_METHOD} dispersionMethod=${DISPERSION_METHOD} filterFlag=${FILTERING}
        cp ${EXON_OUT_FILENAME} counts-table.tsv
      fi
      if [ "${HAS_GENES}" != "1" ] && [ "${HAS_EXONS}" != "1" ]; then

        run_goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
          --q-threshold 1.0 \
          gene-stats.tsv exon-stats.tsv  \
          ${COLUMNS} \
          --output stats.tsv

      elif [ "${HAS_GENES}" != "1" ]; then
        mv gene-stats.tsv stats.tsv
      elif [ "${HAS_EXONS}" != "1" ]; then
        mv exon-stats.tsv stats.tsv
      fi
   fi
}
