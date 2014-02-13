# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# IS_TRANSCRIPT = whether alignments were done against a cDNA reference.
# GROUPS_DEFINITION = description of the groups, in the format group-1=sample_i,sample_j/group-2=sample_k,..
# COMPARE_DEFINITION
# ANNOTATION_FILE = file describing annotations in the Goby annotation format.
# ANNOTATION_TYPES = gene|exon|other, specifies the kind of annotations to calculate counts for.
# USE_WEIGHTS_DIRECTIVE = optional, command line flags to have Goby annotation-to-counts adjust counts with weigths.

function eval {
EVAL=raw-counts
}

function setupWeights {

   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_WEIGHT_ADJUSTMENT}" == "NONE" ]; then

       USE_WEIGHTS_DIRECTIVE=" "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_WEIGHT_ADJUSTMENT}" == "GC_CONTENT" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights gc --adjust-gc-bias ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_BIAS_ADJUSTMENT_FORMULA} "

   elif [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_WEIGHT_ADJUSTMENT}" == "HEPTAMERS" ]; then

       USE_WEIGHTS_DIRECTIVE="--use-weights heptamers "
   else
     dieUponError "weight adjustment  not supported: ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_WEIGHT_ADJUSTMENT}"
   fi

}

function setupAnnotationTypes {
   ANNOTATION_TYPES=""
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_ESTIMATE_COUNTS_GENE}" == "true" ]; then

       ANNOTATION_TYPES="${ANNOTATION_TYPES}gene"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_ESTIMATE_COUNTS_EXON}" == "true" ]; then
       if [ "${ANNOTATION_TYPES}" != "" ]; then
          ANNOTATION_TYPES="${ANNOTATION_TYPES},"
       fi
       ANNOTATION_TYPES="${ANNOTATION_TYPES}exon"
   fi
   if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_ESTIMATE_COUNTS_OTHER}" == "true" ]; then
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
  if [ "${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_ANNOTATION_SOURCE}" == "GENE_EXON_OTHER" ]; then
    # gene exon annotation file.
    ANNOTATION_SOURCE="${REFERENCE_DIRECTORY}/exon-annotations.tsv"
  else
    # CNV annotation file.
    ANNOTATION_SOURCE="${REFERENCE_DIRECTORY}/cnv-annotations.tsv"
  fi
}

. ${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_FILES_PARALLEL_SCRIPT}
. ${RESOURCES_R_SHELL_SCRIPT}

function plugin_alignment_analysis_combine {
   set -x
   set -T
   RESULT_FILE=$1
   shift
   PART_RESULT_FILES=$*

   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_NUM_TOP_HITS}
   Q_VALUE_THRESHOLD=${PLUGINS_ALIGNMENT_ANALYSIS_DIFF_EXP_DESEQ_Q_VALUE_THRESHOLD}

   # following function from DIFF_EXP_DESEQ/parallel.sh
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

      # Run DESeq on gene and exon input:

      if [ "$HAS_GENES" != "1" ]; then

        DESEQ_OUTPUT="output=gene-stats.tsv graphOutput=.png"
        run_R -f ${RESOURCES_DESEQ_SCRIPT_R_SCRIPT} --slave --quiet --no-restore --no-save --no-readline \
            --args ${DESEQ_OUTPUT} input=${GENE_OUT_FILENAME} elementType=GENE  ${SAMPLE_GROUP_MAPPING}


      fi
      if [ "$HAS_EXONS" != "1" ]; then

        DESEQ_OUTPUT="output=exon-stats.tsv graphOutput=.png"
        run_R -f ${RESOURCES_DESEQ_SCRIPT_R_SCRIPT} --slave --quiet --no-restore --no-save --no-readline \
             --args ${DESEQ_OUTPUT} input=${EXON_OUT_FILENAME} elementType=EXON ${SAMPLE_GROUP_MAPPING}    \

      fi
      if [ "${HAS_GENES}" != "1" ] && [ "${HAS_EXONS}" != "1" ]; then

        run-goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
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
