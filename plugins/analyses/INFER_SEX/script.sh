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
  ls -1 $* >${SPLICING_PLAN_RESULT}

}

# This function return the number of parts in the slicing plan. It returns zero if the alignments could not be split.
function plugin_alignment_analysis_num_parts {
   SPLICING_PLAN_FILE=$1

   if [ $? -eq 0 ]; then

        echo `cat ${SPLICING_PLAN_FILE} |wc -l `

   else
        echo 0
   fi

}

function plugin_alignment_analysis_process {

   SLICING_PLAN_FILENAME=$1
   ARRAY_JOB_INDEX=$2
   shift
   shift
   NEXT_LINE=$[${ARRAY_JOB_INDEX} + 1]
   ALIGNMENT=`head -${NEXT_LINE} ${SLICING_PLAN_FILENAME} | tail -1`

     # Note that we override the grid jvm flags to request only 4Gb:
     run_goby ${PLUGIN_NEED_PROCESS_JVM} infer-sex  ${ALIGNMENT} -o ${TAG}-stats-${ARRAY_JOB_INDEX}.tsv

     dieUponError  "Infer sex failed for part, sub-task ${CURRENT_PART} failed."

}

function plugin_alignment_analysis_combine {

   RESULT_FILE=stats.vcf.gz
   shift
   PART_RESULT_FILES=$*

   OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_OUTPUT_FORMAT}
   NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_GOBY_NUM_TOP_HITS}

        run_goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
          --tsv \
          ${PART_RESULT_FILES}  \
          --output ${TMPDIR}/${TAG}-inferred-sex.tsv

}