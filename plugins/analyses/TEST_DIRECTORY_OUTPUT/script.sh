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
 echo "one line" >${SPLICING_PLAN_RESULT}
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

    exit 0
}

function plugin_alignment_analysis_combine {

   RESULT_FILE=out.tsv
   shift
   PART_RESULT_FILES=$*

  mkdir output-directory
  echo "First file" >output-directory/1.sbi
  echo "Second file" >output-directory/1.sbip
  ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Result written to JOB_DIR/results # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part
  echo "Result written to JOB_DIR/results "> out.tsv
}