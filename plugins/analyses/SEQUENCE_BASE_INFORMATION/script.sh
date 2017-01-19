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
          --number-of-bytes 50000000 \
          --output ${SPLICING_PLAN_RESULT} \
          $*

  if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_ANNOTATIONS}" = "true" ]; then
     URL="${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_ANNOTATION_FILE_URL}"
     echo "Fetching annotations from URL: ${URL}"
     mkdir -p ${JOB_DIR}/results-annotated
     ${RESOURCES_FETCH_URL_SCRIPT} ${URL} annotations.tsv
     cp annotations.tsv  ${JOB_DIR}/results-annotated/

  fi

  if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_GERMLINE_VARMAP}" = "true" ]; then
     URL="${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_VARMAP_FILE_URL}"
     echo "Fetching annotations from URL: ${URL}"
     mkdir -p ${JOB_DIR}/results-annotated
     ${RESOURCES_FETCH_URL_SCRIPT} ${URL} true-genotypes.varmap
     cp true-genotypes.varmap  ${JOB_DIR}/results-annotated/
  fi
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

     if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_GERMLINE_VARMAP}" = "true" ]; then
        VARMAP_OPTION=" -x SequenceBaseInformationOutputFormat:random-seed=3764 \
        -x SequenceBaseInformationOutputFormat:sampling-rate=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_VARMAP_SAMPLING_RATE} \
        -x SequenceBaseInformationOutputFormat:true-genotype-map=${JOB_DIR}/results-annotated/true-genotypes.varmap "
     else
        VARMAP_OPTION=" "
     fi
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
           --diploid ${FORCE_DIPLOID} \
            ${VARMAP_OPTIONS} \
           -x SequenceBaseInformationOutputFormat:sampling-rate=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_SAMPLING_RATE} \
           -x SequenceBaseInformationOutputFormat:random-seed=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_RANDOM_SEED} \
           ${ENTRIES_FILES}
     dieUponError  "Compare sequence variations part, sub-task ${CURRENT_PART} failed."

    if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_GERMLINE_VARMAP}" = "true" ]; then
          # This is a genotype run. The results are already annotated. Copy and exit.
          cp ${TAG}-out-${CURRENT_PART}.sbi  ${JOB_DIR}/split-results/
          cp ${TAG}-out-${CURRENT_PART}.sbip  ${JOB_DIR}/split-results/
          exit 0
    fi
    if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_ANNOTATIONS}" = "false" ]; then
          # Add Java 8 to the path so that variation analysis tools will use it:
          export PATH=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin:${PATH}

          mkdir -p ${JOB_DIR}/split-results
          dieUponError  "cannot create split-results directory. sub-task ${CURRENT_PART} failed."
          cp ${TAG}-out-${CURRENT_PART}.sbi  ${JOB_DIR}/split-results/
          cp ${TAG}-out-${CURRENT_PART}.sbip  ${JOB_DIR}/split-results/

          mkdir -p ${JOB_DIR}/split-mutated
          dieUponError  "cannot create split-mutated directory. sub-task ${CURRENT_PART} failed."
          ${RESOURCES_ARTIFACTS_DLVARIATION_JAR}/bin/mutate.sh ${PLUGIN_NEED_PROCESS_JVM}  \
                                            -i ${TAG}-out-${CURRENT_PART}.sbi -o ${TAG}-mutated-${CURRENT_PART}.sbi \
                                            --strategy ${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_STRATEGY}
          dieUponError  "cannot create mutated file. sub-task ${CURRENT_PART} failed."

          cp ${TAG}-mutated-${CURRENT_PART}.sbi   ${JOB_DIR}/split-mutated/
          cp ${TAG}-mutated-${CURRENT_PART}.sbip  ${JOB_DIR}/split-mutated/
    else

         ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java -cp ${RESOURCES_ARTIFACTS_DLVARIATION_JAR}/somatic/target/somatic-bin-native.jar -Xmx${PLUGIN_NEED_PROCESS_JVM}   \
                                                 org.campagnelab.dl.somatic.tools.CombineWithGoldStandard \
                                                 --sampling-fraction ${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_ANNOTATION_SAMPLING_RATE} \
                                                 --annotations ${JOB_DIR}/results-annotated/annotations.tsv \
                                                 -i ${TAG}-out-${CURRENT_PART}.sbi  \
                                                 -o ${TAG}-out-${CURRENT_PART}-annotated.sbi

        cp ${TAG}-out-${CURRENT_PART}-annotated.sbi* ${JOB_DIR}/results-annotated/
    fi

      ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "End discover-sequence-variations for part # ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part
      # Create an empty TSV file
      touch ${TAG}-out-${ARRAY_JOB_INDEX}.tsv

}

function plugin_alignment_analysis_combine {
    # Add Java 8 to the path so that variation analysis tools will use it:
    export PATH=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin:${PATH}

    RESULT_FILE=out.tsv
    shift
    PART_RESULT_FILES=$*
    if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_ANNOTATIONS}" = "false" ]; then


        if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_GERMLINE_VARMAP}" = "true" ]; then
            RECORDS_PER_BUCKET=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_RECORDS_PER_BUCKET}
            ${RESOURCES_ARTIFACTS_DLVARIATION_JAR}/bin/randomize.sh  ${PLUGIN_NEED_COMBINE_JVM} \
                                            -i ${JOB_DIR}/split-results/*.sbi -o out \
                                            --records-per-bucket ${RECORDS_PER_BUCKET} --chunk-size 50
            dieUponError  "cannot Randomize. sub-task concat failed."
        else
            ${RESOURCES_ARTIFACTS_DLVARIATION_JAR}/bin/concat.sh  ${PLUGIN_NEED_COMBINE_JVM} \
                                            -i  ${JOB_DIR}/split-results/*.sbi -o out
            dieUponError  "cannot QuickConcat. sub-task concat failed."

            RECORDS_PER_BUCKET=${PLUGINS_ALIGNMENT_ANALYSIS_SEQUENCE_BASE_INFORMATION_RECORDS_PER_BUCKET}
            ${RESOURCES_ARTIFACTS_DLVARIATION_JAR}/bin/randomize.sh  ${PLUGIN_NEED_COMBINE_JVM} \
                                            -i  ${JOB_DIR}/split-mutated/*.sbi -o mutated-randomized \
                                            --records-per-bucket ${RECORDS_PER_BUCKET} --chunk-size 50
        dieUponError  "cannot Randomize. sub-task concat failed."
        fi



          # Make a backup of the results in case the web app fails to display (https://bitbucket.org/campagnelaboratory/gobyweb/issues/20/alignment-analysis-job-completed-but-error):
        mkdir -p ${JOB_DIR}/results-copy
        cp  ${TMPDIR}/out.sbi* ${TMPDIR}/mutated-randomized*  ${JOB_DIR}/results-copy/
 else
      echo "Producing annotated file."
      ${RESOURCES_ARTIFACTS_DLVARIATION_JAR}/bin/concat.sh  ${PLUGIN_NEED_COMBINE_JVM} \
                                             -i  ${JOB_DIR}/results-annotated/*-out-*-annotated.sbi -o out-annotated
      dieUponError  "cannot QuickConcat. sub-task concat failed."
      cp out-annotated.sbi* ${JOB_DIR}/results-annotated/
 fi
 ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Result written to JOB_DIR/results-copy ${ARRAY_JOB_INDEX}." --index ${CURRENT_PART} --job-type job-part
 }