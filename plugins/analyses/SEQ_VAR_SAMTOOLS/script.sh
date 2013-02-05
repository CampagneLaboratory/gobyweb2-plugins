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

function plugin_alignment_analysis_sequential {

       OUTPUT_FORMAT=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_SAMTOOLS_OUTPUT_FORMAT}
       NUM_TOP_HITS=${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_SAMTOOLS_NUM_TOP_HITS}

       ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Start discover with samtools mpileup" --index ${CURRENT_PART} --job-type job-part

        # Build fasta and fasta index file if necessary in /scratchLocal

        #
        # First try to get the reference as the .fa.gz
        #
        REFERENCE_FASTA_FILEPATH=`find ${REFERENCE_DIRECTORY} -name \*.fa.gz`;
        if [ -s ${REFERENCE_FASTA_FILEPATH} ]; then
            # The .fa.gz exists
            if [ ! -s ${REFERENCE_FASTA_FILEPATH}.fai ]; then
                # fa.gz has no samtools index, create it
                nice ${RESOURCES_SAMTOOLS_EXEC_PATH} faidx  ${REFERENCE_FASTA_FILEPATH}
                if [ ! -s ${REFERENCE_FASTA_FILEPATH}.fai ]; then
                    # Couldn't create the index
                    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Reference fasta and index files could not be created on compute node." --index ${CURRENT_PART} --job-type job-part
                    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed" --index ${CURRENT_PART} --job-type job
                    jobFailedEmail
                    exit
                fi
            fi
        else
            #
            # If that fails, create the .fasta from the .compact-reads file
            #
            REFERENCE_COMPACT_FILE=`find ${REFERENCE_DIRECTORY} -name \*.compact-reads |grep -v toplevel`;
            REFERENCE_BASENAME=`basename ${REFERENCE_COMPACT_FILE} .compact-reads`
            REFERENCE_FASTA_FILEPATH=${REFERENCE_DIRECTORY}/${REFERENCE_BASENAME}.fa
            REFERENCE_FASTA_GZ_FILEPATH=${REFERENCE_DIRECTORY}/${REFERENCE_BASENAME}.fa.gz

            if [ ! -s ${REFERENCE_FASTA_FILEPATH} ]; then

                echo Reference file does not exist on node => generate
                ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Reference fasta and index files are being created on compute node." --index ${CURRENT_PART} --job-type job-part

                goby compact-to-fasta -i ${REFERENCE_COMPACT_FILE} -o ${REFERENCE_FASTA_FILEPATH}
                # Compress it to save disc space
                gzip ${REFERENCE_FASTA_FILEPATH}
                REFERENCE_FASTA_FILEPATH=REFERENCE_FASTA_GZ_FILEPATH
                # create the samtools index
                nice ${RESOURCES_SAMTOOLS_EXEC_PATH} faidx  ${REFERENCE_FASTA_FILEPATH}

                if [ ! -s ${REFERENCE_FASTA_FILEPATH} ] || [ ! -s ${REFERENCE_FASTA_FILEPATH}.fai ]; then
                     # The fasta file or its index could not be created. Fail the job:
                     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Reference fasta and index files could not be created on compute node." --index ${CURRENT_PART} --job-type job-part
                     ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed" --index ${CURRENT_PART} --job-type job
                     jobFailedEmail
                     exit
                fi
            fi
       fi
       echo Proceeding to estimate P-value between groups for variants in entries files.
        #  entries in these files must not have path or .bam extension

       echo ${GROUP1_SAMPLE_LIST} | awk '{for (i=1;i<=NF;i++) print $i}' >  group-1-samples.lst
       echo ${GROUP2_SAMPLE_LIST} | awk '{for (i=1;i<=NF;i++) print $i}' >  group-2-samples.lst
       echo group-1 samples:
       cat group-1-samples.lst
       echo group-2 samples:
       cat group-2-samples.lst

       num_samples_in_group_1=`wc -l group-1-samples.lst |awk '{print $1} ' `
       num_samples_in_group_2=`wc -l group-2-samples.lst |awk '{print $1} ' `
       cat group-1-samples.lst group-2-samples.lst >all-samples.lst

       cat all-samples.lst
       # The bam files and bai files have to be in the directory where we run bcftools. This is the only way for the name of the
       # bam file to match the name in all-samples.lst

       DIR=`pwd`
       ( cd  ${SGE_O_WORKDIR}/source/; cp ${ENTRIES_FILES} *.bai $DIR/ )
       cd $DIR
       ls -lth
       ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Running samtools mpileup" --index 1 --job-type job-part

       if [ "${OUTPUT_FORMAT}" == "GROUP_COMPARISONS" ]; then

          nice ${RESOURCES_SAMTOOLS_EXEC_PATH} mpileup -6 -uf ${REFERENCE_FASTA_FILEPATH} ${ENTRIES_FILES}  | ${BCFTOOLS_EXEC_PATH} view  -c -v -g -1 ${num_samples_in_group_1} -s all-samples.lst - > ${SGE_O_WORKDIR}/${TAG}-samtools.vcf

       elif [ "${OUTPUT_FORMAT}" == "GENOTYPES" ]; then

          nice ${RESOURCES_SAMTOOLS_EXEC_PATH} mpileup -6 -uf ${REFERENCE_FASTA_FILEPATH} ${ENTRIES_FILES}  | ${BCFTOOLS_EXEC_PATH} view  -g -v -A - > ${SGE_O_WORKDIR}/${TAG}-samtools.vcf

       else
           ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Analysis type is not supported with samtools mpileup: ${$OUTPUT_FORMAT}" --index 1 --job-type job-part
           ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed" --index ${CURRENT_PART} --job-type job
       fi

       # bcftools does not set status to 0 in case of success? Just check the file existence and size
       if [ -s ${SGE_O_WORKDIR}/${TAG}-samtools.vcf ]; then
         ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_DIFF_EXP_STATUS} --description "Annotating results" --index 1 --job-type job-part


         ${RESOURCES_ANNOTATE_VCF_EXEC_PATH} ${PLUGINS_ALIGNMENT_ANALYSIS_SEQ_VAR_SAMTOOLS_ANNOTATE_VARIATIONS} \
            ${SGE_O_WORKDIR}/${TAG}-samtools.vcf ${SGE_O_WORKDIR}/${TAG}.vcf.gz

        /bin/mkdir -p ${RESULT_DIR}
        /bin/cp ${SGE_O_WORKDIR}/${TAG}.vcf.gz ${RESULT_DIR}/

       fi

}
