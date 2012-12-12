#!/bin/sh
# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true|false
# READS = reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search

# ${RESOURCES_STAR_EXEC_PATH} = path to gsnap, obtained from the GSNAP_GOBY resource

# ALIGNER_OPTIONS = any STAR options the end-user would like to set

. ${RESOURCES_GOBY_SHELL_SCRIPT}
. ${RESOURCES_SCALA_SHELL_SCRIPT}
. ${RESOURCES_EXTRACT_NONMATCHED_SHELL_SCRIPT}
# . constants.sh
# . auto-options.sh
function plugin_align {

     OUTPUT=$1
     BASENAME=$2
     # set the number of threads to the number of cores available on the server:
     NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
     MIN_SCORE=$((${INPUT_READ_LENGTH}*70/100))
     if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
         MIN_SCORE=$((${MIN_SCORE}*2))
     fi
     echo "Aligning with minScore= ${MIN_SCORE}"
     local SHARED_MEM=NoSharedMemory
     local SHARED_MEM=LoadAndRemove
     ALIGNER_OPTIONS="${ALIGNER_OPTIONS}  --genomeLoad ${SHARED_MEM} --genomeDir ${INDEX_DIRECTORY} --runThreadN ${NUM_THREADS} --outFilterScoreMin ${MIN_SCORE} --outFilterMatchNmin ${MIN_SCORE}"
                           #        --genomeLoad LoadAndRemove

     cd ${TMPDIR}
     goby reformat-compact-reads  --start-position=${START_POSITION} --end-position=${END_POSITION}  ${READS_FILE} -o small-reads.compact-reads
     dieUponError "reformat reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

     local ambiguityOption="--outFilterMultimapNmax ${PLUGINS_ALIGNER_STAR_GOBY_AMBIGUITY_THRESHOLD}"

     if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then

         # Convert the compact-reads slice to FASTQ for paired-end data:  (note that STAR 2.1.1 requires the sequence on one line, so we use a max 10,000 bp per line!)
         run-goby ${PLUGIN_NEED_ALIGN_JVM} compact-to-fasta  -n 10000 -i small-reads.compact-reads -o 1.fastq -p 2.fastq --output-format fastq
         dieUponError "Convert compact-reads to fastq failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

         nice ${RESOURCES_STAR_EXEC_PATH}  ${ALIGNER_OPTIONS} ${PLUGINS_ALIGNER_STAR_GOBY_ALIGNER_OPTIONS} ${ambiguityOption} --readFilesIn 1.fastq 2.fastq
         RETURN_STATUS=$?
     else
         # Convert the compact-reads slice to FASTQ for single-end data:
         run-goby ${PLUGIN_NEED_ALIGN_JVM} compact-to-fasta -n 10000 -i small-reads.compact-reads -o reads.fastq --output-format fastq
         dieUponError "Convert compact-reads to fastq failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

         nice ${RESOURCES_STAR_EXEC_PATH}  ${ALIGNER_OPTIONS} ${PLUGINS_ALIGNER_STAR_GOBY_ALIGNER_OPTIONS} ${ambiguityOption} --readFilesIn reads.fastq
         RETURN_STATUS=$?
     fi

     cat Log.out
     cat Log.progress.out
     echo "STARR Finished with status code=${RETURN_STATUS}"

     if [ ! ${RETURN_STATUS} -eq 0 ]; then
            #cp reads.fastq ${SGE_O_WORKDIR}/split-results/reads-${CURRENT_PART}.fastq
            # Failed, no result to copy
            copy_logs align ${CURRENT_PART} ${NUMBER_OF_PARTS}
            ${QUEUE_WRITER} --tag ${TAG} --index ${CURRENT_PART} --job-type job-part --status ${JOB_PART_FAILED_STATUS} --description "STAR alignment failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
            exit
     fi

     #RESULT_DIR= directory on shared filesystem, send output files to $RESULT_DIR/split-results
     #CURRENT_PART= unique id associated with this part of the job
     mkdir -p ${SGE_O_WORKDIR}/split-results/
     cp  SJ.out.tab ${SGE_O_WORKDIR}/split-results/SpliceJunctionCoverage-${CURRENT_PART}.tsv
     #cp  Aligned.out.sam ${SGE_O_WORKDIR}/split-results/aligned-${CURRENT_PART}.sam

     # ADD MD tags to the sam file with samtools: NB: we sort the BAM file because calmd is terribly slow on non-sorted input
     ${RESOURCES_SAMTOOLS_EXEC_PATH} view -S -b -u Aligned.out.sam |${RESOURCES_SAMTOOLS_EXEC_PATH} sort - sam_sorted
     ${RESOURCES_SAMTOOLS_EXEC_PATH} calmd -u sam_sorted.bam ${REFERENCE_DIRECTORY}/*.fa.gz >Aligned.out.bam
    # cp Aligned.out.bam ${SGE_O_WORKDIR}/split-results/Aligned-${CURRENT_PART}.bam

     # Convert SAM output to Goby:
     run-goby ${PLUGIN_NEED_ALIGN_JVM} sam-to-compact -i Aligned.out.bam -o ${OUTPUT} --read-names-are-query-indices
     dieUponError "SAM conversion to Goby output failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

     if [ "${PLUGINS_ALIGNER_STAR_GOBY_NON_MATCHING}" == "true" ]; then

     	 extract_unmatched_in_plugin_align ${READS} ${OUTPUT} ${NUMBER_OF_PARTS} ${CURRENT_PART} ${SGE_O_WORKDIR}/split-results
     fi
#extra variables:

#${SGE_O_WORKDIR}= directory on shared filesystem, send output files to $RESULT_DIR/split-results
#CURRENT_PART= unique id associated with this part of the job


}

# This function is called after the alignment slices have been combined into one final output.
# It is called with three arguments, the basename of the alignment (present in the directory where this function is
# invoked), the reads filename (full path), the tag useful to create an output.

function plugin_alignment_combine {
    TAG=$1
    READS=$2
    BASENAME=$3
cat > script.awk <<EOT
        BEGIN{
            print("sample\tchromosome\tfirst base of the intron (1-based)\tlast base of the intron (1-based)\tstrand\tintron motif\tannotation\tnumber of uniquely mapping reads crossing the junction\treserved\tmaximum left/right overhang");
            split("non-canonical,GT/AG,CT/AC,GC/AG,CT/GC,AT/AC,GT/AT",motifs,",")
        }
        {
            printf("%s\t",sample)
            for (i=1; i<=NF ;i++) {
                token=\$i
                if (i!=5) {
                    printf("%s",token);
                } else {
                   code=token;
                   # Awk arrays are one-based:
                   printf("%s",motifs[code+1])
                }
                if (i!=NF) printf("\t")
            }
            print("")
        }
EOT

    cat ${SGE_O_WORKDIR}/split-results/SpliceJunctionCoverage-*.tsv |  \
        awk -v sample=${BASENAME} -f script.awk >combined.tsv
   #cp combined.tsv ${SGE_O_WORKDIR}/
    scala ${PLUGIN_NEED_ALIGNMENT_POST_PROCESSING_JVM} goby.jar ${SGE_O_WORKDIR}/CombineSpliceData.scala combined.tsv \
     > ${RESULT_DIR}/${TAG}-SpliceJunctionCoverage-all.tsv
    #${RESOURCES_SAMTOOLS_EXEC_PATH} merge out.bam ${SGE_O_WORKDIR}/split-results/*.bam
    #${RESOURCES_SAMTOOLS_EXEC_PATH} sort out.bam sorted
    #${RESOURCES_SAMTOOLS_EXEC_PATH} index sorted
    #cp sorted.bam*  ${RESULT_DIR}/

    if [ "${PLUGINS_ALIGNER_STAR_GOBY_NON_MATCHING}" == "true" ]; then

            combine_splits  ${SGE_O_WORKDIR}/split-results ${RESULT_DIR}/ ${BASENAME}
    fi

}
