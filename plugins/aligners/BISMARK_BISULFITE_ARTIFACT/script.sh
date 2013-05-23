#!/bin/sh

# This script expects the following variables to be defined:

# PAIRED_END_ALIGNMENT = true|false
# READS = reads file

# INDEX_DIRECTORY = directory that contains the indexed database
# INDEX_PREFIX = name of the indexed database to search

# ${RESOURCES_ILLUMINA_ADAPTERS_FILE_PATH} = path to adapters.txt, obtained from the ILLUMINA_ADAPTERS resource

# ALIGNER_OPTIONS = any GSNAP options the end-user would like to set

function plugin_align {

    OUTPUT=$1
    BASENAME=$2
    
    # set the number of threads to the number of cores available on the server divided by 4:
    ALIGNER_OPTIONS="${ALIGNER_OPTIONS} -p $((`grep physical /proc/cpuinfo | grep id | wc -l` / 4))"

	#set other aligner options
	ALIGNER_OPTIONS="${ALIGNER_OPTIONS} --fastq --bowtie2 --path_to_bowtie ${RESOURCES_ARTIFACTS_BOWTIE2_ARTIFACT_BINARIES} "
	if [ "${LIB_PROTOCOL_PRESERVE_STRAND}" == "false" ]; then
		ALIGNER_OPTIONS="${ALIGNER_OPTIONS} --non_directional"
	fi
	
	#bismark is for methylation analysis so non-bisulfite doesnt make sense
    if [ "${BISULFITE_SAMPLE}" == "false" ]; then
        false
        dieUponError "only bisulfite samples are supported, alignment failed"
    fi

	#take out the slice of the reads file we are currently working with
	SPLIT_READS='split-reads.compact-reads'
	
	goby reformat-compact-reads  --start-position=${START_POSITION} --end-position=${END_POSITION}  ${READS_FILE} -o ${SPLIT_READS}
    dieUponError "reformat reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"


	#trim adapters off of reads
	TRIMMED_READS='trimmed-reads.compact-reads'

    goby trim  --input ${SPLIT_READS} --output ${TRIMMED_READS} --complement \
    		--adapters  ${RESOURCES_ILLUMINA_ADAPTERS_FILE_PATH}  --min-left-length 4
    dieUponError "trim reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
    
    ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
    INDEX_DIRECTORY=$(eval echo \${RESOURCES_ARTIFACTS_BISMARK_ARTIFACT_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})/

    #convert to fastq and set input options
    FASTQ_READS='trimmed-reads-sanger.fastq'
	INPUT_OPTIONS="${INDEX_DIRECTORY}"
    if [ "${PAIRED_END_ALIGNMENT}" == "true" ]; then
        # Bismark expects fastq format files with sequence line in one row - guarantee with --fasta-line-length parameter
    	goby compact-to-fasta --output-format fastq --quality-encoding Sanger --fasta-line-length 1000 --input ${TRIMMED_READS} --output ${FASTQ_READS}-1 --pair-output ${FASTQ_READS}-2
    	
    	INPUT_OPTIONS="${INPUT_OPTIONS} -1 ${FASTQ_READS}-1 -2 ${FASTQ_READS}-2"
    	
    else
        goby compact-to-fasta --output-format fastq --quality-encoding Sanger --fasta-line-length 1000 --input ${TRIMMED_READS} --output ${FASTQ_READS}
    	
    	INPUT_OPTIONS="${INPUT_OPTIONS} ${FASTQ_READS}"
        
    fi
    
    dieUponError "converting to fastq format failed"    


    # Get bismark from the Bismark resource artifact:
    BISMARK_EXEC_PATH=${RESOURCES_ARTIFACTS_BISMARK_ARTIFACT_SCRIPTS}/bismark

    # Run the Bismark alignment
    ${BISMARK_EXEC_PATH} ${ALIGNER_OPTIONS} ${INPUT_OPTIONS}
    dieUponError "alignment of reads with bismark failed"
    
    ls -ltr
    
    ${RESOURCES_SAMTOOLS_EXEC_PATH} view -Sbu "${FASTQ_READS}_bt2_bismark.sam" |
    	${RESOURCES_SAMTOOLS_EXEC_PATH} sort -o - output |
    	${RESOURCES_SAMTOOLS_EXEC_PATH} calmd - ${INDEX_DIRECTORY}/*.fa > 'md-alignment.sam'
    dieUponError "adding MD tag to output failed"
    
    ls -ltr
    head 'md-alignment.sam'
    
    #convert to compact alignment format
	goby sam-to-compact --preserve-all-tags --quality-encoding Sanger --input 'md-alignment.sam' --output ${OUTPUT}
	dieUponError "converting to compact alignment format failed"
	
	ls -ltr
	
	#copy summary report back
	mkdir -p ${SGE_O_WORKDIR}/split-results/mapping-reports
	cp "${FASTQ_READS}_bt2_Bismark_mapping_report.txt" "${SGE_O_WORKDIR}/split-results/${TAG}-mapping-report-${CURRENT_PART}.txt"
    dieUponError "copying back mapping report failed"
    
}

# This function is called after the alignment slices have been combined into one final output.
# It is called with three arguments, the basename of the alignment (present in the directory where this function is
# invoked), the reads filename (full path), the tag useful to create an output.

function plugin_alignment_combine {
    TAG=$1
    READS=$2
    BASENAME=$3

	#tarball mapping-records and copy them back
	mkdir mapping-reports
	cp ${SGE_O_WORKDIR}/split-results/mapping-reports/* mapping-reports/
	tar -zvcf mapping-reports.tar.gz mapping-reports
	
	cp mapping-reports.tar.gz ${RESULT_DIR}/${TAG}-${BASENAME}-mapping-reports.tar.gz

}
