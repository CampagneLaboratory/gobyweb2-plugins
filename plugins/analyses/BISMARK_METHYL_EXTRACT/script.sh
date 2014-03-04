#!/bin/bash
# COPYRIGHT_MAY_GO_HERE
# This script expects the following variables to be defined:


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
pwd
. ${RESOURCES_GOBY_SHELL_SCRIPT}
. ${JOB_DIR}/plugin-constants.sh


function plugin_alignment_analysis_split {
    echo ;
}

# This function return the number of parts in the slicing plan. It returns zero if the alignments could not be split.
function plugin_alignment_analysis_num_parts {

	echo "${PLUGIN_NUM_SPLITS}"

}

function plugin_alignment_analysis_process {

	local BASENAME=${PLUGIN_BASENAME[${CURRENT_PART}]}
	
	#convert compact-alignment back to sam
	local SAM_FILE_UNSORTED=${BASENAME}-unsorted.sam
	local SAM_FILE=${BASENAME}.sam

	BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
    ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
    SEQUENCE_CACHE_DIR=$(eval echo \${RESOURCES_ARTIFACTS_GOBY_INDEXED_GENOMES_SEQUENCE_CACHE_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
    INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

	goby compact-to-sam --genome ${SEQUENCE_CACHE_DIR}/random-access-genome --output ${SAM_FILE_UNSORTED} ${ENTRIES_DIRECTORY}/${BASENAME}

    ${RESOURCES_SAMTOOLS_EXEC_PATH} view -Sbu ${SAM_FILE_UNSORTED} |
    	${RESOURCES_SAMTOOLS_EXEC_PATH} sort -o - output |
    	${RESOURCES_SAMTOOLS_EXEC_PATH} view -h - -o ${SAM_FILE}

	#run bismark methylation extractor

    BISMARK_OPTIONS=''
	if [ ${PLUGIN_PAIRED[${CURRENT_PART}]} == "true" ]; then
	    BISMARK_OPTIONS="${BISMARK_OPTIONS} --paired-end"
	else
        BISMARK_OPTIONS="${BISMARK_OPTIONS} --single-end"
    fi

	${RESOURCES_BISMARK_METHYL_EXTRACT} ${BISMARK_OPTIONS} ${SAM_FILE}
	
	ls -althr

	mkdir -p ${SGE_O_WORKDIR}/tempout/${CURRENT_PART}
	cp * ${SGE_O_WORKDIR}/tempout/${CURRENT_PART}/
	
	${SGE_O_WORKDIR}/ProcessAggregator.groovy CpG_OT_${SAM_FILE}.txt > ${TAG}-CpG_OT-${CURRENT_PART}.tsv
	${SGE_O_WORKDIR}/ProcessAggregator.groovy CHG_OT_${SAM_FILE}.txt > ${TAG}-CHG_OT-${CURRENT_PART}.tsv
	${SGE_O_WORKDIR}/ProcessAggregator.groovy CHH_OT_${SAM_FILE}.txt > ${TAG}-CHH_OT-${CURRENT_PART}.tsv

	${SGE_O_WORKDIR}/ProcessAggregator.groovy CpG_OB_${SAM_FILE}.txt > ${TAG}-CpG_OB-${CURRENT_PART}.tsv
	${SGE_O_WORKDIR}/ProcessAggregator.groovy CHG_OB_${SAM_FILE}.txt > ${TAG}-CHG_OB-${CURRENT_PART}.tsv
	${SGE_O_WORKDIR}/ProcessAggregator.groovy CHH_OB_${SAM_FILE}.txt > ${TAG}-CHH_OB-${CURRENT_PART}.tsv
}

function plugin_alignment_analysis_combine {

	local OUTPUT_FILE="methylation.tsv"
	shift
	local PART_RESULT_FILES=$*
	
	ln -s ${PART_RESULT_FILES} .

	COLUMNS='Chromosome\tPosition\tStrand\tContext'

    for sample in "${PLUGIN_BASENAME[@]}"
    do
        COLUMNS="${COLUMNS}\t${sample}-MC\t${sample}-UC\t${sample}-MR"
    done

    echo -e ${COLUMNS} > tmp-output.tsv

    ${SGE_O_WORKDIR}/CombineAggregator.groovy ${TAG}-CpG_OT-*.tsv | groovy -p -a '\t' -e 'split[0..1].join("\t") + "\t+\tCpG\t" + split[2..-1].join("\t")' >>  tmp-output.tsv
    ${SGE_O_WORKDIR}/CombineAggregator.groovy ${TAG}-CpG_OB-*.tsv | groovy -p -a '\t' -e 'split[0..1].join("\t") + "\t-\tCpG\t" + split[2..-1].join("\t")' >>  tmp-output.tsv

    ${SGE_O_WORKDIR}/CombineAggregator.groovy ${TAG}-CHG_OT-*.tsv | groovy -p -a '\t' -e 'split[0..1].join("\t") + "\t+\tCHG\t" + split[2..-1].join("\t")' >>  tmp-output.tsv
    ${SGE_O_WORKDIR}/CombineAggregator.groovy ${TAG}-CHG_OB-*.tsv | groovy -p -a '\t' -e 'split[0..1].join("\t") + "\t-\tCHG\t" + split[2..-1].join("\t")' >>  tmp-output.tsv

    ${SGE_O_WORKDIR}/CombineAggregator.groovy ${TAG}-CHH_OT-*.tsv | groovy -p -a '\t' -e 'split[0..1].join("\t") + "\t+\tCHH\t" + split[2..-1].join("\t")' >>  tmp-output.tsv
    ${SGE_O_WORKDIR}/CombineAggregator.groovy ${TAG}-CHH_OB-*.tsv | groovy -p -a '\t' -e 'split[0..1].join("\t") + "\t-\tCHH\t" + split[2..-1].join("\t")' >>  tmp-output.tsv
    sed 's/\tNA/\t/g' tmp-output.tsv > ${OUTPUT_FILE}

}