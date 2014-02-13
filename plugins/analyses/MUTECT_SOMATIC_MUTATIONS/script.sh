#!/bin/sh

. ${RESOURCES_GOBY_SHELL_SCRIPT}
. ${PLUGINS_ALIGNMENT_ANALYSIS_MUTECT_SOMATIC_MUTATIONS_FILES_MAPS_IN_BASH3}

function plugin_alignment_analysis_split {

  NUMBER_OF_PARTS=$1
  SPLICING_PLAN_RESULT=$2
  shift
  shift

  # Remove .tmh for all alignments. We don't need them for this type of processing and loading .tmh
  # can slow down loading of alignments whenever we need to make a slice.
  rm ${JOB_DIR}/source/*.tmh

  goby suggest-position-slices \
          --number-of-bytes 50000000 \
          --restrict-per-chromosome \
          --output ${SPLICING_PLAN_RESULT} \
          $*

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


   # Remove .tmh for all alignments. We don't need them for this type of processing and loading .tmh
   # can slow down loading of alignments whenever we need to make a slice.
   rm ${JOB_DIR}/source/*.tmh

   ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
   if [ ! "${ORG}" == "HOMO_SAPIENS" ]; then
        dieUponError "Invalid organism ${ORG}"
   fi

      if [ ! "${PLUGINS_ALIGNMENT_ANALYSIS_MUTECT_SOMATIC_MUTATIONS_COVARIATE_INFO_URL}" == "NONE" ]; then
	    echo "Fetching covariates' details at ${PLUGINS_ALIGNMENT_ANALYSIS_MUTECT_SOMATIC_MUTATIONS_COVARIATE_INFO_URL}"
        ${RESOURCES_FETCH_URL_SCRIPT} ${PLUGINS_ALIGNMENT_ANALYSIS_MUTECT_SOMATIC_MUTATIONS_COVARIATE_INFO_URL} covariates.tsv ${JOB_DIR}/results/
   else
        echo "ERROR: covariates URL not specified"
   fi

   BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
   ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
   SEQUENCE_CACHE_DIR=$(eval echo \${RESOURCES_ARTIFACTS_GOBY_INDEXED_GENOMES_SEQUENCE_CACHE_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
   INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

   WINDOW_LIMITS=`awk -v arrayJobIndex=${ARRAY_JOB_INDEX} '{ if (lineNumber==arrayJobIndex) print " -s "$3" -e "$6; lineNumber++; }' ${SLICING_PLAN_FILENAME}`
   INTERVAL=`awk -v arrayJobIndex=${ARRAY_JOB_INDEX} '{ if (lineNumber==arrayJobIndex) print ""$1":"($2+1)"-"($5+1); lineNumber++; }' ${SLICING_PLAN_FILENAME}`
   IS_PAIRED_END_MATCH="All Compact-Reads files were paired-end = true"
   # Copy cosmic.vcf and dbsnp.vcf to TMPDIR:
   cp ${RESOURCES_ARTIFACTS_MUTECT_HOMO_SAPIENS_DATA_FILES}/* ${TMPDIR}/

   ls -l
   tail -n +2 covariates.tsv > covariates_no_header.tsv

   
   #loads the samples' details and checks if they exist
    while IFS=$'\t' read sample patient gender type kind tissue parents
    do
        echo reading patient_id=$patient, kind_of_sample=$kind
        kind_lc=`echo ${kind} | tr [:upper:] [:lower:]`
        case "${kind_lc}" in
                "germline")
                        FILES=(`ls ${JOB_DIR}/source/${sample}.*`)
                        if [ ${#FILES[@]} -gt 0 ]; then
                            put "GermlineDetails" "${patient}" "$sample"
                        else
                            echo "ERROR: Germline sample was not provided for patient ${patient}"
                        fi
                        ;;
                "somatic")
                         FILES=(`ls ${JOB_DIR}/source/${sample}.*`)
                         if [ ${#FILES[@]} -gt 0 ]; then
                            put "SomaticDetails" "${patient}" "$sample"

                          else
                            echo "ERROR: Somatic sample was not provided for patient ${patient}"
                        fi
                        ;;
                *)
                        echo "ERROR: Invalid kind of sample found ${kind} for patient ${patient}"
        esac
    done < covariates_no_header.tsv

    #for each pair...
    getKeySet "GermlineDetails"
    for id in $keySet
    do
        get "GermlineDetails" $id
        GERMLINE_FILE=$value
        get "SomaticDetails" $id
        SOMATIC_FILE=$value

        #1) concatenate-alignments and produce a slice of Goby Alignments (GA) + add-read-origin-info to GA
               echo "Concatenating and slicing germline sample ${GERMLINE_FILE} "
               run_goby ${PLUGIN_NEED_PROCESS_JVM} concatenate-alignments \
                    ${WINDOW_LIMITS} \
                    --add-read-origin-info \
                    --output ${TMPDIR}/germline-ca-${GERMLINE_FILE} \
                    ${JOB_DIR}/source/${GERMLINE_FILE}.entries
	           dieUponError "Concatenate of goby alignment for germline sample of ${id}, failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
	
		
               echo "Concatenating and slicing somatic sample ${SOMATIC_FILE} "
               run_goby ${PLUGIN_NEED_PROCESS_JVM} concatenate-alignments \
                    ${WINDOW_LIMITS} \
                    --add-read-origin-info \
                    --output ${TMPDIR}/somatic-ca-${SOMATIC_FILE} \
                    ${JOB_DIR}/source/${SOMATIC_FILE}.entries
               dieUponError "Concatenate of goby alignment for somatic sample of ${id}, failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"
	
        #2)convert the GA to BAM

               #convert the Germline slice
               echo "Converting Goby alignment ${TMPDIR}/germline-ca-${GERMLINE_FILE} to BAM"
               run_goby ${PLUGIN_NEED_PROCESS_JVM} compact-to-sam \
                    --genome ${SEQUENCE_CACHE_DIR}/random-access-genome \
                    --output ${TMPDIR}/germline-ca-${GERMLINE_FILE}.bam \
                    ${TMPDIR}/germline-ca-${GERMLINE_FILE}
                dieUponError "Convertion of goby alignment to BAM for germline sample of ${id}, failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

               #convert the Somatic slice
               echo "Converting Goby alignment ${TMPDIR}/somatic-ca-${SOMATIC_FILE} to BAM"
               run_goby ${PLUGIN_NEED_PROCESS_JVM} compact-to-sam \
                    --genome ${SEQUENCE_CACHE_DIR}/random-access-genome \
                    --output ${TMPDIR}/somatic-ca-${SOMATIC_FILE}.bam \
                    ${TMPDIR}/somatic-ca-${SOMATIC_FILE}
                dieUponError "Convertion of goby alignment to BAM  for somatic sample of ${id}, failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

        #3) sort, remove potential PCR duplicates and index Bam files
                ${RESOURCES_SAMTOOLS_EXEC_PATH} sort ${TMPDIR}/germline-ca-${GERMLINE_FILE}.bam ${TMPDIR}/germline-ca-${GERMLINE_FILE}-sorted
                GOBY_OUTPUT=`run_goby 1g cfs \
                 --header-only  ${JOB_DIR}/source/${GERMLINE_FILE}.entries`
                echo ${GOBY_OUTPUT}
                RMDUP_OPTION="-s"
                if echo "$GOBY_OUTPUT" | grep -q "$IS_PAIRED_END_MATCH"; then
                    RMDUP_OPTION="-S"
                fi
                ${RESOURCES_SAMTOOLS_EXEC_PATH} rmdup ${RMDUP_OPTION} ${TMPDIR}/germline-ca-${GERMLINE_FILE}-sorted.bam ${TMPDIR}/germline-ca-${GERMLINE_FILE}-sorted-nodup.bam
                ${RESOURCES_SAMTOOLS_EXEC_PATH} index ${TMPDIR}/germline-ca-${GERMLINE_FILE}-sorted-nodup.bam

                ${RESOURCES_SAMTOOLS_EXEC_PATH} sort ${TMPDIR}/somatic-ca-${SOMATIC_FILE}.bam ${TMPDIR}/somatic-ca-${SOMATIC_FILE}-sorted
                 GOBY_OUTPUT=`run_goby 1g cfs \
                 --header-only  ${JOB_DIR}/source/${SOMATIC_FILE}.entries`
                echo ${GOBY_OUTPUT}
                RMDUP_OPTION="-s"
                if echo "$GOBY_OUTPUT" | grep -q "$IS_PAIRED_END_MATCH"; then
                    RMDUP_OPTION="-S"
                fi
                ${RESOURCES_SAMTOOLS_EXEC_PATH} rmdup ${RMDUP_OPTION} ${TMPDIR}/somatic-ca-${SOMATIC_FILE}-sorted.bam ${TMPDIR}/somatic-ca-${SOMATIC_FILE}-sorted-nodup.bam
                ${RESOURCES_SAMTOOLS_EXEC_PATH} index ${TMPDIR}/somatic-ca-${SOMATIC_FILE}-sorted-nodup.bam


        #4) run MuTect
               echo "Running MuTect..."
               ${RESOURCES_MUTECT_EXEC_PATH} \
                    --analysis_type MuTect \
                    --input_file:normal ${TMPDIR}/germline-ca-${GERMLINE_FILE}-sorted-nodup.bam  \
                    --input_file:tumor ${TMPDIR}/somatic-ca-${SOMATIC_FILE}-sorted-nodup.bam \
                    --out ${TAG}-${id}-stats-${ARRAY_JOB_INDEX}.tsv  \
                    --reference_sequence ${INDEXED_GENOME_DIR}/*toplevel.fasta \
                    --cosmic ${TMPDIR}/cosmic.vcf \
                    --dbsnp ${TMPDIR}/dbsnp.vcf \
                    --intervals ${INTERVAL}  #may need to convert the WINDOW_LIMITS to mutect syntax
                    #--coverage_file <coverage.wig.txt>
                 dieUponError "MuTect analysis on  ${id}, failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

    done

}


function plugin_alignment_analysis_combine {
   RESULT_FILE=results.tsv
   shift
   PART_RESULT_FILES=$*

   for FILE in `ls ${PART_RESULT_FILES}`
    do
       tail -n +2 ${FILE} > ${FILE}_no_header
       mv -f ${FILE}_no_header ${FILE}
    done

    # use goby in fdr mode. To concatenate TSVs produced by Mutect in a single file
    echo "Concating TSVs with Goby in results.tsv"
    run_goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
        --output ${RESULT_FILE} \
        ${PART_RESULT_FILES}
     dieUponError "Concatenation of TSV files with Goby failed" 
}
