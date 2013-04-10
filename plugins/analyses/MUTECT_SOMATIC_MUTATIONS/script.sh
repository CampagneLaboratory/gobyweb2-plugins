
. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_alignment_analysis_split {

  NUMBER_OF_PARTS=$1
  SPLICING_PLAN_RESULT=$2
  shift
  shift
  goby suggest-position-slices \
          --number-of-bytes 50000000 \
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

   ORG=` echo ${ORGANISM} | tr [:lower:] [:upper:]  `
   if [! "${ORG}" == "HOMO_SAPIENS" ]; then
        dieUponError "Invalid organism ${ORG}"
   fi

   BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}' | tr [:lower:] [:upper:] `
   ENSEMBL_RELEASE=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'| tr [:lower:] [:upper:] `
   SEQUENCE_CACHE_DIR=$(eval echo \${RESOURCES_ARTIFACTS_GOBY_INDEXED_GENOMES_SEQUENCE_CACHE_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
   INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORG}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})

   WINDOW_LIMITS=`awk -v arrayJobIndex=${ARRAY_JOB_INDEX} '{ if (lineNumber==arrayJobIndex) print " -s "$3" -e "$6; lineNumber++; }' ${SLICING_PLAN_FILENAME}`
   echo "Fetching covariates' details at ${PLUGINS_ALIGNMENT_ANALYSIS_MUTECT_COVARIATE_INFO_URL}"
   ${RESOURCES_FETCH_URL_SCRIPT} ${PLUGINS_ALIGNMENT_ANALYSIS_MUTECT_COVARIATE_INFO_URL} covariates.tsv ${JOB_DIR}/results/
   echo "Fetching b37_cosmic_v54_120711.vcf "
   ${RESOURCES_FETCH_URL_SCRIPT} http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf b37_cosmic_v54_120711.vcf ${JOB_DIR}/
   echo "Fetching dbsnp_132_b37.leftAligned.vcf.gz "
   ${RESOURCES_FETCH_URL_SCRIPT} http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/dbsnp_132_b37.leftAligned.vcf.gz dbsnp_132_b37.leftAligned.vcf.gz ${JOB_DIR}/

   gunzip ${JOB_DIR}/dbsnp_132_b37.leftAligned.vcf.gz

   tail -n +2 covariates.tsv > ${TMPDIR}/covariates_no_header.tsv

   declare -a GermlineDetails='()'
   declare -a SomaticDetails='()'

   #loads the samples' details and checks if they exist
    while IFS=$'\t' read sample patient gender type kind tissue parents
    do
        echo reading patient_id=$patient, kind_of_sample=$kind
        kind_lc=`echo ${kind} | tr [:upper:] [:lower:]`
        case "${kind_lc}" in
                "germline")
                        FILES=(`ls source/${sample}.*`)
                        if [ ${#FILES[@]} -gt 0 ]; then
                            GermlineDetails["id-${patient}"]="$sample"
                        else
                            echo "ERROR: Germline sample was not provided for patient ${patient}"
                        fi
                        ;;
                "somatic")
                         FILES=(`ls source/${sample}.*`)
                         if [ ${#FILES[@]} -gt 0 ]; then
                            SomaticDetails["id-${patient}"]="$sample"
                          else
                            echo "ERROR: Somatic sample was not provided for patient ${patient}"
                        fi
                        ;;
                *)
                        echo "ERROR: Invalid kind of sample found ${kind} for patient ${patient}"
        esac
    done < ${TMPDIR}/covariates_no_header.tsv

    #for each pair...
    for id in "${!GermlineDetails[@]}"
    do
        echo $id ${GermlineDetails[id]} ${SomaticDetails[id]}

        if [[ ${SomaticDetails[id]} ]]; then

               #1) concatenate-alignments and produce a slice of Goby Alignments (GA) + add-read-origin-info to GA
               echo "Concatenating and slicing germline sample ${GermlineDetails[id]} "
               run-goby concatenate-alignments \
                    ${WINDOW_LIMITS} \
                    ---add-read-origin-info
                    --output ${TMPDIR}/germline-ca-${GermlineDetails[id]} \
                    ${TMPDIR}/${GermlineDetails[id]}.entries

               echo "Concatenating and slicing somatic sample ${SomaticDetails[id]} "
               run-goby concatenate-alignments \
                    ${WINDOW_LIMITS} \
                    ---add-read-origin-info
                    --output ${TMPDIR}/somatic-ca-${SomaticDetails[id]} \
                    ${TMPDIR}/${SomaticDetails[id]}.entries

               #2)convert the GA to BAM

               #convert the Germline slice
               echo "Converting Goby alignment ${TMPDIR}/germline-ca-${GermlineDetails[id]} to BAM"
               run-goby ${PLUGIN_NEED_PROCESS_JVM} compact-to-sam \
                    --genome ${SEQUENCE_CACHE_DIR}/random-access-genome \
                    --output ${TMPDIR}/germline-ca-${GermlineDetails[id]}.bam \
                    ${TMPDIR}/germline-ca-${GermlineDetails[id]}

               #convert the Somatic slice
               echo "Converting Goby alignment ${TMPDIR}/somatic-ca-${SomaticDetails[id]} to BAM"
               run-goby ${PLUGIN_NEED_PROCESS_JVM} compact-to-sam \
                    --genome ${SEQUENCE_CACHE_DIR}/random-access-genome \
                    --output ${TMPDIR}/somatic-ca-${SomaticDetails[id]}.bam \
                    ${TMPDIR}/somatic-ca-${SomaticDetails[id]}

               #3) run MuTect
               echo "Running MuTect with \
                 --input_file:normal ${TMPDIR}/germline-ca-${GermlineDetails[id]}.bam  \
                 --input_file:tumor ${TMPDIR}/somatic-ca-${SomaticDetails[id]}.bam  \
                 --out ${id}-stats.tsv"
               ${MUTECT_EXEC_PATH} \
                    --analysis_type MuTect \
                    --input_file:normal ${TMPDIR}/germline-ca-${GermlineDetails[id]}.bam  \
                    --input_file:tumor ${TMPDIR}/somatic-ca-${SomaticDetails[id]}.bam \
                    --out ${id}-stats.tsv  \
                    --reference_sequence ${INDEXED_GENOME_DIR}/*toplevel.fasta
                    --cosmic ${JOBDIR}/b37_cosmic_v54_120711.vcf
                    --dbsnp ${JOBDIR}/dbsnp_132_b37.leftAligned.vcf
                    #--intervals <intervals_to_process>  #may need to convert the WINDOW_LIMITS to mutect syntax
                    #--coverage_file <coverage.wig.txt>


        else
                echo "ERROR: Tumor sample was not provided for patient ${id}"
        fi
    done

}


function plugin_alignment_analysis_combine {
   RESULT_FILE=stats.vcf.gz
   shift
   PART_RESULT_FILES=$*

    # use goby in fdr mode. To concatenate TSVs produced by Mutect in a single file
    echo "Concating TSVs with Goby in results.tsv"
    run-goby ${PLUGIN_NEED_COMBINE_JVM} fdr \
        --output results.tsv \
        *-stats.tsv

}