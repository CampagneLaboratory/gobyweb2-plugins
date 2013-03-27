. ${RESOURCES_SCALA_SHELL_SCRIPT}

SPLICING_PLAN_FILENAME="${SGE_O_WORKDIR}/splicejunctioncoverage-groups.tsv"
function plugin_alignment_analysis_split {
	local SPLICING_PLAN_RESULT=$2
    cp ${SGE_O_WORKDIR}/splicejunctioncoverage-groups.tsv ${SPLICING_PLAN_RESULT}
    # the split plan was already built as the file splicejunctioncoverage-groups.tsv

}

# This function return the number of parts in the slicing plan. It returns zero if the alignments could not be split.
function plugin_alignment_analysis_num_parts {
	local SPLICING_PLAN_FILE=$1

	if [ $? -eq 0 ]; then
		echo `wc -l < ${SPLICING_PLAN_FILENAME}`
    else

        echo 0
    fi
}

. ${RESOURCES_GOBY_SHELL_SCRIPT}

function plugin_alignment_analysis_process {

	local SPLICING_PLAN_FILENAME=$1
	local CURRENT_PART=$2

	#get basenames for this part
	local SOURCE_LINE=`sed -n ${CURRENT_PART}p < ${SPLICING_PLAN_FILENAME}`
	local SOURCE_PATH=`echo ${SOURCE_LINE} |cut -d " " -f 1`  # first token is scp path for Splice file.
	local SOURCE_GROUP=`echo ${SOURCE_LINE} |cut -d " " -f 2` # Second token is the group

	#copy over splice data
	scp "${SOURCE_PATH}" ${TAG}-stats-${CURRENT_PART}.tsv
	#	 ${SGE_O_WORKDIR}/split-results/"${CURRENT_PART}_SpliceJunctionCoverage.tsv"
    local FORCE_GOBY_SPLICE_USE=${PLUGINS_ALIGNMENT_ANALYSIS_SPLICING_DIFF_EXP_FORCE_GOBY_SPLICE_USE}

	if [ $? -eq 0 -a ! ${FORCE_GOBY_SPLICE_USE} ]; then
		echo "found splice junction coverage data for sample "
	else

		echo "found splice junction coverage data could not be found, extracting from alignment files"

        copy_local  entries ${CURRENT_PART}
        copy_local  header  ${CURRENT_PART}
        copy_local  index   ${CURRENT_PART}
        extract_splicing_info    ${CURRENT_PART} ${TAG}-stats-${CURRENT_PART}.tsv
	fi
}

function copy_local {

  local FILE_TYPE=$1
  local CURRENT_PART=$2
  cp ${SPLICING_PLAN_FILENAME} ./${FILE_TYPE}-filenames.tsv
  sed -i s/-SpliceJunctionCoverage-all.tsv/\\\*.${FILE_TYPE}/g ./${FILE_TYPE}-filenames.tsv

  # Scp the first path from the TSV file, at line SOURCE_LINE to source/${TAG}/
  local SOURCE_LINE=`sed -n ${CURRENT_PART}p < ./${FILE_TYPE}-filenames.tsv`
  local SOURCE_PATH=`echo ${SOURCE_LINE} |cut -d " " -f 1`  # first token is scp path for Splice file.
  # Copy the alignment file of the appropriate type to the current directory:
  scp ${SOURCE_PATH} ./
}

function extract_splicing_info {
   local CURRENT_PART=$1
   local DESTINATION_FILE=$2
   run-goby  ${PLUGIN_NEED_COMBINE_JVM} extract-splicing-events --min-mapping-quality 255  *.entries -o  ${DESTINATION_FILE}
  # scala ${PLUGIN_NEED_COMBINE_JVM} ${SGE_O_WORKDIR}/goby.jar  ${PLUGINS_ALIGNMENT_ANALYSIS_SPLICING_DIFF_EXP_FILES_EXTRACT_SPLICING_SCRIPT} \
   #    *.entries  > ${DESTINATION_FILE}

}

# This function is called after the analysis parts have finished executing
. ${RESOURCES_R_SHELL_SCRIPT}
. constants.sh
. auto-options.sh
function plugin_alignment_analysis_combine {

   RESULT_FILE=$1
   shift
   PART_RESULT_FILES=$*
    # On some weird systems, touching the directory helps update its file list:
    touch ${SGE_O_WORKDIR}/split-results/*

    scala ${PLUGIN_NEED_COMBINE_JVM} ${SGE_O_WORKDIR}/goby.jar  ${PLUGINS_ALIGNMENT_ANALYSIS_SPLICING_DIFF_EXP_FILES_PROCESS_SCRIPT} \
                                                ${PART_RESULT_FILES} > counts.tsv
    cp counts.tsv ${SGE_O_WORKDIR}/
    cp ${SGE_O_WORKDIR}/sampleGroups.tsv .
    dieUponError  "Cannot copy sample to group mapping information to local directory."

    # The following for debugging:
    #head -10000 counts.tsv >1.tsv
    #cp 1.tsv counts.tsv


    # Run DESeq or EdgeR to estimate p-values:
    if [ "${PLUGINS_ALIGNMENT_ANALYSIS_SPLICING_DIFF_EXP_STAT_ENGINE}" == "DESEQ" ]; then

     run-R -f ${RESOURCES_DESEQ_SCRIPT_R_SCRIPT} --slave --quiet --no-restore --no-save \
           --no-readline --args input=counts.tsv elementType=SPLICE \
           output=out1.tsv graphOutput=.png sampleGroupMapping=sampleGroups.tsv
    else

     run-R -f ${RESOURCES_EDGE_R_SCRIPT_R_SCRIPT} --slave --quiet --no-restore --no-save \
                --no-readline --args input=counts.tsv elementType=SPLICE \
                output=out1.tsv mdsPlotOutput=mds.png smearPlotOutput=smear.png sampleGroupMapping=sampleGroups.tsv \
                normalizationMethod=TMM dispersionMethod=tagwise filterFlag=true
      
    fi

    dieUponError  "Calling statistics with R script failed."

    # 08/22/2012 - no need to do this R script now preserves rownames - there are no guarantees that the
    # order of the input counts.tsv file is maintained in out1.tsv resulting in mixing up genomic locations
    # if this cut and paste operation is carried out.
    # R does not preserve rownames if they they contain some characters. Rather than trying to guess
    # what characters are allowed, we just copy and paste the columns here (order is preserved):
    #cut -f 1 counts.tsv >ids.tsv
    #cut -f 2- out1.tsv >data.tsv
    #paste ids.tsv data.tsv >out2.tsv

    cp out1.tsv  ${SGE_O_WORKDIR}/
    scala ${PLUGIN_NEED_COMBINE_JVM} ${SGE_O_WORKDIR}/goby.jar  ${PLUGINS_ALIGNMENT_ANALYSIS_SPLICING_DIFF_EXP_FILES_POST_PROCESS_SCRIPT} \
        ${REFERENCE_DIRECTORY}/exon-annotations.tsv \
        out1.tsv > junctions.tsv

}