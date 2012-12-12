
. ${RESOURCES_GOBY_SHELL_SCRIPT}

#args list
#1: Path to the reads file
#2: Path to the alignment basename
#3: index of the alignment part being processed
#4: Directory where split results should be stored until combine is called
#3: Output file

function extract_unmatched_in_plugin_align {

    local READS_FILE=$1           # The complete reads file.
    local ALIGNMENT_BASENAME=$2
    local NUMBER_OF_PARTS=$3
    local CURRENT_PART=$4
    local SPLIT_DIR=$5   # Where split results should be stored until combine is called

    #export unmatched reads

    extract_unmatched_reads "${READS_FILE}" "${ALIGNMENT_BASENAME}" "${ALIGNMENT_BASENAME}-unmatched.compact-reads" ${NUMBER_OF_PARTS} ${CURRENT_PART}
    dieUponError "extract_unmatched_reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

    #copy slice to shared filesystem

    mkdir -p ${SPLIT_DIR}
    cp "${ALIGNMENT_BASENAME}-unmatched.compact-reads" "${SPLIT_DIR}/unmatched${CURRENT_PART}.compact-reads"
}

function extract_unmatched_reads {

    local READS_FILE=$1
    local ALIGNMENT=$2
    local OUTPUT_FILE=$3
    local NUMBER_OF_PARTS=$4
    local CURRENT_PART=$5


    # export unmatched reads
    # create filter file
    goby alignment-to-read-set "${ALIGNMENT}" --non-matching-reads --non-ambiguous-reads -s unmatched
    dieUponError "creating filter failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

    # filter out unmatched reads
    goby reformat-compact-reads "${READS_FILE}" -f "${ALIGNMENT}-unmatched.filter" -o "${OUTPUT_FILE}"
    dieUponError "reformat reads failed, sub-task ${CURRENT_PART} of ${NUMBER_OF_PARTS}, failed"

}

function combine_splits {
    local SPLIT_DIR=$1
    local RESULT_DIR=$2
    local ALIGNMENT_BASENAME=$3
    # copy files to local file system

    mkdir "unmatched-slices"
    cp ${SPLIT_DIR}/unmatched*.compact-reads "./unmatched-slices"

    # concat files together

    local UNMATCHED_SLICE_FILENAMES=""
    for file in `ls ./unmatched-slices/*`
    do
        UNMATCHED_SLICE_FILENAMES="${UNMATCHED_SLICE_FILENAMES} $file"
    done

    goby concatenate-compact-reads --quick-concat --output "${ALIGNMENT_BASENAME}-unmatched.compact-reads" ${UNMATCHED_SLICE_FILENAMES}
    dieUponError "concatenate unmapped reads failed, sub-task combine, failed"

    # send them to $RESULT_DIR with final filename

    cp "${ALIGNMENT_BASENAME}-unmatched.compact-reads" "${RESULT_DIR}/${ALIGNMENT_BASENAME}-unmatched.compact-reads"
    dieUponError "copying unmapped reads failed, sub-task combine, failed"
    # TODO: make it possible to add -unmatched.compact-reads to output of align plugin job
}
