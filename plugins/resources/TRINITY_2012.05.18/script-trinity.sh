

#args list
#1: Path to input compact-reads file
#2: Output file
#3: Other options for trinity

function run_trinity {

	local INPUT=$1
	local OUTPUT=$2
	shift
	shift
	
	
	local OTHER_OPTIONS=$*

	#trinity doesnt support compact-reads, convert to fastq
	local TEMPFILE=`mktemp readsXXXX`
	run-goby 4g compact-to-fasta --output-format fasta --input ${INPUT} --output ${TEMPFILE}

	#run trinity on converted file
	local TEMPDIR=`mktemp -d trinity_outXXXX`
	${RESOURCES_ARTIFACTS_TRINITY_TRINITY_2012_05_18}/Trinity.pl --bflyHeapSpaceMax 4G --seqType fa --kmer_method inchworm --single ${TEMPFILE} --CPU 8 --output ${TEMPDIR} ${OTHER_OPTIONS}
	
	#copy trinity output file to specified location
	cp ${TEMPDIR}/Trinity.fasta ${OUTPUT}
	
}




