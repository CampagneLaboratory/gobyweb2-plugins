

function setup_trinity {

	if [ ! -e Trinity.pl ]
	then
		tar -x -f ${RESOURCES_TRINITY_TRINITY_TAR}
	fi
	
}

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
	
	setup_trinity
	
	#run trinity on converted file
	local TEMPDIR=`mktemp -d trinity_outXXXX`
	./Trinity.pl --bflyHeapSpaceMax 4G --seqType fa --kmer_method inchworm --single ${TEMPFILE} --CPU 8 --output ${TEMPDIR} ${OTHER_OPTIONS}
	
	#copy trinity output file to specified location
	cp ${TEMPDIR}/Trinity.fasta ${OUTPUT}
	
}




