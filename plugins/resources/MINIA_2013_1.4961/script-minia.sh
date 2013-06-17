


#args list
#1: Path to input compact-reads file
#2: Output file

function run_minia {

	local INPUT=$1
	local OUTPUT=$2
	shift
	shift
	

	#minia doesnt support compact-reads, convert to fasta
	local TEMPFILE=`mktemp readsXXXX`
	run-goby 4g compact-to-fasta --output-format fasta --fasta-line-length 10000 --input ${INPUT} --output ${TEMPFILE}
	
	#run minia on converted file
	#./minia input kmer_length min_abundance estimated_genome_size prefix
	${RESOURCES_ARTIFACTS_MINIA_EXECUTABLE}/minia ${TEMPFILE} 25 3 3000000000 unused
	
	# copy minia output file to specified location
	cp contigs.fa ${OUTPUT}
	
}




