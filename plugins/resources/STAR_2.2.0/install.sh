# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'EXECUTABLE' )
            wget ftp://ftp2.cshl.edu/gingeraslab/tracks/STARrelease/Alpha/STAR_2.2.0g.tgz
            gzip -c -d  STAR_2.2.0g.tgz |tar -xvf -
            (cd STAR_2.2.0g; make)
            cp STAR_2.2.0g/STAR ${installation_path}/
            return 0
            ;;

        'INDEX' )
            ORGANISM=$3
            BUILD_NUMBER=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"

            STAR_DIR=${RESOURCES_ARTIFACTS_STAR_EXECUTABLE}
            INDEX_DIR=index
            mkdir -p ${INDEX_DIR}
            GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_GENOMES_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
            FAI_INDEXED_GENOME_DIR=$(eval echo \${RESOURCES_ARTIFACTS_FAI_INDEXED_GENOMES_SAMTOOLS_FAI_INDEX_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})


            NUM_THREADS=`grep physical  /proc/cpuinfo |grep id|wc -l`
            SJDB=gencode.v14.annotation.gtf.sjdb
            SPLICE_SITES_OPTION=""
            if [ "$ORGANISM" = "HOMO_SAPIENS" ]; then

                rm -f ${SJDB}
                wget ftp://ftp2.cshl.edu/gingeraslab/tracks/STARrelease/STARgenomes/SpliceJunctionDatabases/${SJDB}
                sed -e 's/^chr//' ${SJDB}  >${SJDB}.fixed
                SPLICE_SITES_OPTION=" --sjdbFileChrStartEnd ${SJDB}.fixed --sjdbOverhang 49"
            elif [ "$ORGANISM" = "MUS_MUSCULUS" ]; then
                rm -f ${SJDB}
                wget ftp://ftp2.cshl.edu/gingeraslab/tracks/STARrelease/STARgenomes/SpliceJunctionDatabases/${SJDB}
                cp ${SJDB}  >${SJDB}.fixed
                SPLICE_SITES_OPTION=" --sjdbFileChrStartEnd ${SJDB}.fixed --sjdbOverhang 49"
            fi
            #INPUT_FASTA_NO_GZ=genome.fasta

            INPUT_FASTA_NO_GZ=${FAI_INDEXED_GENOME_DIR}/genome-toplevel.fasta
            nice ${STAR_DIR}/STAR --runMode genomeGenerate --genomeDir ${INDEX_DIR} \
            --genomeFastaFiles ${INPUT_FASTA_NO_GZ}  --runThreadN ${NUM_THREADS} ${SPLICE_SITES_OPTION}
            STATUS=$?
            if [ ${STATUS} != 0 ]; then
             return ${STATUS}
            fi
            # Keep only the first ID to make short chromosome names:
            cut -f 1 -d " " ${INDEX_DIR}/chrName.txt > tmp-chrName
            mv tmp-chrName  ${INDEX_DIR}/chrName.txt
            cp -r ${INDEX_DIR} ${installation_path}/
            echo "Finished indexing Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID}"
            exit 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}


function get_attribute_values() {

    id=$1
    out=$2

    echo "get_attribute_values for ID=${id}"
    set +xv
    echo # get environment variables for GobyWeb job from SGE work directory:
    . ${SGE_O_WORKDIR}/constants.sh

    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
    ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
    echo >>${out} "organism=${ORGANISM}"
    echo >>${out} "reference-build=${BUILD_NUMBER}"
    echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"
    set -xv
    echo "Printing result from ${out}:"
    cat ${out}
    exit 0
}