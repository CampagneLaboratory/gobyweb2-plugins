    # Installation script for ENSEMBL_GENOMES
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'TOPLEVEL' )
            ORGANISM=$3
            GENOME_REFERENCE_ID=$4
            ENSEMBL_RELEASE=$5
            echo "Organism=${ORGANISM} Reference-build=${GENOME_REFERENCE_ID} ENSEMBL_RELEASE=${ENSEMBL_RELEASE}"
            ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`
            # get genome from local NFS, or wget from ensembl servers if not present
            if [ -e ~gobyweb/genomes/${ORGANISM}/${GENOME_REFERENCE_ID}/${ENSEMBL_RELEASE}/ ]; then

                cp ~gobyweb/genomes/${ORGANISM}/${GENOME_REFERENCE_ID}/${ENSEMBL_RELEASE}/*  copied.dna.toplevel.fa.gz
            else

                wget ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/${ORG_LOWERCASE}/dna/\*.dna.toplevel.fa.gz
                # Check the NFS cache again, put the wget file in the cache if no other process did it:
                if [ ! -e ~gobyweb/genomes/${ORGANISM}/${GENOME_REFERENCE_ID}/${ENSEMBL_RELEASE}/ ]; then
                    mkdir -p ~gobyweb/genomes/${ORGANISM}/${GENOME_REFERENCE_ID}/${ENSEMBL_RELEASE}/
                    cp  *.dna.toplevel.fa.gz  ~gobyweb/genomes/${ORGANISM}/${GENOME_REFERENCE_ID}/${ENSEMBL_RELEASE}/dna.toplevel.fa.gz
                fi
            fi
            cp *dna.toplevel.fa.gz ${installation_path}/genome-toplevel.fasta.gz
            if [ -e ${installation_path}/genome-toplevel.fasta.gz ]; then

             exit 0
            else

             exit 1
            fi
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

    echo get_attribute_values for ID=${id}

    # get environment variables for GobyWeb job from SGE work directory:
    . ${SGE_O_WORKDIR}/constants.sh

    BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
    ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
    echo >>${out} "organism=${ORGANISM}"
    echo >>${out} "reference-build=${BUILD_NUMBER}"
    echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"

    echo "Printing result from ${out}:"
    cat ${out}
    return 0
}