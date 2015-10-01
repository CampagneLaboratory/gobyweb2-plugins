# Installation script for KALLISTO index
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'INDEX' )

                (
                  ORGANISM=$3
                  BUILD_NUMBER=$4
                  ENSEMBL_RELEASE=$5
                  VERSION="0.42.3"
                  set -x
                  ENSEMBL_TRANSCRIPTS_DIR=$(eval echo \${RESOURCES_ARTIFACTS_ENSEMBL_TRANSCRIPTS_TOPLEVEL_${ORGANISM}_${BUILD_NUMBER}_${ENSEMBL_RELEASE}})
                  cp ${ENSEMBL_TRANSCRIPTS_DIR}/*.fasta.gz transcripts.fasta.gz
                  gunzip  transcripts.fasta.gz
                  ${RESOURCES_ARTIFACTS_KALLISTO_BINARIES}/bin/kallisto index -i transcripts_index transcripts.fasta
                  cp transcripts_index ${installation_path}/
                )
            if [ -e ${installation_path}/transcripts_index ]; then
               return 0
            else
               return 127
            fi
            return 127
            ;;



        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}
