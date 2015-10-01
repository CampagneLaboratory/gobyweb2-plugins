# Installation script for SALMON
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'INDEX' )

                (
                VERSION="0.4.2"
                    set -x
                    . ${RESOURCES_ARTIFACTS_SALMON_BINARIES}/setup.sh
                    cp ${RESOURCES_ARTIFACTS_ENSEMBL_TRANSCRIPTS_TOPLEVEL}/*.fa.gz transcripts.fa.gz
                    gunzip  transcripts.fa.gz
                    ${RESOURCES_ARTIFACTS_SALMON_BINARIES}/bin/salmon index -t transcripts.fasta -i transcripts_index
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
