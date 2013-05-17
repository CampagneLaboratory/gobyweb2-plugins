# Installation script for Bedtools

function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            VERSION="2.17.0"
            ${RESOURCES_FETCH_URL_SCRIPT} https://bedtools.googlecode.com/files/BEDTools.v${VERSION}.tar.gz

            tar -zxvf BEDTools.v${VERSION}.tar.gz
            cd bedtools-${VERSION}
            make clean
            make all
            cp bin/* ${installation_path}/

            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
