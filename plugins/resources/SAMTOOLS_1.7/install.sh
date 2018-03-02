# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
            VERSION="1.7"
            ${RESOURCES_FETCH_URL_SCRIPT} https://github.com/samtools/samtools/releases/download/${VERSION}/samtools-${VERSION}.tar.bz2 samtools-${VERSION}.tar.bz2
            bunzip2  samtools-${VERSION}.tar.bz2
            tar -xvf samtools-${VERSION}.tar
            cd samtools-${VERSION}
            make
            cp samtools ${installation_path}/
            cp bcftools/bcftools ${installation_path}/
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
