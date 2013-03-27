# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            VERSION="0.2.6"
            ${RESOURCES_FETCH_URL_SCRIPT} http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2/download tabix-0.2.6.tar.bz2
            bunzip2  tabix-${VERSION}.tar.bz2
            tar -xvf tabix-${VERSION}.tar
            cd tabix-${VERSION}
            make
            cp tabix ${installation_path}/
            cp bgzip ${installation_path}/
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
