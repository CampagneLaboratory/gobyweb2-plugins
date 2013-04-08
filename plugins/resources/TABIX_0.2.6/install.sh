# Installation script for Tabix version 0.2.6
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            VERSION="0.2.6"
            ${RESOURCES_FETCH_URL_SCRIPT} "http://downloads.sourceforge.net/project/samtools/tabix/tabix-0.2.6.tar.bz2?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fsamtools%2Ffiles%2Ftabix%2F&ts=1364484044&use_mirror=iweb" tabix-${VERSION}.tar.bz2
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
