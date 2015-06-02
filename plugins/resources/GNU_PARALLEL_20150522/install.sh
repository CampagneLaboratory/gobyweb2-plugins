# Installation script for Gnu Parallel

function plugin_install_artifact {

    id=$1
    installation_path=$2

    set -e
    VERSION=20150522
    case ${id} in

        'BINARIES' )

       ${RESOURCES_FETCH_URL_SCRIPT} http://ftp.gnu.org/gnu/parallel/parallel-${VERSION}.tar.bz2 parallel.tar.bz2

            bunzip2 parallel.tar.bz2
            tar -xvf parallel.tar
            cd parallel-${VERSION}
            ./configure --prefix=${installation_path} && make && make install
            if [ -e ${installation_path}/bin/parallel ]; then
             exit 0
            else
             exit 1
            fi
            ;;


        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
