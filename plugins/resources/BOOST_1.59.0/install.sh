# Installation script for SALMON
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )

                (
                VERSION="0.4.2"
                    set -x

                    ${RESOURCES_FETCH_URL_SCRIPT} http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz/download boost.tar.gz
                    gzip  -c -d  boost.tar.gz|tar -xvf -
                    cd boost_1_59_0
                    ./bootstrap.sh --prefix=${installation_path}
                    # install will `fail' with status 1 because we don't install optional components. See http://stackoverflow.com/questions/12906829/failed-updating-58-targets-when-trying-to-build-boost-what-happened
                    # ignore the 1 status.
                    ./b2 install || echo "success"
                )
            if [ -e ${installation_path}/lib ]; then
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
