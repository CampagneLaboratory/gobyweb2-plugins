# Installation script for SALMON
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )

                (
                VERSION="0.5.0"
                    set -x
                    yum install -y autoconf268.noarch
                    ${RESOURCES_FETCH_URL_SCRIPT} https://github.com/COMBINE-lab/salmon/archive/v${VERSION}.tar.gz Salmon-src-dist.tar.gz
                    ls -ltrh
                    gzip  -c -d Salmon-src-dist.tar.gz |tar -xvf -
                    cd  salmon-${VERSION}

# the doc told us to, but we don't or the file CMakeLists.txt cannot be found
#                    mkdir build
#                    cd build
                    rm -fr  CMakeCache.txt CMakeFiles
                    export BOOST_INCLUDEDIR=${RESOURCES_ARTIFACTS_BOOST_LIB_BINARIES}/include
                    export BOOST_LIBRARYDIR=${RESOURCES_ARTIFACTS_BOOST_LIB_BINARIES}/lib
                    cmake  -DFETCH_BOOST=TRUE -DCMAKE_INSTALL_PREFIX=${installation_path}
                    make
                    make install
# Create a setup script, which will set LD_LIBRARY_PATH:
cat >${installation_path}/setup.sh <<EOT
export LD_LIBRARY_PATH=${installation_path}/lib:\${LD_LIBRARY_PATH}
EOT
            chmod +x ${installation_path}/setup.sh
                )
            if [ -e ${installation_path}/bin/salmon ]; then
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
