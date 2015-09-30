# Installation script for SALMON
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )

                (
                VERSION="0.4.2"
                    set -x
                    yum install -y autoconf268.noarch
                    # relying on fetching boost is not very stable: the Fossies archive could ban your IP if you download boost too often
                    # This is not a problem when gobyweb caches downloads, but inside a container the cache does not contain the archive..
                    # yum install -y boost.i686
                    ${RESOURCES_FETCH_URL_SCRIPT} http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.gz/download boost.tar.gz
                    gzip  -c -d  boost.tar.gz|tar -xvf -
                    cd boost_1_59_0
                    ./bootstrap.sh --prefix=${installation_path}
                    ./b2 install
                    cd ..
                    ${RESOURCES_FETCH_URL_SCRIPT} https://github.com/COMBINE-lab/salmon/archive/v${VERSION}.tar.gz Salmon-src-dist.tar.gz
                    ls -ltrh
                    gzip  -c -d Salmon-src-dist.tar.gz |tar -xvf -

                    cd  salmon-${VERSION}

# the doc told us to, but we don't or the file CMakeLists.txt cannot be found
#                    mkdir build
#                    cd build
                    rm -fr  CMakeCache.txt CMakeFiles
                    cmake  -DCMAKE_INSTALL_PREFIX=${installation_path} -DBOOST_ROOT=${RESOURCES_ARTIFACT_BOOST_LIB_BINARIES} # -DFETCH_BOOST=TRUE
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
