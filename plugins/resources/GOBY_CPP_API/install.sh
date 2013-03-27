# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'LIBRARIES' )
                set -x
                . ${RESOURCES_ARTIFACTS_PROTOBUF_CPP_LIBRARIES}/setup.sh
                (
                VERSION="2.1.2"
                    ${RESOURCES_FETCH_URL_SCRIPT} http://chagall.med.cornell.edu/goby/releases/archive/release-goby_${VERSION}/goby-cpp.zip

                    unzip goby-cpp.zip

                    cd  2.1.2/cpp/
                    chmod +x autogen.sh
                    ./autogen.sh
                    ./configure --prefix=${installation_path}
                    make
                    make install
                )
cat >${installation_path}/setup.sh <<EOT
export LOCAL_LIB=${installation_path}/
export PKG_CONFIG_PATH=/usr/lib/pkgconfig:${installation_path}/lib/pkgconfig:\${PKG_CONFIG_PATH}
export PATH=${installation_path}/bin:\${PATH}
export LD_LIBRARY_PATH=${installation_path}/lib:\${LD_LIBRARY_PATH}
EOT
            chmod +x ${installation_path}/setup.sh
            if [ -e ${installation_path}/lib/libgoby.a ]; then
                return 0
            else
               return 127
            fi

            ;;



        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}
