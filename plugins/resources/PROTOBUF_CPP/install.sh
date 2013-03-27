# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'LIBRARIES' )

                cat >${installation_path}/setup.sh <<EOT
export LOCAL_LIB=${installation_path}/
export PKG_CONFIG_PATH=/usr/lib/pkgconfig:${installation_path}/lib/pkgconfig:\${PKG_CONFIG_PATH}
export PATH=${installation_path}/bin:\${PATH}
export LD_LIBRARY_PATH=${installation_path}/lib:\${LD_LIBRARY_PATH}
EOT
                chmod +x ${installation_path}/setup.sh
                ${installation_path}/setup.sh

                mkdir -p ${installation_path}/lib/pkgconfig/
                mkdir -p ${installation_path}/bin/

                (${RESOURCES_FETCH_URL_SCRIPT} http://ftp.gnu.org/gnu/autoconf/autoconf-2.68.tar.gz
                    tar zxvf autoconf-2.68.tar.gz
                    cd autoconf-2.68
                    ./configure --prefix=${installation_path}
                    make
                    make install)

                (${RESOURCES_FETCH_URL_SCRIPT} http://protobuf.googlecode.com/files/protobuf-2.4.1.tar.gz
                      tar zxvf protobuf-2.4.1.tar.gz
                      cd protobuf-2.4.1
                      ./configure --prefix=${installation_path}
                      make
                      make install)

                (${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp.csx.cam.ac.uk/pub/software/programming/pcre/pcre-8.21.tar.gz
                      tar zxvf pcre-8.21.tar.gz
                      cd pcre-8.21
                      ./configure --prefix=${installation_path}
                      make
                      make install
                )

            if [ -e ${installation_path}/lib/libprotobuf.a ]; then
                return 0
            else
               return 127
            fi



            return 0
            ;;



        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}
