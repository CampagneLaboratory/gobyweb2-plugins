# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
                set -x
                . ${RESOURCES_ARTIFACTS_PROTOBUF_CPP_LIBRARIES}/setup.sh
                (
                VERSION="286"
                    ${RESOURCES_FETCH_URL_SCRIPT} http://last.cbrc.jp/archive/last-${VERSION}.zip

                    unzip last-${VERSION}.zip

                    cd  last-${VERSION}
                    make
                    mkdir ${installation_path}/bin/
                    cp lastdb ${installation_path}/bin/
                    cp lastal ${installation_path}/bin/
                    cp lastex ${installation_path}/bin/
                    cp -r scripts ${installation_path}/
                    cp -r examples ${installation_path}/
                )
            if [ -e ${installation_path}/bin/lastal ]; then
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
