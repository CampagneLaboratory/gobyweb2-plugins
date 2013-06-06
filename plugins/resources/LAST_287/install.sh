# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )

                (
                VERSION="287"
                    ${RESOURCES_FETCH_URL_SCRIPT} http://last.cbrc.jp/last/index.cgi/archive/${VERSION}.zip ${VERSION}.zip

                    unzip ${VERSION}.zip

                    cd  last-${VERSION}
                    echo "${VERSION}" > src/version.hh
                    make
                    mkdir ${installation_path}/bin/
                    cp src/lastdb ${installation_path}/bin/
                    cp src/lastal ${installation_path}/bin/
                    cp src/lastex ${installation_path}/bin/
                    cp -r scripts ${installation_path}/
                    cp -r examples ${installation_path}/
                )
            if [ -e ${installation_path}/bin/lastal ]; then
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
