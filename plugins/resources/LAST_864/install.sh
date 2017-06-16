# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )

                (
                VERSION="864"
                    ${RESOURCES_FETCH_URL_SCRIPT} http://last.cbrc.jp/last-${VERSION}.zip ${VERSION}.zip
                    unzip ${VERSION}.zip

                    cd  last-${VERSION}
                    # enable multi-threading compilation.
                    make install prefix=${installation_path}
                    #mkdir ${installation_path}/bin/
                    #cp src/lastdb ${installation_path}/bin/
                    #cp src/lastal ${installation_path}/bin/
                    #cp src/last-pair-probs ${installation_path}/bin/
                    #cp src/last-split ${installation_path}/bin/
                    cp -r scripts ${installation_path}/
                    cp -r examples ${installation_path}/
                    # older versions of last defined scripts with the .sh or .py extension. Create links to
                    # maintain backward compatibility:
                    cd ${installation_path}/scripts/
                    ln -s  last-map-probs last-map-probs.py
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
