# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'DISTRIBUTION' )

            VERSION="2.0.6"
            ${RESOURCES_FETCH_URL_SCRIPT} http://dist.groovy.codehaus.org/distributions/groovy-binary-${VERSION}.zip
            unzip groovy-binary-${VERSION}.zip
            cp -r groovy-${VERSION}/* ${installation_path}/

            exit 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}

