# Installation script for Groovy resource
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'DISTRIBUTION' )

            VERSION="2.0.6"
            ${RESOURCES_FETCH_URL_SCRIPT} https://dl.bintray.com/groovy/maven/groovy-binary-${VERSION}.zip groovy-binary-${VERSION}.zip
            cp groovy-binary-${VERSION}.zip ${installation_path}/
            (cd ${installation_path} && rm -rf groovy-binary-${VERSION} && unzip groovy-binary-${VERSION}.zip && mv groovy-${VERSION}/* ./ && rm groovy-binary-${VERSION}.zip)
            exit 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}

