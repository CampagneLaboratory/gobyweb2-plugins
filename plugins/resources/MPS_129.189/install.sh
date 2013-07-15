# Installation script for MPS 3.0 EAP
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'DISTRIBUTION' )
            VERSION="3.0"
            BUILD="EAP-129.189"
            ${RESOURCES_FETCH_URL_SCRIPT} http://download.jetbrains.com/mps/30/MPS-${VERSION}-${BUILD}.tar.gz MPS-${VERSION}-${BUILD}.tar.gz
            gzip -c -d MPS-${VERSION}-${BUILD}.tar.gz |tar -xvf -
            cp -r MPS\ ${VERSION}/* ${installation_path}/
            return 0
            ;;
         'LANGUAGES' )
                return 0;
            ;;
        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
