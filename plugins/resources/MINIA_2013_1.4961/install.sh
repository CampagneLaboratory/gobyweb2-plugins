# Installation script for Minia
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'EXECUTABLE' )
            VERSION=1.4961
            ${RESOURCES_FETCH_URL_SCRIPT} http://minia.genouest.org/files/minia-${VERSION}.tar.gz
            gzip -c -d  minia-1.4961.tar.gz |tar -xvf -
            (cd minia-${VERSION}; make)
            cp minia-${VERSION}/minia ${installation_path}/
            return 0
            ;;


        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}
