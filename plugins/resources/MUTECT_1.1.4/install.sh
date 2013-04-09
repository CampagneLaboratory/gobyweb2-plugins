# Installation script for MuTect version 1.1.4
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            VERSION="1.1.4"
            ${RESOURCES_FETCH_URL_SCRIPT} "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/muTect-1.1.4-bin.zip" muTect-${VERSION}-bin.zip
            unzip muTect-${VERSION}-bin.zip
            cp muTect-${VERSION}.jar muTect.jar
            cp muTect.jar ${installation_path}/
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}