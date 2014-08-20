
# Installation script for Trimmomatic, see http://www.usadellab.org/cms/?page=trimmomatic
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            VERSION="0.32"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${VERSION}.zip  Trimmomatic.zip
            unzip Trimmomatic.zip
            cd Trimmomatic-${VERSION}
            cp trimmomatic-${VERSION}.jar ${installation_path}/
            return 1
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
