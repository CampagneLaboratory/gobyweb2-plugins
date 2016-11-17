# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="1.1-SNAPSHOT"
            ${RESOURCES_FETCH_URL_SCRIPT} https://www.dropbox.com/s/05m0ox95hpxchp0/somatic-${VERSION}-bin.jar

            mv somatic-${VERSION}-bin.jar ${installation_path}/somatic-bin.jar

            if [ -e ${installation_path}/somatic-bin.jar ]; then
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