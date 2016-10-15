# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="1.0.3-SNAPSHOT"
            ${RESOURCES_FETCH_URL_SCRIPT} https://www.dropbox.com/s/hvs1zmmi1l7ko3d/model-training-${VERSION}-bin.jar

            mv model-training-${VERSION}-bin.jar ${installation_path}/model-training-bin.jar

            if [ -e ${installation_path}/model-training-bin.jar ]; then
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