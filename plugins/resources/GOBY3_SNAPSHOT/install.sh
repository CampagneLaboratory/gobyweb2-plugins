# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="3.0.3-SNAPSHOT"
            ${RESOURCES_FETCH_URL_SCRIPT} https://www.dropbox.com/s/afprsn3zbrza0tl/goby-${VERSION}.zip
            unzip goby-${VERSION}.zip
            export JAVA_HOME=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}
            mv goby-${VERSION}/* ${installation_path}
            chmod +x ${installation_path}/goby
            if [ -e ${installation_path}/goby.jar ]; then
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