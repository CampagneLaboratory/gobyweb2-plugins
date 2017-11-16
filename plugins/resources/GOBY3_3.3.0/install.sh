# Installation script for GOBY version 3.2.6
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="3.3.0"
            ${RESOURCES_FETCH_URL_SCRIPT} http://chagall.med.cornell.edu/goby/releases/release-goby_${VERSION}/goby.zip
            unzip goby.zip
            export JAVA_HOME=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}
            mv goby-${VERSION}/goby ${installation_path}
            mv goby-${VERSION}/*.jar ${installation_path}
            mkdir ${installation_path}/models
            mkdir ${installation_path}/config
            mv goby-${VERSION}/models/* ${installation_path}/models
            mv goby-${VERSION}/config/* ${installation_path}/config
            chmod +x ${installation_path}/*
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