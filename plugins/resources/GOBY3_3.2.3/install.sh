# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="3.2.3"
            ${RESOURCES_FETCH_URL_SCRIPT} http://chagall.med.cornell.edu/goby/releases/archive/release-goby_${VERSION}.tgz
            tar -zxvf release-goby_${VERSION}.tgz
            cd release-goby_${VERSION}
            unzip goby_${VERSION}-bin.zip
            export JAVA_HOME=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}
            mv goby-${VERSION}/goby.jar ${installation_path}
            mkdir ${installation_path}/models
            mv goby-${VERSION}/models/* ${installation_path}/models
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