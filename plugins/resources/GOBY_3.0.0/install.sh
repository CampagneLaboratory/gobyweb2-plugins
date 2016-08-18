# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )

            VERSION="3.0.0.3"
            mkdir -p $HOME/.ssh && touch $HOME/.ssh/known_hosts
            ssh-keyscan github.com >> $HOME/.ssh/known_hosts
            ssh-keyscan bitbucket.org >> $HOME/.ssh/known_hosts
            git clone git@bitbucket.org:campagnelaboratory/goby.git
            cd goby
            git checkout tags/${VERSION}
            export JAVA_HOME=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}
            ant -f build.xml jar
            mv goby.jar ${installation_path}
            mkdir ${installation_path}/models
            mv ./models/* ${installation_path}/models
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