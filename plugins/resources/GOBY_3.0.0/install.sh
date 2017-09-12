# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )

            git clone git://github.com/CampagneLaboratory/genotype-dl-models.git
            cd goby3
            export JAVA_HOME=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}
            mvn install
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