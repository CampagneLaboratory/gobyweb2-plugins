# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2
    set -x
    case ${id} in

        'GOBY' )

            VERSION="3.0.0"

            git clone git@bitbucket.org:campagnelaboratory/goby.git
            git checkout tags/3.0.0
            cd goby
            ant -f build.xml jar
            mv goby.jar ${installation_path}
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