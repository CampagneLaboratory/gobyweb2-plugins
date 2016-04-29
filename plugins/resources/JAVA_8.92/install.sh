function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
            VERSION="8u92"
            VERSIONALT="1.80_92"
            set -x
            ${RESOURCES_FETCH_URL_SCRIPT} http://download.oracle.com/otn-pub/java/jdk/$VERSION-b14/jre-$VERSION-linux-x64.tar.gz java.tar.gz
            tar -zxvf java.tar.gz -C /${installation_path}
            cd jre$VERSIONALT
            if [ -e ${installation_path}/bin/kallisto ]; then
               return 0
            else
               return 127
            fi
            return 127
            ;;



        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}
