function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'LINUX_BINARIES' )
            VERSION="8u141"
            VERSIONALT="1.8.0_141"
            set -x
            ${RESOURCES_FETCH_URL_SCRIPT} http://download.oracle.com/otn-pub/java/jdk/8u141-b15/336fa29ff2bb4ef291e347e091f7f4a7/jdk-${VERSION}-linux-x64.tar.gz "gpw_e24=http%3A%2F%2Fwww.oracle.com%2F;oraclelicense=accept-securebackup-cookie"
            tar -zxvf java.tar.gz -C ${installation_path}
            mv ${installation_path}/jdk$VERSIONALT/* ${installation_path}/
            cd ${installation_path}/
            if [ -e ${installation_path}/bin/java ]; then
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