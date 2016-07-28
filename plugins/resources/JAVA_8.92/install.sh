function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
            VERSION="8u45"
            VERSIONALT="1.8.0_45"
            set -x
            ${RESOURCES_FETCH_URL_SCRIPT} http://download.oracle.com/otn-pub/java/jdk/$VERSION-b14/jdk-$VERSION-linux-x64.tar.gz java.tar.gz "gpw_e24=http%3A%2F%2Fwww.oracle.com%2F;oraclelicense=accept-securebackup-cookie"
#           ${RESOURCES_FETCH_URL_SCRIPT} http://www.vim.org/scripts/download_script.php?src_id=6273 java.tar.gz ""  '"gpw_e24=http%3A%2F%2Fwww.oracle.com%2F;"'
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