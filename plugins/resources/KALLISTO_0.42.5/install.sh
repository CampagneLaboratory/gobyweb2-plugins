function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
            VERSION="0.42.5"
            set -x
            ${RESOURCES_FETCH_URL_SCRIPT} https://github.com/pachterlab/kallisto/archive/v$VERSION.tar.gz kallisto.tar.gz
            tar -zxvf kallisto.tar.gz
            cd kallisto-$VERSION
            mkdir build
            cd build
            cmake -DCMAKE_INSTALL_PREFIX=${installation_path}  ..
            make
            make install
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
