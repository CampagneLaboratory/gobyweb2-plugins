# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            VERSION="0.1.10"
            ${RESOURCES_FETCH_URL_SCRIPT} http://sourceforge.net/projects/vcftools/files/vcftools_${VERSION}.tar.gz/download vcftools_${VERSION}.tar.gz
            gunzip  vcftools_${VERSION}.tar.gz
            tar -xvf vcftools_${VERSION}.tar
            cd vcftools_${VERSION}
            export PREFIX=${installation_path}
            make

            RETURN_STATUS=$?
            if [ ! $RETURN_STATUS -eq 0 ]; then
                return 127
            else
                return 0
            fi
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
