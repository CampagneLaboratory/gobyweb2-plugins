# Installation script for EdgeR version 2.6.12 which is tied to Bioconductor version 2.11
# edgeR dependencies: Limma
# edgeR script requires: Cairo

function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'BINARIES' )
            . ${RESOURCES_ARTIFACTS_R_BINARIES}/setup.sh
            RUN_R=${RESOURCES_ARTIFACTS_R_BINARIES}/bin/R

            VERSION="1.5-2"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.r-project.org/src/contrib/Cairo_1.5-2.tar.gz
            ${RUN_R} CMD INSTALL Cairo_${VERSION}.tar.gz --library=${installation_path}/

            VERSION="3.12.3"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/limma_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL limma_${VERSION}.tar.gz --library=${installation_path}/

            VERSION="2.6.12"
            ${RESOURCES_FETCH_URL_SCRIPT} http://bioconductor.org/packages/2.10/bioc/src/contrib/edgeR_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL edgeR_${VERSION}.tar.gz --library=${installation_path}/


cat>${installation_path}/setup.sh<<EOT
export R_LIBS=${installation_path}:${R_LIBS}
EOT
       chmod +x ${installation_path}/setup.sh
       # check package installation by checking that all index files exist

       fileList="${installation_path}/Cairo/DESCRIPTION ${installation_path}/limma/DESCRIPTION ${installation_path}/edgeR/DESCRIPTION"

       for file in $fileList; do
           if [ ! -e $file ]; then
               echo "${file} does not exist"
               return 127
           fi
       done
        return 0
       ;;

       *)  echo "Resource artifact id not recognized: "+$id
            return 99
       ;;

    esac

    exit 1

}
