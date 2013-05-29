# Installation script for DESEQ version 1.8.3 which is tied to Bioconductor version 2.10
# make sure R_HOME is set to the correct dependent version of R
# get R version library and make installation path?
# DESeq dependencies and their related dependencies: Biogenerics, Biobase, Locfit, Genefilter

# DESeq script dependency Cairo

function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )

            . ${RESOURCES_ARTIFACTS_R_BINARIES}/setup.sh
            RUN_R=${RESOURCES_ARTIFACTS_R_BINARIES}/bin/R
            # Cairo
            VERSION="1.5-2"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.r-project.org/src/contrib/Cairo_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL Cairo_${VERSION}.tar.gz --library=${installation_path}/

            # BiocGenerics
            VERSION="0.2.0"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/BiocGenerics_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL BiocGenerics_${VERSION}.tar.gz --library=${installation_path}

            # Biobase
            VERSION="2.16.0"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/Biobase_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL Biobase_${VERSION}.tar.gz --library=${installation_path}

            # locfit
            VERSION="1.5-9.1"
            ${RESOURCES_FETCH_URL_SCRIPT}  http://cran.fhcrc.org/src/contrib/locfit_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL locfit_${VERSION}.tar.gz --library=${installation_path}

            # DBI
            VERSION="0.2-7"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.fhcrc.org/src/contrib/DBI_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL DBI_${VERSION}.tar.gz --library=${installation_path}

            # RSQLite
            VERSION="0.11.3"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.fhcrc.org/src/contrib/RSQLite_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL RSQLite_${VERSION}.tar.gz --library=${installation_path}

            # IRanges
            VERSION="1.14.4"
            ${RESOURCES_FETCH_URL_SCRIPT}  http://www.bioconductor.org/packages/2.10/bioc/src/contrib/IRanges_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL IRanges_${VERSION}.tar.gz --library=${installation_path}

            # AnnotationDbi
            VERSION="1.18.4"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/AnnotationDbi_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL AnnotationDbi_${VERSION}.tar.gz --library=${installation_path}

            # xtable
            VERSION="1.7-1"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.fhcrc.org/src/contrib/xtable_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL xtable_${VERSION}.tar.gz --library=${installation_path}

            # XML
            VERSION="3.96-1.1"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.fhcrc.org/src/contrib/XML_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL XML_${VERSION}.tar.gz --library=${installation_path}

            # annotate
            VERSION="1.34.1"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/annotate_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL annotate_${VERSION}.tar.gz --library=${installation_path}

            # genefilter
            VERSION="1.38.0"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/genefilter_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL genefilter_${VERSION}.tar.gz --library=${installation_path}

            # RColorBrewer
            VERSION="1.0-5"
            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.fhcrc.org/src/contrib/RColorBrewer_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL RColorBrewer_${VERSION}.tar.gz --library=${installation_path}

            # geneplotter
            VERSION="1.34.0"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/geneplotter_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL geneplotter_${VERSION}.tar.gz --library=${installation_path}

            # DESeq
            VERSION="1.8.3"
            ${RESOURCES_FETCH_URL_SCRIPT} http://www.bioconductor.org/packages/2.10/bioc/src/contrib/DESeq_${VERSION}.tar.gz
            ${RUN_R} CMD INSTALL DESeq_${VERSION}.tar.gz --library=${installation_path}

            cat>${installation_path}/setup.sh<<EOT
export R_LIBS=${installation_path}:${R_LIBS}
EOT

            chmod +x ${installation_path}/setup.sh

       # check package installation by checking that all index files exist
           fileList="${installation_path}/Cairo/DESCRIPTION ${installation_path}/BiocGenerics/DESCRIPTION ${installation_path}/Biobase/DESCRIPTION
  ${installation_path}/locfit/DESCRIPTION  ${installation_path}/DBI/DESCRIPTION  ${installation_path}/RSQLite/DESCRIPTION  ${installation_path}/IRanges/DESCRIPTION
   ${installation_path}/AnnotationDbi/DESCRIPTION ${installation_path}/xtable/DESCRIPTION ${installation_path}/XML/DESCRIPTION
   ${installation_path}/annotate/DESCRIPTION ${installation_path}/genefilter/DESCRIPTION ${installation_path}/RColorBrewer/DESCRIPTION
   ${installation_path}/geneplotter/DESCRIPTION ${installation_path}/DESeq/DESCRIPTION"

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
