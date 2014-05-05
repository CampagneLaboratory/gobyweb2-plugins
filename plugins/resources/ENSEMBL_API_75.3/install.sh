# Installation script for ENSEMBL_API

function plugin_install_artifact {

    id=$1
    installation_path=$2
    echo "Processing ${id}"

    VERSION="75"
    ENSEMBL_ROOT_URL="https://github.com/Ensembl"
    case ${id} in

        'INSTALL_DIR')
            mkdir src
            cd src

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl/archive/release/${VERSION}.zip ensembl-${VERSION}.zip
            unzip ensembl-${VERSION}.zip
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            mv ensembl-release-${VERSION} ensembl
            rm ensembl-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl-compara/archive/release/${VERSION}.zip ensembl-compara-${VERSION}.zip
            unzip ensembl-compara-${VERSION}.zip
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            mv ensembl-compara-release-${VERSION} ensembl-compara
            rm  ensembl-compara-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl-variation/archive/release/${VERSION}.zip ensembl-variation-${VERSION}.zip
            unzip ensembl-variation-${VERSION}.zip
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            mv ensembl-variation-release-${VERSION} ensembl-variation
            rm ensembl-variation-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT}  ${ENSEMBL_ROOT_URL}/ensembl-funcgen/archive/release/${VERSION}.zip  ensembl-functgenomics-${VERSION}.zip
            unzip ensembl-functgenomics-${VERSION}.zip
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            mv ensembl-funcgen-release-${VERSION} ensembl-functgenomics
            rm ensembl-functgenomics-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl-tools/archive/release/${VERSION}.zip ensembl-tools-${VERSION}.zip
            unzip ensembl-tools-${VERSION}.zip
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            mv ensembl-tools-release-${VERSION} ensembl-tools
            rm ensembl-tools-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz
            gzip -c -d bioperl-1.2.3.tar.gz |tar -xf -
            if [ ! $? -eq 0 ]; then
                    return 1
            fi

            mkdir ${installation_path}/bioperl

            perl Makefile.PL PREFIX=${installation_path}/bioperl INSTALLSITELIB=${installation_path}/bioperl/lib
            make
            make test
            make install

            cd ..

            cp -r src ${installation_path}/
cat >${installation_path}/setup.sh <<EOF
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/bioperl-1.2.3
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-compara/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-variation/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-funcgen/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-tools/modules
export PERL5LIB
EOF
            chmod +x  ${installation_path}/setup.sh

            ls -ltr ${installation_path}/src/

            # Check that all the pieces have been installed or fail:
            if [ ! -e ${installation_path}/src/ensembl ]; then
                    return 1
            fi
            if [ ! -e ${installation_path}/src/ensembl-tools ]; then
                    return 1
            fi
            if [ ! -e ${installation_path}/src/ensembl-variation ]; then
                    return 1
            fi
            if [ ! -e ${installation_path}/src/ensembl-functgenomics ]; then
                    return 1
            fi
            if [ ! -e ${installation_path}/src/ensembl-compara ]; then
                    return 1
            fi
            return 1 # Force fail until we know this works.
            return 0
            ;;

           'VEP_CACHE')

                ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`
                . ${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/setup.sh
 #export PERL5LIB=${RESOURCES_ARTIFACTS_VCF_TOOLS_BINARIES}/lib/perl5/site_perl:${PERL5LIB}
                perl ${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-tools/scripts/variant_effect_predictor/INSTALL.pl \
                --CACHEDIR ${installation_path} --species ${ORG_LOWERCASE}

    return 1
    exit 1
                if [ -e ${installation_path}/${ORG_LOWERCASE} ]; then
                    return 0
                else
                    return 1
                fi
                ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}



function get_attribute_values() {

    id=$1
    out=$2

       BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
       ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
       echo >>${out} "organism=${ORGANISM}"
       echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"

       echo "Printing result from ${out}:"
       cat ${out}
       return 0
}