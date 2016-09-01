# Installation script for ENSEMBL_API

function get_version() {
 echo "75"
}

function plugin_install_artifact {

    id=$1
    installation_path=$2
    echo "Processing ${id}"

    VERSION=`get_version`
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

            ${RESOURCES_FETCH_URL_SCRIPT} http://archive.ubuntu.com/ubuntu/pool/universe/b/bioperl/bioperl_1.6.923.orig.tar.gz
            gzip -c -d bioperl_1.6.923.orig.tar.gz |tar -xf -
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            rm bioperl_1.6.923.orig.tar.gz

            cd ..
            cp -r src ${installation_path}/
            if [ ! $? -eq 0 ]; then
                    return 1
            fi
            ln -s ${installation_path}/src/BioPerl-1.6.923 ${installation_path}/bioperl

            if [ ! $? -eq 0 ]; then
                    return 1
            fi
cat >${installation_path}/setup.sh <<EOF
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/BioPerl-1.6.923
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-compara/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-variation/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-functgenomics/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-tools/modules
export ENSEMBL_API_VERSION=${VERSION}
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

            return 0
            ;;

           'VEP_CACHE')
#uncomment for testing with the SDK command line resource install:
#ORGANISM=homo_sapiens
#GENOME_REFERENCE_ID=1000GENOMES.37
                ORGANISM=$3
                ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`
                ${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp.ensembl.org/pub/release-${VERSION}/variation/VEP/${ORG_LOWERCASE}_vep_${VERSION}.tar.gz

                mkdir -p ${installation_path}/.vep
                gzip -c -d  ${ORG_LOWERCASE}_vep_*.tar.gz | (cd ${installation_path}/.vep ; tar -xf -)

                if [ -e ${installation_path}/.vep/${ORG_LOWERCASE} ]; then
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
#ORGANISM=homo_sapiens
#GENOME_REFERENCE_ID=1000GENOMES.37
       VERSION=`get_version`
       BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
       ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
       echo >>${out} "organism=${ORGANISM}"
       echo >>${out} "ensembl-version-number=${VERSION}"
       cat ${out}
       return 0
}

