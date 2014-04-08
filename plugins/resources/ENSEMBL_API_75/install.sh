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

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl-compara/archive/release/${VERSION}.tar.gz ensembl-compara-${VERSION}.zip
            unzip ensembl-compara-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl-variation/archive/release/${VERSION}.tar.gz ensembl-variation-${VERSION}.zip
            unzip ensembl-variation-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT}  ${ENSEMBL_ROOT_URL}/ensembl-funcgen/archive/release/${VERSION}.zip  ensembl-functgenomics-${VERSION}.zip
            unzip ensembl-functgenomics-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} ${ENSEMBL_ROOT_URL}/ensembl-tools/archive/release/${VERSION}.tar.gz  ensembl-tools-${VERSION}.zip
            unzip ensembl-tools-${VERSION}.zip

            ${RESOURCES_FETCH_URL_SCRIPT} http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz
            gzip -c -d bioperl-1.2.3.tar.gz |tar -xf -

            cd ..

            cp -r src ${installation_path}/
cat >${installation_path}/setup.sh <<EOF
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/bioperl-1.2.3
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-compara/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-variation/modules
PERL5LIB=\${PERL5LIB}:\${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/src/ensembl-functgenomics/modules
export PERL5LIB
EOF
            chmod +x  ${installation_path}/setup.sh
            return 0
            ;;

           'VEP_CACHE')

                ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`

                ${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp.ensembl.org/pub/release-${VERSION}/variation/VEP/${ORG_LOWERCASE}_vep_${VERSION}.tar.gz

                gzip -c -d  ${ORG_LOWERCASE}_vep_*.tar.gz | (cd ${installation_path} ; tar -xf -)

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