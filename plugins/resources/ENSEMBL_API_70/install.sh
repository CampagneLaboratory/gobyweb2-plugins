# Installation script for ENSEMBL_API

function plugin_install_artifact {

    id=$1
    installation_path=$2
    echo "Processing ${id}"
    case ${id} in

        'INSTALL_DIR')
            mkdir src
            cd src
            #wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
           # ${RESOURCES_ARTIFACTS_MAVEN_DISTRIBUTION}/bin/mvn
            ${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
            gzip -c -d ensembl-api.tar.gz| tar -xf -


            ${RESOURCES_FETCH_URL_SCRIPT} http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz
            #cp ~gobyweb/url-cache/bioperl-1.2.3.tar.gz .
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
                VERSION="70"
                . ${SGE_O_WORKDIR}/constants.sh
                . ${SGE_O_WORKDIR}/auto-options.sh
                ORG_LOWERCASE=`echo  ${ORGANISM}| tr '[:upper:]' '[:lower:]'`

                if [ -e ~/url-cache/${ORG_LOWERCASE}_vep_*.tar.gz ]; then
                    cp ~/url-cache/${ORG_LOWERCASE}_vep_*.tar.gz .
                else
                # replace with  ${RESOURCES_FETCH_URL_SCRIPT} ftp://ftp.ensembl.org/pub/release-70/variation/VEP/${ORG_LOWERCASE}_vep_\*.tar.gz

                    wget ftp://ftp.ensembl.org/pub/release-70/variation/VEP/${ORG_LOWERCASE}_vep_\*.tar.gz
                    cp ${ORG_LOWERCASE}_vep_*.tar.gz   ~/url-cache/
                fi
                mkdir -p ${installation_path}/VEP_CACHE/
                gzip -c -d  ${ORG_LOWERCASE}_vep_*.tar.gz | (cd ${installation_path}/ ; tar -xf -)

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

       . ${SGE_O_WORKDIR}/constants.sh

       BUILD_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $1}'`
       ENSEMBL_VERSION_NUMBER=`echo ${GENOME_REFERENCE_ID} | awk -F\. '{print $(NF)}'`
       echo >>${out} "organism=${ORGANISM}"
       echo >>${out} "ensembl-version-number=${ENSEMBL_VERSION_NUMBER}"

       echo "Printing result from ${out}:"
       cat ${out}
       return 0
}