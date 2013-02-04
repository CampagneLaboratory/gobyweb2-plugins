# Installation script for ENSEMBL_API

function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'INSTALL_DIR')
            mkdir src
            cd src
            #wget ftp://ftp.ensembl.org/pub/ensembl-api.tar.gz
            cp ~gobyweb/url-cache/ensembl-api.tar.gz .
            gzip -c -d ensembl-api.tar.gz| tar -xvf -


            #wget http://bioperl.org/DIST/old_releases/bioperl-1.2.3.tar.gz
            cp ~gobyweb/url-cache/bioperl-1.2.3.tar.gz .
            gzip -c -d bioperl-1.2.3.tar.gz |tar -xvf -

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

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
