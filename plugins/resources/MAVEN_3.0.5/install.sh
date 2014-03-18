# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'DISTRIBUTION' )
            VERSION="3.0.5"
            ${RESOURCES_FETCH_URL_SCRIPT} http://mirrors.ibiblio.org/apache/maven/maven-3/3.0.5/binaries/apache-maven-${VERSION}-bin.tar.gz
            gzip -c -d apache-maven-${VERSION}-bin.tar.gz |tar -xvf -
            cp -r apache-maven-${VERSION}/* ${installation_path}/
        # How to use : perl /scratchLocal/campagne/ARTIFACT_REPOSITORY/artifacts/ENSEMBL_API/INSTALL_DIR/75/src/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl

            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
