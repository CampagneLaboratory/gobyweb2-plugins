# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
<<<<<<< HEAD
=======
            VERSION="${PLUGIN_VERSION}"
>>>>>>> 1809ee0ac8d84a4e697cfb6f114e4e36f22bd503
            . ${RESOURCES_ARTIFACTS_ENSEMBL_API_INSTALL_DIR}/setup.sh

        # How to use : perl /scratchLocal/campagne/ARTIFACT_REPOSITORY/artifacts/ENSEMBL_API/INSTALL_DIR/70/src/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl

            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
