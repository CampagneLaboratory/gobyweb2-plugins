# Installation script for MuTect Homo Sapiens Data 1.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'DATA' )
            VERSION="1.1.4"
            ${RESOURCES_FETCH_URL_SCRIPT} "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf" b37_cosmic_v54_120711.vcf
            ${RESOURCES_FETCH_URL_SCRIPT} "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/dbsnp_132_b37.leftAligned.vcf.gz" dbsnp_132_b37.leftAligned.vcf.gz
            unzip dbsnp_132_b37.leftAligned.vcf.gz
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}