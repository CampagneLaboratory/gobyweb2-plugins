# Installation script for MuTect Homo Sapiens Data 1.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'DATA' )

            ${RESOURCES_FETCH_URL_SCRIPT} "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/b37_cosmic_v54_120711.vcf" b37_cosmic_v54_120711.vcf
            cp b37_cosmic_v54_120711.vcf ${installation_path}/
            ${RESOURCES_FETCH_URL_SCRIPT} "http://www.broadinstitute.org/cancer/cga/sites/default/files/data/tools/mutect/dbsnp_132_b37.leftAligned.vcf.gz" dbsnp_132_b37.leftAligned.vcf.gz
            gunzip dbsnp_132_b37.leftAligned.vcf.gz
            cp dbsnp_132_b37.leftAligned.vcf ${installation_path}/

            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}