# Installation script
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'MODELS' )

            git clone git@bitbucket.org:campagnelaboratory/genotype-dl-models.git
            mkdir -p ${installation_path}/
            mv ./genotype-dl-models/* ${installation_path}/
            if [ -e ${installation_path}/* ]; then
                return 0
            else
               return 127
            fi
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1