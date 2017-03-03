# Installation script for dlvariation version
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="1.2.3"
            ${RESOURCES_FETCH_URL_SCRIPT} http://chagall.med.cornell.edu/goby/releases/dlvariation/release-dlvariation_${VERSION}.zip
            unzip  release-dlvariation_${VERSION}.zip
            mv release-dlvariation_${VERSION}/* ${installation_path}/
            # Make version independent jars:
            cp ${installation_path}/somatic/target/somatic-${VERSION}-bin-native.jar ${installation_path}/somatic/target/somatic-bin-native.jar
            cp ${installation_path}/gpus/target/gpus-${VERSION}.jar ${installation_path}/gpus/target/gpus-bin-native.jar
            cp ${installation_path}/genotype/target/genotype-${VERSION}-bin-native.jar ${installation_path}/genotype/target/genotype-bin-native.jar

            if [ -e ${installation_path}/somatic/target/somatic-bin-native.jar ]; then
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

}