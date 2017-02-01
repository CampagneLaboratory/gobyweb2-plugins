# Installation script for GOBY version 3.0.0
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in

        'JAR' )
            VERSION="1.2.2-SNAPSHOT"
            ${RESOURCES_FETCH_URL_SCRIPT} https://www.dropbox.com/s/x0vf5nan1fg212v/release-dlvariation_${VERSION}.zip
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