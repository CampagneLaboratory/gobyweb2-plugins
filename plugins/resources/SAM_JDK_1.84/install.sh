# Installation script for STAR version 2.2.0g
function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'JAR' )
            VERSION="1.84"
            ${RESOURCES_FETCH_URL_SCRIPT} http://sourceforge.net/projects/picard/files/sam-jdk/${VERSION}/sam-${VERSION}.jar/download sam-${VERSION}.jar

            cp sam-${VERSION}.jar ${installation_path}/sam-jdk.jar
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}
