function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'SCALA_RUNTIME_2_9_2' )

            ${RESOURCES_FETCH_URL_SCRIPT} "http://www.scala-lang.org/downloads/distrib/files/scala-2.9.2.tgz"
            #${JOB_DIR}/fetch_url "http://www.scala-lang.org/downloads/distrib/files/scala-2.9.2.tgz"
            gzip -c -d scala-2.9.2.tgz	| tar -x -f -
            cp -r scala-2.9.2/* "${installation_path}/"
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}