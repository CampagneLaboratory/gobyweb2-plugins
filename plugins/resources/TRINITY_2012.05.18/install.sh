function plugin_install_artifact {

    id=$1
    installation_path=$2

    case ${id} in


        'TRINITY_2012_05_18' )

            #${RESOURCES_FETCH_URL_SCRIPT} "http://hivelocity.dl.sourceforge.net/project/trinityrnaseq/trinityrnaseq_r2012-05-18.tar.gz" trinityrnaseq_r2012-05-18.tar.gz
            ${JOB_DIR}/fetch_url "http://hivelocity.dl.sourceforge.net/project/trinityrnaseq/trinityrnaseq_r2012-05-18.tar.gz" trinityrnaseq_r2012-05-18.tar.gz
            tar -zxvf trinityrnaseq_r2012-05-18.tar.gz
            cp -r trinityrnaseq_r2012-05-18/* "${installation_path}/"
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac

    exit 1

}