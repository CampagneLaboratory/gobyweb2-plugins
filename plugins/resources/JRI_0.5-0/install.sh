# Installation script for R programming language

function plugin_install_artifact {
    id=$1
    installation_path=$2

    case ${id} in
        'BINARIES' )
            VERSION="0.9-6"

            . ${RESOURCES_ARTIFACTS_R_BINARIES}/setup.sh
            RUN_R=${RESOURCES_ARTIFACTS_R_BINARIES}/bin/R

            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.r-project.org/src/contrib/rJava_${VERSION}.tar.gz rJava.tar.gz
            ${RUN_R} CMD INSTALL rJava.tar.gz --library=${installation_path}/

cat>${installation_path}/setup.sh<<EOT
export R_LIBS=${installation_path}:${R_LIBS}
EOT
        # Source the environment:
            chmod +x ${installation_path}/setup.sh
            . ${installation_path}/setup.sh

cat > script.R <<EOT
library(rJava)
.jinit()
s <- .jnew("java/io/File","mytestfile")
.jcall(s,"Z","createNewFile")
q()
EOT
            export _JAVA_OPTIONS="-Xms256m -Xmx256m"
            # Run the script, it should create mytestfile if all was installed correctly:
            ${RUN_R} CMD BATCH --no-save --no-restore script.R

            # Remove the JAVA_OPTIONS because they would even override the command line when memory options are
            # specified explicitly:
            unset _JAVA_OPTIONS

            # check that the test and installation were successful:
            if [ -e mytestfile ]; then
                 # OK we are all good
                 return 0
            else
                 # The installation failed. Unable to call Java from R
                 return 127
            fi
            return 0
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            return 99
            ;;

    esac


  exit 1

}





