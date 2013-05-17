# Installation script for R programming language

function plugin_install_artifact {
    id=$1
    installation_path=$2

    case ${id} in
        'BINARIES' )
            VERSION="2.15.3"

            ${RESOURCES_FETCH_URL_SCRIPT} http://cran.us.r-project.org/src/base/R-2/R-${VERSION}.tar.gz

            tar -xzvf R-${VERSION}.tar.gz

            cd R-${VERSION}
            ./configure --prefix=${installation_path}
            make
            make install
            #make install-tests

cat>${installation_path}/setup.sh<<EOT
export LOCAL_LIB=${installation_path}/
export PATH=${installation_path}/bin:\${PATH}
export LD_LIBRARY_PATH=${installation_path}/lib:\${LD_LIBRARY_PATH}
export R_HOME=${installation_path}
export R_LIBS=${installation_path}/lib64/R/library
EOT

            chmod +x ${installation_path}/setup.sh
            if [ -e ${installation_path}/bin/R ]; then
                                               return 0
                                           else
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





