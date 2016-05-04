function plugin_install_artifact {
    set -x

    id=$1
    installation_path=$2

    case ${id} in

        'BINARIES' )
            VERSION="6.1.0"
            set -x
            #${RESOURCES_FETCH_URL_SCRIPT} http://mirrors-usa.go-parts.com/gcc/releases/gcc-$VERSION/gcc-$VERSION.tar.gz gcc.tar.gz
            #${RESOURCES_FETCH_URL_SCRIPT} http://mirrors.concertpass.com/gcc/releases/gcc-$VERSION/gcc-$VERSION.tar.gz gcc.tar.gz
            ${RESOURCES_FETCH_URL_SCRIPT} file://~/usr/gcc/gcc.tar.gz gcc.tar.gz
            tar -zxvf gcc.tar.gz
            cd gcc-$VERSION
            SOURCE_DIR=$(pwd)

            ## from here: script adapted from https://gist.github.com/pmalek/5da2e947a51aa773da80

            SOURCE_DIR_NAME=${SOURCE_DIR##*/}
            cd ${SOURCE_DIR}/

            # download prerequisites for gcc's configure
            ./contrib/download_prerequisites

            # check the downloads worked
            if [ ! -d gmp ]; then
              echo "There is no gmp or mpfr or mpc directory in $(pwd). Please run download_prerequisites."
              return 127
            fi
            if [ ! -d mpfr ]; then
              echo "There is no mpfr or mpc directory in $(pwd). Please run download_prerequisites."
              return 127
            fi
            if [ ! -d gmp ]; then
              echo "There is no mpc directory in $(pwd). Please run download_prerequisites."
              return 127
            fi
<<"TRYAUTOINSTALL"
            # install gmp
            cd ${SOURCE_DIR}/gmp
            mkdir $installation_path/gmp
            ./configure --enable-shared --enable-static --prefix=$installation_path/gmp && \
            make && make check && make install
            make distclean
            echo "gmp Installed"

            # install mpfr
            cd ${SOURCE_DIR}/mpfr
            mkdir $installation_path/mpfr
            ./configure --enable-shared --enable-static --prefix=$installation_path/mpfr --with-gmp=$installation_path/gmp && \
            make && make check && make install
            make distclean
            echo "mpfr Installed"

            # install mpc
            cd ${SOURCE_DIR}/mpc
            mkdir $installation_path/mpc
            ./configure --enable-shared --enable-static --prefix=$installation_path/mpc --with-gmp=$installation_path/gmp --with-mpfr=$installation_path/mpfr && \
            make && make check && make install
            make distclean
            echo "mpc Installed"
TRYAUTOINSTALL
            echo $installation_path/lib/ >> /etc/ld.so.conf
            echo $installation_path/lib64/ >> /etc/ld.so.conf
            ldconfig

            cd ${SOURCE_DIR}/
            ulimit -s 32768 # for gcc tests
            mkdir -p $installation_path/auto-load/usr/lib
            mkdir -p gcc-build && cd gcc-build
            ${SOURCE_DIR}/configure --enable-shared --prefix=$installation_path  --enable-bootstrap --enable-languages=c,c++ --enable-libgomp --enable-threads=posix --with-fpmath=sse --disable-multilib MAKEINFO=missing
            # ${SOURCE_DIR}/configure --enable-shared --prefix=$installation_path  --enable-bootstrap --enable-languages=c,c++ --enable-libgomp --enable-threads=posix --with-gmp=$installation_path/gmp --with-mpfr=$installation_path/mpfr \
             --with-mpc=$installation_path/mpc --with-fpmath=sse --disable-multilib MAKEINFO=missing
            # force to create sym links instead of hard links (not allowed in the artifact repo mounted on docker)
            alias ln='ln -s'
            echo "alias ln='ln -s'" >> $HOME/.bashrc
            make
            make install

            # find mv .py command is due to this ldconfig error (gcc copies some .py files into /usr/local/lib64/)
            # ldconfig: /usr/local/lib/../lib64/libstdc++.so.6.0.20-gdb.py is not an ELF file - it has the wrong magic bytes at the start.
            find . -iname "*.py" -exec mv {} $installation_path/auto-load/usr/lib/ \; && ldconfig

            if [ -e ${installation_path}/bin/gcc ] && [ -e ${installation_path}/bin/g++ ]; then
               echo "gcc $VERSION Installed!"
               return 0
            else
               return 127
            fi
            return 127
            ;;

        *)  echo "Resource artifact id not recognized: "+$id
            exit 99
            ;;

    esac

    exit 1

}