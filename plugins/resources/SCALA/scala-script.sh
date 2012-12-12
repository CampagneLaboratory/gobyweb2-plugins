# #!/bin/sh
# . constants.sh
# . auto-options.sh

# This script will either find a global Scala installation matching the SCALA plugin version (via
# the SCALA_HOME environment variable) or install the requested version in the SGE work directory.
function install {
    local filename=`basename ${RESOURCES_SCALA_TGZ_ARCHIVE} .tgz`
    local VERSION="${filename#scala-}"
    #echo version="${VERSION}"
    if [ -x "${SCALA_HOME}/bin/scala" ]; then
        # found existing SCALA_HOME installation with correct version
        return
    fi
    if [ -x "scala-${VERSION}/bin/scala" ]; then
        # Found scala of correct version in current directory
        SCALA_HOME="`pwd`/scala-${VERSION}/"
        export SCALA_HOME
        return
    fi
    if [ -z "${SCALA_HOME+xxx}" ]; then
        # SCALA_HOME is not defined, we decompress the distribution
        #echo Installing Scala resource
        gzip -c -d ${RESOURCES_SCALA_TGZ_ARCHIVE}	| tar -x -f -
        SCALA_HOME="`pwd`/scala-${VERSION}/"
        export SCALA_HOME
    fi
    echo "Scala found or installed at ${SCALA_HOME}"

}
# Install as soon as the script is sourced:
install

function scala {
 # arg 1: amount of memory (e.g., 3g)
 # arg 2: classpath
     install
     if [ "${SCALA_HOME+xxx}" ]; then
        local MEM=$1
        local CLASSPATH_SCALA=$2
        shift 2
        export JAVA_HOME=/softlib/exe/x86_64/pkg/sun_jdk/6.0.2/dist
        ${SCALA_HOME}/bin/scala -J-Xmx${MEM} -cp ${CLASSPATH_SCALA} $*
     fi
}
