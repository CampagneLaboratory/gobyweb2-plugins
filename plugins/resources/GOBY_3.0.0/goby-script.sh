# COPYRIGHT_MAY_GO_HERE
function run_goby {
   set -x
   set -T
   set_home
   memory="$1"
   shift
   mode_name="$1"
   shift
   # set both minimum and max right away so failure occurs early and parallel GC is happier.
   java -Xms${memory} -Xmx${memory} -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
                                             -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function goby {
   set -x
   set -T
   set_home
   mode_name="$1"
   shift
     # set both minimum and max right away so failure occurs early and parallel GC is happier.
   java ${GRID_JVM_FLAGS} -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
                                             -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function set_home {
   WORKING_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
    if [[ $OSTYPE == "cygwin" ]]; then
        WORKING_DIR=`cygpath -m "${WORKING_DIR}"`
    fi

    export GOBY_HOME=${WORKING_DIR}
}