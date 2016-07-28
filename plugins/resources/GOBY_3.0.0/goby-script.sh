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
   ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java -Xms${memory} -Xmx${memory} -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
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
   ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java ${GRID_JVM_FLAGS} -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
                                             -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function set_home {
    export GOBY_HOME=${RESOURCES_ARTIFACTS_GOBY_JAR}
}