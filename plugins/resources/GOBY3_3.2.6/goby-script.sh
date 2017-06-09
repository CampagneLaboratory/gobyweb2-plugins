JVM_MEM_OPTIONS="-XX:+UseCompressedClassPointers -XX:MaxMetaspaceSize=512m -XX:CompressedClassSpaceSize=512m"
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
   ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java -Xms${memory} -Xmx${memory} ${JVM_MEM_OPTIONS} \
                       -Dlogback.configurationFile=${RESOURCES_ARTIFACTS_GOBY3_JAR}/config/goby-logback.xml  \
                       -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY3_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function goby {
   set -x
   set -T
   set_home
   mode_name="$1"
   shift
     # set both minimum and max right away so failure occurs early and parallel GC is happier.
   ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java ${GRID_JVM_FLAGS} ${JVM_MEM_OPTIONS} \
                        -Dlogback.configurationFile=${RESOURCES_ARTIFACTS_GOBY3_JAR}/config/goby-logback.xml  \
                        -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY3_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function set_home {
    export GOBY_HOME=${RESOURCES_ARTIFACTS_GOBY3_JAR}
}

function run_goby_wrapper {
   set -x
   set -T
   set_home
   memory="$1"
   shift
   mode_name="$1"
   shift
   JVM_MEM_OPTIONS="-XX:+UseCompressedClassPointers -XX:MaxMetaspaceSize=512m -XX:CompressedClassSpaceSize=512m"
   # set both minimum and max right away so failure occurs early and parallel GC is happier.
   export JAVA_HOME=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}
   export PATH=${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin:${PATH}
   ${RESOURCES_ARTIFACTS_GOBY3_JAR}/goby ${memory} ${mode_name} $*
}