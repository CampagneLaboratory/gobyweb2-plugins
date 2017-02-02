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
   SLF4J_CONFIG=${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}/config/goby-logback.xml
   # set both minimum and max right away so failure occurs early and parallel GC is happier.
   ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java -Xms${memory} -Xmx${memory} ${JVM_MEM_OPTIONS} \
                       -Djava.ext.dirs=${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}/ \
                       -Dlogback.configurationFile=${SLF4J_CONFIG} \
                       -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function goby3 {
   set -x
   set -T
   set_home
   mode_name="$1"
   shift
      SLF4J_CONFIG=${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}/config/goby-logback.xml
     # set both minimum and max right away so failure occurs early and parallel GC is happier.
   ${RESOURCES_ARTIFACTS_JAVA_LINUX_BINARIES}/bin/java ${GRID_JVM_FLAGS} ${JVM_MEM_OPTIONS} \
                       -Djava.ext.dirs=${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}/ \
                       -Dlogback.configurationFile=${SLF4J_CONFIG} \
                       -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}/goby.jar \
                       --mode ${mode_name} $*
}


function set_home {
    export GOBY_HOME=${RESOURCES_ARTIFACTS_GOBY3_INDELS_JAR}
}
