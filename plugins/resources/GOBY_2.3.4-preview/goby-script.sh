# COPYRIGHT_MAY_GO_HERE
function run_goby {
   set -x
   set -T
   memory="$1"
   shift
   mode_name="$1"
   shift
    if [ "${RJAVA_HOME:-notset}" == "notset" ]; then
       R_OPTIONS=" "
    else
       R_OPTIONS="-Djava.library.path=${RJAVA_HOME}"
    fi

   # set both minimum and max right away so failure occurs early and parallel GC is happier.
   java -Xms${memory} -Xmx${memory} ${R_OPTIONS} \
   -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
                      -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_GOBY_GOBY_JAR} \
                       --mode ${mode_name} $*
}


function goby {
   set -x
   set -T
   mode_name="$1"
   shift
    if [ "${RJAVA_HOME:-notset}" == "notset" ]; then
       R_OPTIONS=" "
    else
       R_OPTIONS="-Djava.library.path=${RJAVA_HOME}"
    fi

     # set both minimum and max right away so failure occurs early and parallel GC is happier.
   java ${GRID_JVM_FLAGS} ${R_OPTIONS} \
   -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
                       -Dgoby.configuration=file:${TMPDIR}/goby.properties \
                       -jar ${RESOURCES_GOBY_GOBY_JAR} \
                       --mode ${mode_name} $*
}
