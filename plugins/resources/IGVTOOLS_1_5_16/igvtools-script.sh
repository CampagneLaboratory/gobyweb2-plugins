# COPYRIGHT_MAY_GO_HERE
# Run IGVTools with the given amount of memory
# Usage: igvtools 1g [arguments]

function igvtools {
   set -x
   set -T
   memory="$1"
   shift

   # set both minimum and max right away so failure occurs early and parallel GC is happier.
   java -Xms${memory} -Xmx${memory} -Djava.awt.headless=true -Dlog4j.debug=true -Dlog4j.configuration=file:${TMPDIR}/log4j.properties \
                      -jar ${RESOURCES_IGVTOOLS_JAR} \
                      $*
}
