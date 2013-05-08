# #!/bin/sh
# . constants.sh
# . auto-options.sh


function scala {
 # arg 1: amount of memory (e.g., 3g)
 # arg 2: classpath
     export SCALA_HOME="${RESOURCES_ARTIFACTS_SCALA_SCALA_RUNTIME_2_9_2}"
     local MEM=$1
     local CLASSPATH_SCALA=$2
     shift 2
     export JAVA_HOME=/softlib/exe/x86_64/pkg/sun_jdk/6.0.2/dist
     ${SCALA_HOME}/bin/scala -J-Xmx${MEM} -cp ${CLASSPATH_SCALA} "$@"

}
