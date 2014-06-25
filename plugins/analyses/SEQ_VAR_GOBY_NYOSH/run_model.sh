 

export MPS_HOME=${RESOURCES_ARTIFACTS_MPS_DISTRIBUTION}
MPS_LIBS=`cat ${RESOURCES_MPS_JARS_LIST} |awk '{ORS=":"; print $1}'`
NYOSH_SUPPORT_LIBS="$RESOURCES_ARTIFACTS_MPS_SUPPORT_LIBS/*"
CLASSPATH=${MPS_LIBS}${NYOSH_SUPPORT_LIBS}:${JOB_DIR}/plugin.jar:${JOB_DIR}
MODEL=AlignmentAnalysisPlugin
NYOSH_SCRIPTNAME=AlignmentAnalysisScript
CLASSNAME=${MODEL}.${NYOSH_SCRIPTNAME}
java ${PLUGIN_NEED_DEFAULT_JVM_OPTIONS} -classpath ${CLASSPATH} -Dlog4j.configuration=file:${RESOURCES_GOBYWEB_SERVER_SIDE_LOG4J_PROPERTIES} ${CLASSNAME} "$@"
