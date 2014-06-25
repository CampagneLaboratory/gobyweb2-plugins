
# These functions work if the following variables are set in the current environment:
#  BROKER_HOSTNAME
#  BROKER_PORT
#  PLUGIN_NEED_DEFAULT_JVM_OPTIONS
#  TAG
#
# Otherwise, they do not fail, but they silently don't publish anything.
# If the functions are invoked within an oge_wrapper script or a plugin's script.sh, such variables should be available by default as they are defined in the constants.sh

function isEnabled {
    if [ -n "$BROKER_HOSTNAME" ] && [ -n "$BROKER_PORT" ] && [ -n "$PLUGIN_NEED_DEFAULT_JVM_OPTIONS" ] && [ -n "$TAG" ]
    then
        return 0
    else
        return 1
    fi
}

function publish {

    if isEnabled; then
        STARTTIME=$(($(date +%s%N)/1000000))
        if [ $# -eq 5 ]; then
            java ${PLUGIN_NEED_DEFAULT_JVM_OPTIONS} -cp ${RESOURCES_MERCURY_LIB} \
                -Dlog4j.configuration=file:${RESOURCES_MERCURY_LOG_PROPERTIES} \
                org.campagnelab.mercury.cli.JobInterface --broker-hostname ${BROKER_HOSTNAME} --broker-port ${BROKER_PORT} \
                --job-tag ${TAG} \
                --category "$1" \
                --description "$2" \
                --phase "$3" \
                --index $4 \
                --num-of-parts $5 \
                --jndi-config "${JOB_DIR}/mercury.properties"
         else
            java ${PLUGIN_NEED_DEFAULT_JVM_OPTIONS} -cp ${RESOURCES_MERCURY_LIB} \
                -Dlog4j.configuration=file:${RESOURCES_MERCURY_LOG_PROPERTIES} \
                org.campagnelab.mercury.cli.JobInterface --broker-hostname ${BROKER_HOSTNAME} --broker-port ${BROKER_PORT} \
                --job-tag ${TAG} \
                --category "$1" \
                --description "$2" \
                --phase "$3" \
                --jndi-config "${JOB_DIR}/mercury.properties"
         fi
        ENDTIME=$(($(date +%s%N)/1000000))
        echo "The message was published in $(($ENDTIME - $STARTTIME)) ms."
    fi
}

function trace {
    publish "TRACE" "$@"
}

function debug {
    publish "DEBUG" "$@"
}

function info {
    publish "INFO" "$@"
}

function error {
    publish "ERROR" "$@"
}

function fatal {
    publish "FATAL" "$@"
}