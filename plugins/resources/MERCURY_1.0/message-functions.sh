
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
        echo "log enabled"
        return 0
    else
        echo "log not enabled"
        return 1
    fi
}

function publish {
    java ${PLUGIN_NEED_DEFAULT_JVM_OPTIONS} -cp ${RESOURCES_MERCURY_LIB}:${RESOURCES_MERCURY_PROPERTIES} \
        -Dlog4j.configuration=file:${RESOURCES_MERCURY_LOG_PROPERTIES} \
        org.campagnelab.mercury.cli.JobInterface --broker-hostname ${BROKER_HOSTNAME} --broker-port ${BROKER_PORT} \
        --job-tag ${TAG} \
        --text-message "$1"
}

function trace {
    echo "Publish trace message"
    if isEnabled; then
        publish "$1" "TRACE"
    fi
}

function debug {
    echo "Publish debug message"
    if isEnabled; then
        publish "$1" "DEBUG"
    fi
}

function info {
    echo "Publish info message"
    if isEnabled; then
        publish "$1" "INFO"
    fi
}

function error {
    echo "Publish error message"
    if isEnabled; then
        publish "$1" "ERROR"
    fi
}

function fatal {
    echo "Publish fatal message"
    if isEnabled; then
        publish "$1" "FATAL"
    fi
}