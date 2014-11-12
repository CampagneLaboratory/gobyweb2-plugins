
#these wrapper functions allow to convert Mercury messages to the old file-based mechanism.
#plugin designed to exploit the messaging system offered by Mercury can include this file to have backward compatibility with GobyWeb

function trace {
   ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_START_STATUS} --description "$2" --index 0 --job-type job
}

function debug {
    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_START_STATUS} --description "$2" --index 0 --job-type job
}

function info {
    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_START_STATUS} --description "$2" --index 0 --job-type job
}

function error {
    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "$2" --index 0 --job-type job
}

function fatal {
    ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "$2" --index 0 --job-type job
}