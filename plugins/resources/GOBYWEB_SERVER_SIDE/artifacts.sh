
function install_plugin_artifacts {
    echo "Installing plugin resources"
    if [ -e ${JOB_DIR}/artifacts-install-requests.pb ]; then
         . ${JOB_DIR}/constants.sh
         . ${JOB_DIR}/auto-options.sh

        if [ ! -z "${CURRENT_PART}" ]; then
          CURRENT_PART=1
        fi
        ${RESOURCES_GOBYWEB_SERVER_SIDE_QUEUE_WRITER_WRAPPER} --tag ${TAG} --status u --description "Installing artifacts on `hostname`" --index ${CURRENT_PART} --job-type job

        ${JOB_DIR}/artifact-manager.sh \
                --ssh-requests  ${JOB_DIR}/artifacts-install-requests.pb \
                --install

        if [ $? != 0 ]; then
             ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed: unable to install required software/data artifacts." --index ${CURRENT_PART} --job-type job
             exit 120
        fi

       echo "Expose environment variables for artifacts.."
       cd ${TMPDIR}
       rm -f exports.sh
       ${JOB_DIR}/artifact-manager.sh \
           --bash-exports --ssh-requests  ${JOB_DIR}/artifacts-install-requests.pb \
           --output exports.sh
       if [ $? != 0 ]; then
             ${RESOURCES_GOBYWEB_SERVER_SIDE_QUEUE_WRITER_WRAPPER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed: unable to expose artifact environment." --index ${CURRENT_PART} --job-type job
             exit 121
       fi
       expose_artifact_environment_variables
    fi

}

function expose_artifact_environment_variables {

        . ${JOB_DIR}/constants.sh
        . ${JOB_DIR}/auto-options.sh
        . ${TMPDIR}/exports.sh
}