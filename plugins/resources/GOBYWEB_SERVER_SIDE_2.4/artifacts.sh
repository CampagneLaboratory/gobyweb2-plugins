


function install_plugin_mandatory_artifacts {
    echo "Installing plugin mandatory resources"
    install_plugin_artifacts_internal "only-mandatory"
}

function install_plugin_artifacts {
    echo "Installing plugin resources"
    install_plugin_artifacts_internal "all"
}

function install_plugin_artifacts_internal {
    if [ -e ${JOB_DIR}/artifacts-install-requests.pb ]; then
       set +xv
       . ${JOB_DIR}/constants.sh
       . ${JOB_DIR}/auto-options.sh
       local installation_type=$1
       if [ ! -z "${CURRENT_PART}" ]; then
          CURRENT_PART=1
       fi
        # Groovy is not available before artifact installation..
        # ${RESOURCES_GOBYWEB_SERVER_SIDE_QUEUE_WRITER_WRAPPER} --tag ${TAG} --status u --description "Installing artifacts on `hostname`" --index ${CURRENT_PART} --job-type job
       mkdir -p ${TMPDIR}/steplogs
       ${JOB_DIR}/artifact-manager.sh \
                --ssh-requests  ${JOB_DIR}/artifacts-install-requests.pb \
                --install --installation-type ${installation_type} --log-dir ${TMPDIR}/steplogs
       local STATUS=$?
       java -Xms40m -Xmx250m -cp ${JOB_DIR}/goby/serverside-dependencies.jar::${JOB_DIR}/goby/stepslogger.jar \
                              org.campagnelab.stepslogger.StepsLogTool \
                              --action view \
                              ${TMPDIR}/steplogs/*
       if [ ${STATUS} != 0 ]; then

            ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed: unable to install required software/data artifacts." --index ${CURRENT_PART} --job-type job
            exit 120
       fi



       echo "Expose environment variables for artifacts.."
       cd ${TMPDIR}
       rm -f exports.sh
       ${JOB_DIR}/artifact-manager.sh \
           --bash-exports --ssh-requests  ${JOB_DIR}/artifacts-install-requests.pb \
           --output ${TMPDIR}/exports-${installation_type}.sh
       if [ $? != 0 ]; then
             ${RESOURCES_GOBYWEB_SERVER_SIDE_QUEUE_WRITER_WRAPPER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed: unable to expose artifact environment." --index ${CURRENT_PART} --job-type job
             exit 121
       fi
       cat ${TMPDIR}/exports-${installation_type}.sh >> ${JOB_DIR}/exports.sh
       cat ${TMPDIR}/exports-${installation_type}.sh >> ${TMPDIR}/exports.sh
       expose_artifact_environment_variables
    fi

}

function expose_artifact_environment_variables {
        set +x
        . ${JOB_DIR}/constants.sh
        . ${JOB_DIR}/auto-options.sh
        . ${TMPDIR}/exports.sh
        set -x
}