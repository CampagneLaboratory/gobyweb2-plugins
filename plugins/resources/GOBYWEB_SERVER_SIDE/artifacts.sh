
function install_plugin_artifacts {
    echo "Installing plugin resources"

    mkdir -p ${ARTIFACT_REPOSITORY_DIR}
    REPO_MANAGER_OPTIONS="--repo-dir-quota  1000000000"
    java -Dlog4j.configuration=file:${GOBY_DIR}/log4j.properties \
            -Djava.io.tmpdir=${TMPDIR} -jar artifact-manager.jar \
            --ssh-requests  ${SGE_O_WORKDIR}/artifacts-install-requests.pb \
            --repository ${ARTIFACT_REPOSITORY_DIR} ${REPO_MANAGER_OPTIONS} \
            --install

    if [ $? !=0 ]; then
         ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed: unable to install required software/data artifacts." --index ${CURRENT_PART} --job-type job
         exit 120
    fi
    echo "Expose environment variables for artifacts:"
    rm -f exports.sh
    java -Dlog4j.configuration=file:${GOBY_DIR}/log4j.properties \
        -Djava.io.tmpdir=${TMPDIR} -jar artifact-manager.jar \
        --repository ${ARTIFACT_REPOSITORY_DIR} ${REPO_MANAGER_OPTIONS} \
        --bash-exports --ssh-requests  ${SGE_O_WORKDIR}/artifacts-install-requests.pb \
        --output exports.sh

    cat exports.sh
    . exports.sh
    cp exports.sh ${TMPDIR}/
}