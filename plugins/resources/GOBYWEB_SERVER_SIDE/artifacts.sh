
function install_plugin_artifacts {
    echo "Installing plugin resources"
    if [ -e ${SGE_O_WORKDIR}/artifacts-install-requests.pb ]; then
        mkdir -p ${ARTIFACT_REPOSITORY_DIR}
        REPO_MANAGER_OPTIONS="--repo-dir-quota  1000000000"
        RUN_ARTIFACT_MANAGER="java -Dlog4j.configuration=file:${GOBY_DIR}/log4j.properties \
                 -Djava.io.tmpdir=${TMPDIR} \
                 -cp ${GOBY_DIR}/serverside-dependencies.jar:${GOBY_DIR}/artifact-manager.jar \
                 org.campagnelab.gobyweb.artifacts.ArtifactManager "

        ${RUN_ARTIFACT_MANAGER} \
                --ssh-requests  ${SGE_O_WORKDIR}/artifacts-install-requests.pb \
                --repository ${ARTIFACT_REPOSITORY_DIR} ${REPO_MANAGER_OPTIONS} \
                --install

        if [ $? != 0 ]; then
             ${SGE_O_WORKDIR}/groovy ${QUEUE_WRITER} --tag ${TAG} --status ${JOB_PART_FAILED_STATUS} --description "Job failed: unable to install required software/data artifacts." --index 1 --job-type job
             exit 120
        fi
        echo "Expose environment variables for artifacts:"
        rm -f exports.sh
        ${RUN_ARTIFACT_MANAGER} \
            --repository ${ARTIFACT_REPOSITORY_DIR} ${REPO_MANAGER_OPTIONS} \
            --bash-exports --ssh-requests  ${SGE_O_WORKDIR}/artifacts-install-requests.pb \
            --output exports.sh

        cat exports.sh
        . exports.sh
        cp exports.sh ${TMPDIR}/
    fi


}