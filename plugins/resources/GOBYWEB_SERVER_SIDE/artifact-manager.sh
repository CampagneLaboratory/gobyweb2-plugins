#!/bin/bash
    if [ ! -x "${JOB_DIR}" ]; then
       echo "You must define JOB_DIR before running this job."
       exit 1;
    fi
    . ${JOB_DIR}/constants.sh;
    . ${JOB_DIR}/auto-options.sh;
    mkdir -p ${ARTIFACT_REPOSITORY_DIR}

    RUN_ARTIFACT_MANAGER=" java -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                     -Djava.io.tmpdir=${TMPDIR} \
                     -cp ${JOB_DIR}/goby/serverside-dependencies.jar:${JOB_DIR}/goby/artifact-manager.jar \
                      org.campagnelab.gobyweb.artifacts.ArtifactManager \
                      --repository ${ARTIFACT_REPOSITORY_DIR} --repo-dir-quota  1000000000 "
    ${RUN_ARTIFACT_MANAGER} "$@"
