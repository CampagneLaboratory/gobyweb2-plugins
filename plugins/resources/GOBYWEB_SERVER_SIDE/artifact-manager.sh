#!/bin/bash
RUN_ARTIFACT_MANAGER="java -Dlog4j.configuration=file:${JOB_DIR}/goby/log4j.properties \
                 -Djava.io.tmpdir=${TMPDIR} \
                 -cp ${JOB_DIR}/goby/serverside-dependencies.jar:${JOB_DIR}/goby/artifact-manager.jar \
                 org.campagnelab.gobyweb.artifacts.ArtifactManager "
${RUN_ARTIFACT_MANAGER} "$@"